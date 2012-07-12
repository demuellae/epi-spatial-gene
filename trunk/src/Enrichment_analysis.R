options(stringsAsFactors=F); options(width=Sys.getenv("COLUMNS"));

# VARIOUS FUNCTIONS NEEDED FOR THE ENRICHMENT ANALYSIS

# speed-improved Fisher test (a version of the R function that is reduced to the bare essentials)
hashTable = list()
fast.fisher = function (x, y = NULL, workspace = 2e+05, hybrid = FALSE, control = list(), 
    or = 1, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95, 
    simulate.p.value = FALSE, B = 2000, cache=F) 
{
    if (nrow(x)!=2 | ncol(x)!=2) stop("Incorrect input format for fast.fisher")
    if (max(is.na(x))!=0) {
    
      RVAL <- list(p.value=NA, estimate=NA, null.value=NA, alternative=NA, method="fast.fisher", data.name=NA)
      attr(RVAL, "class") <- "htest"
      return(RVAL)
    }   
    if (cache) {
      key = paste(x,collapse="_")
      cachedResult = hashTable[[key]]
      if (!is.null(cachedResult)) {
        return(cachedResult)
      }
    }
    # ---- START: cut version of fisher.test ----
    DNAME <- deparse(substitute(x))
    METHOD <- "fast.fisher"
    nr <- nrow(x)
    nc <- ncol(x)
    PVAL <- NULL
    if ((nr == 2) && (nc == 2)) {
        m <- sum(x[, 1])
        n <- sum(x[, 2])
        k <- sum(x[1, ])
        x <- x[1, 1]
        lo <- max(0, k - n)
        hi <- min(k, m)
        NVAL <- or
        names(NVAL) <- "odds ratio"
        support <- lo:hi
        logdc <- dhyper(support, m, n, k, log = TRUE)
        dnhyper <- function(ncp) {
            d <- logdc + log(ncp) * support
            d <- exp(d - max(d))
            d/sum(d)
        }
        mnhyper <- function(ncp) {
            if (ncp == 0) 
                return(lo)
            if (ncp == Inf) 
                return(hi)
            sum(support * dnhyper(ncp))
        }
        pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
            if (ncp == 1) {
                if (upper.tail) 
                  return(phyper(x - 1, m, n, k, lower.tail = FALSE))
                else return(phyper(x, m, n, k))
            }
            if (ncp == 0) {
                if (upper.tail) 
                  return(as.numeric(q <= lo))
                else return(as.numeric(q >= lo))
            }
            if (ncp == Inf) {
                if (upper.tail) 
                  return(as.numeric(q <= hi))
                else return(as.numeric(q >= hi))
            }
            d <- dnhyper(ncp)
            if (upper.tail) 
                sum(d[support >= q])
            else sum(d[support <= q])
        }
        if (is.null(PVAL)) {
            PVAL <- switch(alternative, less = pnhyper(x, or), 
                greater = pnhyper(x, or, upper.tail = TRUE), 
                two.sided = {
                  if (or == 0) 
                    as.numeric(x == lo)
                  else if (or == Inf) 
                    as.numeric(x == hi)
                  else {
                    relErr <- 1 + 10^(-7)
                    d <- dnhyper(or)
                    sum(d[d <= d[x - lo + 1] * relErr])
                  }
                })
            RVAL <- list(p.value = PVAL)
        }
        mle <- function(x) {
            if (x == lo) 
                return(0)
            if (x == hi) 
                return(Inf)
            mu <- mnhyper(1)
            if (mu > x) 
                uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
            else if (mu < x) 
                1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                  1))$root
            else 1
        }
        ESTIMATE <- mle(x)
        #names(ESTIMATE) <- "odds ratio"
        RVAL <- c(RVAL, estimate = ESTIMATE, null.value = NVAL)
    }
    RVAL <- c(RVAL, alternative = alternative, method = METHOD, data.name = DNAME)
    attr(RVAL, "class") <- "htest"
    # ---- END: cut version of fisher.test ----    
    if (cache) hashTable[[key]] <<- RVAL # write to global variable
    return(RVAL)                                                                         
}

# load BED files into a list
loadBedFiles = function(directory=character(0),filepaths=character(0),includeExtension="bed") {
  extensionPattern = list() 
  extensionPattern = paste(".",includeExtension[1],"|.",toupper(includeExtension[1]),sep="")
  if(length(includeExtension) >1){
    for (extensionInx in seq(2,length(includeExtension))){
      extensionPattern = paste(extensionPattern,".|",includeExtension[extensionInx],"|.",toupper(includeExtension[extensionInx]),sep="")
    }
  }
  filenames = character(0)
  if (length(directory)>0) {    
    #filenames = paste(directory,dir(directory,pattern=extensionPattern),sep="/")
    filenames = dir(directory,pattern=extensionPattern, full.names=TRUE)
  }
  if (length(filepaths)>0) {
    filenames = c(filenames,filepaths)
  }  
  result = list()
  for (filename in filenames) {
    print(paste("Loading:",filename))
    temp = numeric(0)
    temp = tryCatch(read.table(filename,header=T,sep="\t",comment.char="",quote=""),error=function(x) { return(numeric(0)) })
    if ((length(temp)!=0) & (sum(names(temp)=="X.chrom")) {
      # bed files starting with a header row: "#chrom  chromStart  chromEnd  name"
      names(temp)[names(temp)=="X.chrom"] = "chrom"
    } else {
      # bed files starting with a track row: "track name='E2A_GSM546517_Pre-Pro-B-Cells' description='E2A_GSM546517_Pre-Pro-B-Cells.bed' color=255,0,0"
      temp = read.table(filename,header=F,sep="\t",comment.char="track",quote="")
      colnames = "chrom\tchromStart\tchromEnd\tname\tscore\tstrand"
      if (mode(temp[,1])=="numeric") colnames = paste("bin\t",colnames,sep="")
      curNames = unlist(strsplit(colnames,"\t"))
      names(temp) = 1:ncol(temp)
      curSelection = 1:min(length(curNames),ncol(temp))
      names(temp)[curSelection] = curNames[curSelection]
    }
    if (ncol(temp)<4) temp = data.frame(temp,name=paste("region_",1:nrow(temp),sep=""))
    curLabel = gsub(extensionPattern,"",filename)
    curLabel = unlist(strsplit(curLabel,"/"))
    acceptedChromNames = unlist(strsplit("chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chrX,chrY,chrZ",","))
    accepted = temp$chrom%in%acceptedChromNames
    if (sum(accepted)/nrow(temp)<0.95) print("WARNING: a large number of regions were discarded because of incorrect chromosome names")   
    result[[curLabel[length(curLabel)]]] = temp[accepted,c("chrom","chromStart","chromEnd","name")]
  }
  return(result)
}

# export bed files in a format that is ready for bigBed conversion
writeBedFiles = function(chromatinAnnotations,directory) {
  bigBedConversionScript = paste(directory,"/0_conversionScript.txt",sep="")
  write(bigBedConversionScript,bigBedConversionScript,append=F)
  write("cd ~/Expression_analysis/Mouse_tissue_differentiation",bigBedConversionScript,append=T)  
  bigBedUrls = paste(directory,"/0_genomeBrowserUrls.txt",sep="")
  write(bigBedUrls,bigBedUrls,append=F)  
  for (curName in names(chromatinAnnotations)) {
    print(curName)
    curTab = chromatinAnnotations[[curName]]
    filename = paste(directory,"/",curName,".bed",sep="")
    write.table(curTab[order(curTab$chrom,curTab$chromStart,curTab$chromEnd),],filename,sep="\t",na=".",quote=F,row.names=F,col.names=F)
    conversionCmd = paste("~/Broad_EPP/tools/bedClip",filename,"~/Broad_EPP/mm9.chrom.sizes",gsub(".bed",".tmp",filename))
    write(conversionCmd,bigBedConversionScript,append=T)
    conversionCmd = paste("~/Broad_EPP/tools/bedToBigBed",gsub(".bed",".tmp",filename),"~/Broad_EPP/mm9.chrom.sizes",gsub(".bed",".bb",filename))
    write(conversionCmd,bigBedConversionScript,append=T)
    browserUrl = paste('http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&hgt.customText=track type=bigBed name="',curName,'" description="',curName,'" visibility=4 useScore=0 color=0,120,0 bigDataUrl=http://www.broadinstitute.org/~cbock/regions/',curName,'.bb',sep="")
    write(browserUrl,bigBedUrls,append=T)
  }
  conversionCmd = paste("rsync -avz -e ssh ",directory,"/ --include='*.bb' --exclude='*' cbock@copper.broadinstitute.org:/home/radon01/cbock/public_html/regions/",sep="")
  write(conversionCmd,bigBedConversionScript,append=T)
}

# identify region-based overlap with a list of BED files
identifyOverlap = function(curTab,annotations,includeLabels=F) {
  result = curTab
  # annotate with region overlap data
  for (annotationColName in names(annotations)) {
    print(paste("Processing:",annotationColName))
    result[,annotationColName] = rep(NA,nrow(result))
    curAnnotation = annotations[[annotationColName]]
    curAnnotation = curAnnotation[order(curAnnotation$chrom,curAnnotation$chromStart,curAnnotation$chromEnd),]
    # perform overlap analysis separately for each chromosome using the fast findInterval() function
    for (curChrom in unique(result$chrom[!is.na(result$chrom)])) {
      curSelection = curAnnotation[curAnnotation$chrom==curChrom,]
      if (nrow(curSelection)==0) next
      overlapTab = data.frame(left=findInterval(result[!is.na(result$chrom)&result$chrom==curChrom,"chromStart"],sort(curSelection[,"chromEnd"])),right=findInterval(result[!is.na(result$chrom)&result$chrom==curChrom,"chromEnd"],curSelection[,"chromStart"]))
      isOverlapping = overlapTab$left!=overlapTab$right
      isOverlapping[is.na(isOverlapping)] = F
      if (sum(isOverlapping)==0) next
      if (includeLabels) {        
        temp = apply(overlapTab[isOverlapping,],1,function(x) { paste(unique(curSelection[(x["left"]+1):x["right"],"name"]),collapse=",") })
        result[!is.na(result$chrom)&result$chrom==curChrom,][isOverlapping,annotationColName] = temp
      } else {      
        result[!is.na(result$chrom)&result$chrom==curChrom,][isOverlapping,annotationColName] = 1      
      }
    }
    overlapCount = sum(!is.na(result[,ncol(result)]))
    print(paste("Overlap: ",overlapCount," out of ",nrow(result)," regions (",round(100*overlapCount/nrow(result),2),"%)",sep=""))
  }
  return(result)
}

# helper function for translating Ensembl gene identifiers into canonical gene names
getGeneNames = function(ensemblIds,genome,geneFilesLocation=".",commaSeparatedInput=F) {
  # configuration  
  geneAnnotationDatabase = list(); 
  geneAnnotationDatabase[["hg18"]] = geneAnnotationDatabase[["hg19"]] = "hsapiens_gene_ensembl"
  geneAnnotationDatabase[["mm8"]] = geneAnnotationDatabase[["mm9"]] = "mmusculus_gene_ensembl"
  geneAnnotationDatabase[["danRer6"]] = "drerio_gene_ensembl"
  geneAnnotationColumn = list()
  geneAnnotationColumn[["hg18"]] = geneAnnotationColumn[["hg19"]] = "hgnc_symbol"
  geneAnnotationColumn[["mm8"]] = geneAnnotationColumn[["mm9"]] = "mgi_symbol"  
  geneAnnotationColumn[["danRer6"]] = "zfin_symbol"
  # retrieval of mapping table
  geneMappingTempfile = paste(geneFilesLocation,"/geneMappingTempfile_",genome,".txt",sep="")
  if (file.exists(geneMappingTempfile)) { 
    print(paste("Using existing version of",geneMappingTempfile))
    geneMapping = read.table(geneMappingTempfile,header=T,sep="\t",comment.char="",quote="",stringsAsFactors=F) 
  } else {
    print("Retrieving gene annotations from BioMart")
    library(biomaRt) # mart = useMart('ensembl'); listDatasets(mart)
    ensembl = useMart("ensembl", dataset=geneAnnotationDatabase[[genome]]); # listAttributes( mart=ensembl)
    geneMapping = getBM(attributes=c("ensembl_gene_id",geneAnnotationColumn[[genome]]), mart=ensembl)
    write.table(geneMapping,geneMappingTempfile,sep="\t",quote=F,row.names=F)
  }
  # Ensembl identifier to gene name mapping
  geneHash = geneMapping[,geneAnnotationColumn[[genome]]]
  names(geneHash) = geneMapping$ensembl_gene_id
  if (commaSeparatedInput) {
    splitted = strsplit(ensemblIds,",")
    geneNames = sapply(splitted,function(x) { temp = unique(geneHash[x]); return(paste(temp,collapse=",")) })
  } else {
    geneNames = geneHash[ensemblIds]
  }
  #geneNames = lapply(ensemblIds,function(X) { paste(sort(setdiff(unique(geneHash[X]),"")),collapse=",") } ) 
  return(geneNames)
}

# prepare functional analysis
getRegionId = function(curType,curTable,chrom="chrom",chromStart="chromStart",chromEnd="chromEnd") {
  if (curType == "meth") return(do.call("paste",curTable[,c(chrom,chromStart,chromEnd)]))
  if (curType == "expr") return(curTable$ensemblId)
  if (curType == "cpg") return(do.call("paste",curTable[,c("chrom","chromstart")]))
}

# parse a regionId and build a bed format table
regionId2regionTable = function(curType,curRegionIds, sep=" ") {
  if (curType == "meth" | curType == "expr") {
    curRegionIds = curRegionIds[!is.na(curRegionIds)]
    temp = unlist(strsplit(curRegionIds,sep))
    temp = matrix(temp,ncol=3,nrow=length(temp)/3,byrow=T)
    temp = data.frame(chrom=temp[,1],chromStart=as.numeric(temp[,2]),chromEnd=as.numeric(temp[,3]),name=curRegionIds)
    temp = temp[order(temp$chrom,temp$chromStart),]
    row.names(temp) = 1:nrow(temp)
    return(temp)
  } else if (curType == "cpg") {
    temp = unlist(strsplit(curRegionIds," "))
    temp = matrix(temp,ncol=2,nrow=length(temp)/2,byrow=T)
    temp = data.frame(chrom=temp[,1],chromStart=as.numeric(temp[,2]),chromEnd=as.numeric(temp[,2])+2,name=curRegionIds)
    temp = temp[order(temp$chrom,temp$chromStart),]
    row.names(temp) = 1:nrow(temp)
    return(temp)    
  }
}

# perform split on comma and condense to unique list of identifiers
deconvoluteGeneSet = function(geneSet) {
  result = unique(unlist(strsplit(geneSet,",")))
  return(result[!is.na(result)])
}

# obtain a representative genomic region based on a set of identifier-linked regions
obtainRegionsById = function(curTable,selectionMode="middle") {
  curTable = curTable[order(curTable$chrom,curTable$chromStart),]
  splitted = split(curTable,curTable[,"name"])
  if (selectionMode=="middle") combined = lapply(splitted, function(X) { return(X[round(1+(nrow(X)-1)/2),,drop=F]) })
  if (selectionMode=="random") combined = lapply(splitted, function(X) { return(X[sample(1:nrow(X),1),,drop=F]) })  
  if (selectionMode=="range") combined = lapply(splitted, function(X) { return(data.frame(chrom=X$chrom[1],chromStart=min(X$chromStart),chromEnd=max(X$chromEnd),name=X$name[1])) })
  result = do.call(rbind, lapply(combined, data.frame))
  return(result)
}                 

# high-speed lookup in a named list using a temporary hash
library(hash)
lookupGenes2 = function(curId,curLookupTable) {
  n = 5000
  temp = hash(keys=names(curLookupTable[1:n]),values=curLookupTable[1:n])
  result = temp[curId]
  return(result)
}

# high-speed lookup in a named list (only one gene per item, no support for comma-separated values)
lookupGenes = function(curRegionId,curLookupTable) {
  batchSize = 10000
  result = character(0)
  for (i in 0:(length(curRegionId)/batchSize)) {
    lower = (i*batchSize+1)
    higher = min((i+1)*batchSize,length(curRegionId))
    if (lower <= higher) result = c(result,curLookupTable[curRegionId[lower:higher]])
  }
  return(result)
}



# PREPARATIONS FOR THE ENRICHMENT ANALYSIS (not sure if needed)
 
# conversion of region IDs into Ensembl identifiers
#config_file <- ""
config_dir  <- "."
annotation_dir  <- "../../data/bedfile"
bed_file <- ""

# most probably this part creates a mapping function for ensembl genes
if(FALSE){
geneTable_promoterNarrow = loadBedFiles(filepaths=c("analysis_config/Genes_Ensembl_promoter.bed")) # -5kb to +1kb
ensemblId2region_promoterNarrow = obtainRegionsById(geneTable_promoterNarrow[[1]],selectionMode="middle")
ensemblId2region_promoterNarrow$regionId = getRegionId("meth",ensemblId2region_promoterNarrow)
geneTable_promoterWide = loadBedFiles(filepaths=c("analysis_config/mm9.52.Ensembl_TSS_20kb.bed")) # -10kb to +10kb
ensemblId2region_promoterWide = obtainRegionsById(geneTable_promoterWide[[1]],selectionMode="middle")
ensemblId2region_promoterWide$regionId = getRegionId("meth",ensemblId2region_promoterWide)
#geneTable_locus = loadBedFiles(filepaths=c("analysis_config/mm9.51.Ensembl_locus.bed")) # txstart-50kb to txend+50kb
geneTable_locus = loadBedFiles(filepaths=c("analysis_config/Genes_Ensembl_locus.bed")) # txstart-50kb to txend+50kb
ensemblId2region_locus = obtainRegionsById(geneTable_locus[[1]],selectionMode="range")
ensemblId2region_locus$regionId = getRegionId("meth",ensemblId2region_locus)
regionId2ensemblId_promoterNarrow = list()
regionId2ensemblId_promoterWide = list()
regionId2ensemblId_locus = list()
for (curType in c("meth","cpg")) {
  if (curType=="meth") curTab = meanTab[[curType]][,c("chrom","chromStart","chromEnd")]
  if (curType=="cpg") curTab = data.frame(chrom=meanTab[[curType]]$chrom,chromStart=meanTab[[curType]]$chromstart,chromEnd=meanTab[[curType]]$chromstart+1)
  print(paste(curType,"promoter-narrow"))
  temp = identifyOverlap(curTab,geneTable_promoterNarrow,includeLabels=T)
  regionId2ensemblId_promoterNarrow[[curType]] = temp[,4]
  names(regionId2ensemblId_promoterNarrow[[curType]]) = getRegionId(curType,meanTab[[curType]])
  print(paste(curType,"promoter-wide"))
  temp = identifyOverlap(curTab,geneTable_promoterWide,includeLabels=T)
  regionId2ensemblId_promoterWide[[curType]] = temp[,4]
  names(regionId2ensemblId_promoterWide[[curType]]) = getRegionId(curType,meanTab[[curType]])
  print(paste(curType,"locus"))
  temp = identifyOverlap(curTab,geneTable_locus,includeLabels=T)
  regionId2ensemblId_locus[[curType]] = temp[,4]
  names(regionId2ensemblId_locus[[curType]]) = getRegionId(curType,meanTab[[curType]])
}  
ensemblIds = c(unlist(regionId2ensemblId_promoterNarrow),unlist(regionId2ensemblId_promoterWide),unlist(regionId2ensemblId_locus))
ensemblIds = ensemblIds[!is.na(ensemblIds)]
ensemblIds = unlist(strsplit(ensemblIds,","))
ensemblIds = unique(ensemblIds)
ensemblId2geneName = getGeneNames(ensemblIds,"mm9","analysis_config")

# load data for chromatin analysis
chromatinAnnotations = loadBedFiles("analysis_config/regions_mm9")
writeBedFiles(chromatinAnnotations,"analysis_config/regions_mm9/bigBed")
}



# KEY FUNCTIONS OF THE ENRICHMENT ANALYSIS

# perform chromatin enrichment analysis based on preloaded curAnnotation
performChromatinAnalysis = function(cases,background,curAnnotation,listGenes=F,pvalThreshold=0.001) {
  # prepare table with cases, controls and annotation
  cases = cases[!is.na(cases)]
  background = union(background[!is.na(background)],cases)
  annotatedTable = curAnnotation[curAnnotation$name%in%background,]
  classVal = annotatedTable$name%in%cases
  # calculate enrichment in cases relative to controls
  library(qvalue)
  overlap = character(0)
  pValues = numeric(0)
  oddsRatios = numeric(0)
  supportVals = numeric(0)
  dataCols = names(annotatedTable)[5:ncol(annotatedTable)]
  for (curCol in dataCols) {
    m1 = sum(classVal & !is.na(annotatedTable[,curCol]))
    m2 = sum(classVal) - m1 # == length(intersect(cases,setdiff(background,geneSet)))
    m3 = sum(!is.na(annotatedTable[,curCol])) - m1 # == length(intersect(setdiff(background,cases),geneSet))
    m4 = nrow(annotatedTable) - m1 - m2 - m3 # == length(intersect(setdiff(background,cases),setdiff(background,geneSet)))
    curMatrix = matrix(c(m1,m2,m3,m4),nrow=2,ncol=2,dimnames=list(c("cases","controls"),c("gene set","other")))    
    fisherTestResults = fast.fisher(curMatrix,alternative="greater")
    pValues = c(pValues,fisherTestResults[["p.value"]])
    oddsRatios = c(oddsRatios,fisherTestResults[["estimate"]])
    supportVals = c(supportVals,m1)
    if (listGenes & fisherTestResults[["p.value"]] <= pvalThreshold ) {
      temp = annotatedTable[classVal & !is.na(annotatedTable[,curCol]),"name"]
      #regionId2ensemblId_locus[["meth"]][temp]
      overlap = c(overlap,paste(temp,collapse=", "))
    } else {
      overlap = c(overlap,"")
    }    
  }
  logOdds = log(oddsRatios,2)
  logOdds = ifelse(logOdds< -10,-10,logOdds)
  logOdds = ifelse(logOdds>10,10,logOdds)
  qvalueObject = try(qvalue(signif(pValues,8)))
  if ("qvalues" %in% names(qvalueObject)) { qValues = qvalueObject$qvalues } else { qValues = pmin(pValues*length(pValues),1) }
  rank_pValues = rank(pValues,ties.method="random")
  rank_logOdds = rank(-logOdds,ties.method="random")
  rank_supportVals = rank(-supportVals,ties.method="random")
  maxRank = pmax(rank_pValues,rank_logOdds,rank_supportVals)
  meanRank = rank(rank_pValues+rank_logOdds+rank_supportVals,ties.method="random")
  if (listGenes) {
    result = data.frame(label=dataCols,logOdds=logOdds,support=supportVals,pvalue=pValues,qvalue=qValues,meanRank,maxRank,overlap=overlap)
  } else {
    result = data.frame(label=dataCols,logOdds=logOdds,support=supportVals,pvalue=pValues,qvalue=qValues,meanRank,maxRank)#,rank_pValues,rank_logOdds,rank_supportVals)
  }
  return(result)
}

# perform functional enrichment analysis based on MSigDB gene sets
performTranscriptomeAnalysis = function(cases,background,geneSets,inputFormat="ensemblId") {
  cases = unlist(strsplit(cases,","))
  background = unlist(strsplit(background,","))
  if (inputFormat == "ensemblId") {
    cases = ensemblId2geneName[cases]
    background = ensemblId2geneName[background]
  }  
  library(qvalue)
  cases = unique(toupper(cases)); background = unique(toupper(background))
  background = union(cases,background)
  pValues = numeric(0)
  oddsRatios = numeric(0)
  supportVals = numeric(0)
  for (i in 1:length(geneSets)) {
    geneSetName = names(geneSets)[i]
    geneSet = intersect(toupper(geneSets[[geneSetName]]),background)
    m1 = length(intersect(geneSet,cases))
    m2 = length(cases) - m1 # == length(intersect(cases,setdiff(background,geneSet)))
    m3 = length(geneSet) - m1 # == length(intersect(setdiff(background,cases),geneSet))
    m4 = length(background) - m1 - m2 - m3 # == length(intersect(setdiff(background,cases),setdiff(background,geneSet)))
    curMatrix = matrix(c(m1,m2,m3,m4),nrow=2,ncol=2,dimnames=list(c("cases","controls"),c("gene set","other")))    
    fisherTestResults = fast.fisher(curMatrix,alternative="greater")
    pValues = c(pValues,fisherTestResults[["p.value"]])
    oddsRatios = c(oddsRatios,fisherTestResults[["estimate"]])
    supportVals = c(supportVals,m1)
  }
  logOdds = log(oddsRatios,2)
  logOdds = ifelse(logOdds< -10,-10,logOdds)
  logOdds = ifelse(logOdds>10,10,logOdds)
  qvalueObject = try(qvalue(signif(pValues,8)))
  if ("qvalues" %in% names(qvalueObject)) { qValues = qvalueObject$qvalues } else { qValues = pmin(pValues*length(pValues),1) }
  rank_pValues = rank(pValues,ties.method="random")
  rank_logOdds = rank(-logOdds,ties.method="random")
  rank_supportVals = rank(-supportVals,ties.method="random")
  maxRank = pmax(rank_pValues,rank_logOdds,rank_supportVals)
  meanRank = rank(rank_pValues+rank_logOdds+rank_supportVals,ties.method="random")
  result = data.frame(geneSet=names(geneSets),logOdds=logOdds,support=supportVals,pvalue=pValues,qvalue=qValues,meanRank,maxRank)#,rank_pValues,rank_logOdds,rank_supportVals)
  return(result)
}    

# this function is a wrapper around the normal Wilcoxon test that suppresses warnings and returns NA instead of raising an error if a problem occurs
wilcoxonTestFailsafe = function(...) { 
  options(warn=-1) # temporarily disable warnings related to inexact p-values in the presence of ties
  obj = try(wilcox.test(...), silent=TRUE) 
  options(warn=1)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

# perform EpiGRAPH analysis using a precalculated annotation table
performGenomeAnalysis = function(cases,background,curAnnotation) {
  library(qvalue)
  # prepare table with cases, controls and annotations
  cases = cases[!is.na(cases)]
  background = union(background[!is.na(background)],cases)
  chromCol = names(curAnnotation)[grep("_chrom$",names(curAnnotation))]
  chromStartCol = names(curAnnotation)[grep("_chromStart$",names(curAnnotation))]
  chromEndCol = names(curAnnotation)[grep("_chromEnd$",names(curAnnotation))]
  curAnnotation$UseAtt_regionId = do.call("paste",curAnnotation[,c(chromCol,chromStartCol,chromEndCol)])
  annotatedTable = curAnnotation[curAnnotation$UseAtt_regionId%in%background,]
  classVal = annotatedTable$UseAtt_regionId%in%cases
  dataCols = names(annotatedTable)[setdiff(1:ncol(annotatedTable),grep("UseAtt",names(annotatedTable)))]
  dataCols = setdiff(dataCols,grep("EpiAnd_DnaMet",dataCols))
  curResult = data.frame()
  for (curCol in dataCols) {    
    minSupport = min(sum(!is.na(annotatedTable[classVal,curCol])),sum(!is.na(annotatedTable[!classVal,curCol])))
    minSupportRatio = minSupport/min(sum(classVal),sum(!classVal))
    pval = wilcoxonTestFailsafe(annotatedTable[classVal,curCol],annotatedTable[!classVal,curCol])
    if (is.na(pval) | minSupportRatio<0.8) next
    caseMean = mean(annotatedTable[classVal,curCol],na.rm=T)
    refMean = mean(annotatedTable[!classVal,curCol],na.rm=T)
    curRecord = data.frame(curCol,caseMean,refMean,minSupport,pval)
    names(curRecord) = c("attributeName","meanValueCases","meanValueControls","minSupport","pValues")
    curResult = rbind(curResult,curRecord)
  }
  if (nrow(curResult)==0) {
    return(NULL)
  } else {
    qvalueObject = try(qvalue(signif(curResult$pValues,8)))
    if ("qvalues" %in% names(qvalueObject)) { qValues = qvalueObject$qvalues } else { qValues = pmin(curResult$pValues*length(curResult$pValues),1) }
    absDiff = abs(curResult$meanValueCases - curResult$meanValueControls)
    curOffset = median(abs(c(curResult$meanValueCases,curResult$meanValueControls)))/1000
    relDiffLog = log2((absDiff+curOffset)/(abs(curResult$meanValueControls)+curOffset))
    rank_pValues = rank(curResult$pValues,ties.method="random")
    rank_relDiffLog = rank(-relDiffLog,ties.method="random")  
    rank_supportVals = rank(-curResult$minSupport,ties.method="random")
    maxRank = pmax(rank_pValues,rank_relDiffLog,rank_supportVals)
    meanRank = rank(rank_pValues+rank_relDiffLog+rank_supportVals,ties.method="random")
    result = data.frame(curResult,qValues,absDiff,relDiffLog,meanRank,maxRank)
    return(result)
  }  
}

# import data from MSigDB
loadGeneSets = function(filename,minGenesPerSet=5,maxGenesPerSet=Inf) {
  geneSetsRaw = scan(filename,what=character(0),sep="\n",quote="",strip.white=T,blank.lines.skip=T,multi.line=F)
  temp = strsplit(geneSetsRaw,"\t| /// | // | ///| //| /")
  geneSets = lapply(temp,function(x) { if (length(x) <= 2) result = NA else result = x[3:length(x)]; return(sort(unique(setdiff(result,"")))) })
  names(geneSets) = lapply(temp,function(x) { return(x[1]) })
  geneSetSizes = sapply(geneSets,length)
  return(geneSets[geneSetSizes>=minGenesPerSet & geneSetSizes<=maxGenesPerSet])
}

# import data from EpiGRAPH
loadEpigraphAnnotations = function(curType,regionType="") {
  curAnnotations = numeric(0)
  for (i in 1:5) {
    filename = NA
    if (curType=="meth" & (regionType=="1kb_tiling" | regionType=="Tiling_regions_1kb")) filename = "EpiGRAPH_annotations_meth_tiling_XX.txt" 
    if (curType=="meth" & regionType=="Promoter_ensembl") filename = "EpiGRAPH_annotations_meth_promoter_XX.txt"
    if (curType=="expr") filename = "EpiGRAPH_annotations_expr_XX.txt"    
    filename = paste("analysis_config/",gsub("XX",i,filename),sep="")
    if (!file.exists(filename)) break
    epigraph = read.table(filename,header=T,sep="\t",comment.char="",quote="") 
    if (length(curAnnotations)==0) {
      curAnnotations = epigraph
    } else {
      sharedCols = intersect(names(curAnnotations),names(epigraph))
      curAnnotations = rbind(curAnnotations[,sharedCols],epigraph[,sharedCols])
    }
  }
  curAnnotations = curAnnotations[,setdiff(names(curAnnotations),names(curAnnotations)[grep("EpiAnd_DnaMet_Meth",names(curAnnotations))])]  
  return(curAnnotations)
}

# prepare enrichment analysis
#chromatinAnnotations = loadBedFiles(annotation_dir)
chromatinDir = c("analysis_config/regions_mm9",
	       	"dataset/mm9/histone/encode/caltech","dataset/mm9/histone/encode/licr", "dataset/mm9/histone/encode/psu", "dataset/mm9/histone/encode/sydh",
	       	"dataset/mm9/tfbs/encode/caltech","dataset/mm9/tfbs/encode/licr", "dataset/mm9/tfbs/encode/psu", "dataset/mm9/tfbs/encode/sydh",
		"dataset/mm9/uw")
chromatinAnnotations = loadBedFiles(chromatinDir, includeExtension=c("bed","narrowPeak","broadPeak"))
regionType="1kb_tiling"
#geneSets = loadGeneSets(paste(annotation_dir, "/msigdb.v3.0.symbols.gmt", sep=""))
geneSets = loadGeneSets("analysis_config/regions_mm9/msigdb.v3.0.symbols.gmt")
epigraphAnnotations = list()
epigraphAnnotations[["meth"]] = loadEpigraphAnnotations("meth",regionType)
epigraphAnnotations[["expr"]] = loadEpigraphAnnotations("expr",regionType)
annotatedTableChromatin = list()
for (curType in c("meth","expr")) {  # !!! ,"cpg")) {
  print(curType)
  curBackground = character(0)
  for (curName in names(regionItemSets[[curType]])) curBackground = unique(c(curBackground,regionItemSets[[curType]][[curName]]))
  curBackground = curBackground[!is.na(curBackground)]
  if (length(curBackground) > 0) {
    annotatedTableChromatin[[curType]] = identifyOverlap(regionId2regionTable(curType,curBackground),chromatinAnnotations,includeLabels=F)
  } else {
    annotatedTableChromatin[[curType]] = character(0)
  }
} 

save.image(paste("session_enrichmentAnalysis.",analysis,regionType,".bin",sep=""))



performBicluster = FALSE
sepSymbol=":"
biclusterDir="./biclusters/"
if (!file.exists(biclusterDir)) dir.create(biclusterDir)
# BICLUSTERING using BCBimax with parameters chosen by grid search giving best GO enrichment
if(exists("performBicluster")){
  if(performBicluster){
    library(biclust)
    library(biomaRt)

    #reading the file matrix.csv is comma separated with 3 lines removed from  matrix_1224667147767.txt        
    intM = read.table(file="eurexpress/matrix.csv", header=TRUE, row.names=2, check.names=FALSE, sep=",")
    # removing the first column as it have assay information. intM is now numeric 
    intM = intM[,-1:0]
    numCluster = 25
    intM.binary <- as.matrix(intM)
    intM.binary[,] <- 0
    intM.binary[intM > 0] <- 1
    intMClust = biclust(intM.binary, method=BCBimax(), minr = ceiling(5000/(1.5 * numCluster)), minc=ceiling(800/(1.5 * numCluster)) , number=numCluster)


    ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
    genetab = getBM(attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "entrezgene", "mgi_symbol" ), mart=ensembl)
    genetabNAN = genetab[!is.na(genetab$mgi_symbol),]
    #background.table <- genetabNAN[genetabNAN$entrezgene %in% entrezName(rownames(intM)),] # through entrezName
    sepSymbol=":"
    geneUniverse = rownames(intM)
    background.table =  genetabNAN[genetabNAN$mgi_symbol %in% geneUniverse,] 
    background = data.frame(regionId = paste("chr",background.table$chromosome_name,sepSymbol,background.table$start_position,sepSymbol,background.table$end_position,sep=""), ensemblId = background.table$ensembl_gene_id)
    filename  <- paste(biclusterDir,"background.txt",sep="")
    write.table(file=filename,x=background, quote=F, row.names=F)
    for(clust in seq(1,numCluster)){
      geneSample  <- rownames(intM)[intMClust@RowxNumber[,clust]]
      bimax.table =  genetabNAN[genetabNAN$mgi_symbol %in% geneSample,] 
      bimax = data.frame(regionId = paste("chr",bimax.table$chromosome_name,sepSymbol,bimax.table$start_position,sepSymbol,bimax.table$end_position,sep=""), ensemblId = bimax.table$ensembl_gene_id)
      filename  <- paste(biclusterDir,"bimax",numCluster,"C.",clust,".txt",sep="")
      write.table(file=filename,x=bimax, quote=F, row.names=F)

    }
  }

}

# MANUALLY PERFORMING AN ENRICHMENT ANALYSIS
# configure analysis based on region IDs
a = "manual_enrichment_analyses Jul 10th"
outputDir = "../result/12jul/"
if (!file.exists(outputDir)) dir.create(outputDir)
outputDir = "../result/12jul/manual_enrichment_analyses/"
if (!file.exists(outputDir)) dir.create(outputDir)
regionIdList = list()
ensemblIdList = list()
#input gene set
#inputTable = read.table("analysis_config/Manual_enrichment_Jul10__Gain_methylation_composite_no_oldsmalls_singlemiss.txt",header=T,sep="\t",comment.char="",quote="") 
#inputTable = read.table(paste(biclsuterDir,"/bimax25C1.txt", sep="", header=T)
filename = paste(biclusterDir, "biclusters/bimax25C.1.txt", sep="")
inputTable = read.table(filename, sep="", header=T)
regionIdList[["Manual_enrichment_Jul10__Gain_methylation_composite_no_oldsmalls_singlemiss"]] = inputTable[,"regionId"]
ensemblIdList[["Manual_enrichment_Jul10__Gain_methylation_composite_no_oldsmalls_singlemiss"]] = inputTable[,"ensemblId"]

# load regionId background
#inputTable = read.table("analysis_config/Manual_jointBackground_regionIds.txt",header=T,sep="\t",comment.char="",quote="")
filename = paste(biclusterDir, "biclusters/background.txt", sep="")
inputTable = read.table(filename, sep="", header=T)
regionIdList[["background"]] = inputTable[,"regionId"]
regionIdAll = unname(unlist(regionIdList))

# retrieve ensemblId background
print("Retrieving gene annotations from BioMart")
library(biomaRt)
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneMapping = getBM(attributes=c("ensembl_gene_id","mgi_symbol"), mart=ensembl)
ensemblIdList[["background"]] = toupper(unique(geneMapping[,"ensembl_gene_id"]))

# load data for chromatin analysis
annotatedTableChromatin = list()
annotatedTableChromatin[[curType]] = identifyOverlap(regionId2regionTable(curType,regionIdAll,sep=sepSymbol),chromatinAnnotations,includeLabels=F)

# perform enrichment analysis
maxDigits = 3
for (i in 1:length(names(regionIdList))) {
  for (j in 1:length(names(regionIdList))) {
    if (i >=j ) next
    first = names(regionIdList)[i]
    second = names(regionIdList)[j]
    curLabel = paste(first,"_vs_",second,sep="")    
    if (first %in% names(regionIdList)) {
      print(paste("Chromatin analysis for:",curLabel))
      # prepare region-based analysis
      casesRegionId = sort(unique(regionIdList[[first]]))
      backgroundRegionId = sort(unique(union(regionIdList[[first]],regionIdList[[second]])))
      # chromatin analysis
      temp = performChromatinAnalysis(casesRegionId,backgroundRegionId,annotatedTableChromatin[[curType]],listGenes=T)      
      filename = paste(outputDir,"/EnrichedChromatin_",curLabel,".txt",sep="")
      write.table(format(temp,trim=T,digits=maxDigits),filename,sep="\t",quote=F,row.names=F)            
    }
    if(FALSE){
    if (first %in% names(ensemblIdList)) {
      print(paste("Gene set analysis for:",curLabel))
      # prepare gene-based analysis
      casesEnsemblId = sort(unique(deconvoluteGeneSet(ensemblIdList[[first]])))
      backgroundEnsemblId = sort(unique(union(deconvoluteGeneSet(ensemblIdList[[first]]),deconvoluteGeneSet(ensemblIdList[[second]]))))
      curLabel = paste(first,"_vs_",second,sep="")    
      print(curLabel)
      # chromatin analysis
      temp = performTranscriptomeAnalysis(casesEnsemblId,backgroundEnsemblId,geneSets)
      filename = paste(outputDir,"/EnrichedGeneSets_",curLabel,".txt",sep="")
      write.table(format(temp,trim=T,digits=maxDigits),filename,sep="\t",quote=F,row.names=F)            
    }
    }
  }
}

