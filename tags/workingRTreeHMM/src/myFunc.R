goSignificantCluster  <- function(cluster, geneMat, entrezName, pvalueCutoff = 0.5 ){ 
  library(biclust)
  library(GOstats)
  library(org.Mm.eg.db)
  geneUniverse  <- entrezName[rownames(intM)]
  geneUniverse.entrez  <- unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1])))

  sigBicluter <- 0
  Allcluster  <- NULL 
  sigGO  <- NULL
  GOcategories  <- c("MF")
  #for(i in seq(1,20)){
  for(i in seq(1,cluster@Number)){
    print(i)
    geneSample  <- entrezName[rownames(geneMat)[cluster@RowxNumber[,i]]]
    geneSample.entrez  <- 	unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1])))
    for(goCategory in GOcategories){
      params <- new("GOHyperGParams", geneIds = as.vector(geneSample.entrez) , universeGeneIds = as.vector(geneUniverse.entrez), annotation="org.Mm.eg.db", ontology = goCategory, pvalueCutoff = pvalueCutoff, conditional = TRUE, testDirection = "over")
      hgOver <- hyperGTest(params)
      Allcluster  <- c(Allcluster, hgOver)	
      df  <-  summary(hgOver)
      if (dim(df)[1] > 0 ){
	sigBicluter <- sigBicluter + 1
      }
      #print(sigBicluter)
      sigGO  <- c(sigGO, sigCategories(hgOver))
    }

  }
  sigBicluter  <- (sigBicluter/(cluster@Number))
  sigGO  <- unique(sigGO)
  out  <- NULL
  out$sigBicluter  <- sigBicluter
  out$sigGO  <- length(sigGO)
  out$Allcluster  <-  Allcluster
  return(out)
}
goClust <- function(geneSample.entrez, geneMat=geneMat, cluster=cluster, entrezName=entrezName, pvalueCutoff, geneUniverse.entrez=geneUniverse.entrez )
{
  library(GOstats)
  library(org.Mm.eg.db)
  GOcategories  <- c("BP","MF", "CC")
  #goCatSigCnt <- 0
  Allcluster <- NULL
  Allcluster <- NULL
  for(goCategory in GOcategories){
    params <- new("GOHyperGParams", geneIds = as.vector(geneSample.entrez) , universeGeneIds = as.vector(geneUniverse.entrez), annotation="org.Mm.eg.db", ontology = goCategory, pvalueCutoff = pvalueCutoff, conditional = TRUE, testDirection = "over")
    hgOver <- hyperGTest(params)
    Allcluster  <- c(Allcluster, hgOver)	
    #df  <-  summary(hgOver)
    #if (dim(df)[1] > 0 )
      #goCatSigCnt  <- goCatSigCnt + 1
  }
  return(Allcluster)
}

goClustMF <- function(geneSample.entrez, geneMat=geneMat, cluster=cluster, entrezName=entrezName, pvalueCutoff, geneUniverse.entrez=geneUniverse.entrez )
{
  library(GOstats)
  library(org.Mm.eg.db)
  library(GO.db)
  #goCatSigCnt <- 0
  Allcluster <- NULL
  Allcluster <- NULL
  print(geneSample.entrez)
    params <- new("GOHyperGParams", geneIds = as.vector(geneSample.entrez) , universeGeneIds = as.vector(geneUniverse.entrez), annotation="org.Mm.eg.db", ontology = "MF", pvalueCutoff = pvalueCutoff, conditional = TRUE, testDirection = "over")
    hgOver <- hyperGTest(params)
    #df  <-  summary(hgOver)
    #if (dim(df)[1] > 0 )
      #goCatSigCnt  <- goCatSigCnt + 1
  return(hgOver)
}
goSignificantClusterMulticore  <- function(cluster, geneMat, entrezName, pvalueCutoff = 0.5){ 
  library(biclust)
  library(multicore)
options(cores = 25)
  geneUniverse  <- entrezName[rownames(intM)]
  geneUniverse.entrez  <- unique(unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1]))))

  sigBicluter <- 0
  Allcluster  <- NULL 
  sigGO  <- NULL
  geneSample.entrez = list()
  for(i in seq(1,cluster@Number)){
  #for(i in seq(1,20)){
    geneSample  <- entrezName[rownames(geneMat)[cluster@RowxNumber[,i]]]
    #geneSample.entrez =  c(geneSample.entrez,unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1]))))
    geneSample.entrez[i] = list(unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1]))))
  }
  #clustseq <- seq(1,cluster@Number)
  #detach(package:GO.db)
  #detach(package:org.Mm.eg.db)
  #detach(package:GOstats)
  #detach(package:Category)
  #detach(package:RSQLite)
  #detach(package:AnnotationDbi)
  print(geneSample.entrez)
  #Allcluster  <- mclapply(geneSample.entrez, goClustMF, geneMat=geneMat, cluster=cluster, pvalueCutoff=pvalueCutoff, geneUniverse.entrez=geneUniverse.entrez, mc.cores=32, mc.preschedule=F) 
  #Allcluster  <- lapply(geneSample.entrez, goClustMF, geneMat=geneMat, cluster=cluster, pvalueCutoff=pvalueCutoff, geneUniverse.entrez=geneUniverse.entrez) 
  Allcluster  <-lapply(geneSample.entrez, goClust, geneMat=geneMat, cluster=cluster, pvalueCutoff=pvalueCutoff, geneUniverse.entrez=geneUniverse.entrez) 
  #sigBicluter  <- (sigBicluter/(cluster@Number))
  #sigGO  <- unique(sigGO)
  #out$sigBicluter  <- sigBicluter
  #out$sigGO  <- length(sigGO)
  #out$Allcluster  <-  Allcluster
  return(Allcluster)
}

numSigCluster <- function(out, Number, pThershold=1e-2 )
{
  sigBicluster  <- 0
  for(i in seq(1, Number)){
    pVal  <- pvalues(out$Allcluster[[i]])
    #print(pVal[1])	
    #print(sum(pVal  < pThershold))
    if (sum(pVal  < pThershold) > 0)
      sigBicluster  <- sigBicluster + 1
    #print(sigBicluster)
  }
  return(sigBicluster)
}

treeApply <- function(x, func, post=NULL, pre=NULL, ...) {

  require(XML)
  ans <- NULL

  value <- func(x)

  if(length(value))
    ans <- list(value=value)

  # If there are any component, do a recursive apply on those also.
  # If the result is non-null
  if( length(x[["component"]]) > 0 ) {
    tmp <- lapply(x[["component"]], treeApply, func, ...)
    if(length(tmp) > 0)
      ans$children - tmp
  }

  # invoke the post-processing of component hook.
  if(length(post)) {
    post(x)
  }

  invisible(ans)
}

#creates an regression over biclustering methods, number of clusters and parameter of bicluster.
biClustSearch  <- function(intM, entrezName){
  outC1  <- NULL
  numClusterTest  <-  c(30,110,200)
  minrs <- c(5, 10, 25, 30 ,40)
  mincs <- c(5, 10, 25, 30 ,40)
  pvalueCutoff  <- 1e-3
  #biclustMethods  <- ("BCCC", "BCXmotifs", "BCPlaid", "BCSpectral", "BCBimax", "BCQuest")
  biclustMethods  <- c("BCCC", "BCBimax")
  for(numCluster in numClusterTest){
    #for(biclustmethod in biclustMethods){
    for(minr in minrs){
      for(minc in mincs){
	#if (biclustmethod == "BCBimax"){	
	#intMClust <-biclust(intM.binary, method=BCBimax(), minr=ceiling(5000/(1.5 * numCluster)), minc=ceiling(800/(1.5 * numCluster)) , number=numCluster)
	intMClust <-biclust(intM.binary, method=BCBimax(), minr=minr, minc=minc, number=numCluster)
	#} else if(biclustmethod == "BCCC"){
	#intMClust <-biclust(as.matrix(intM), method=BCCC(), delta=1.5,  alpha=1, number=numCluster)
	#intMClust <-biclust(as.matrix(intM), method=BCCC(), delta=1.5,  alpha=1, number=numCluster)
	#}
	out  <- goSignificantCluster(intMClust, intM, entrezName, pvalueCutoff=pvalueCutoff )
	while(out$sigBicluter == 1){
	  pvalueCutoff  <- pvalueCutoff/10
	  out  <- goSignificantCluster(intMClust, intM, entrezName, pvalueCutoff=pvalueCutoff )
	}
	out$runParam  <- list(biclustmethod,pvalueCutoff, numCluster )
	outC1  <- cbind(outC1 , out)
      }

    }
  }
  return(outC1)
}
biClustSearchMulticore  <- function(intM, entrezName, numClusterTest=c(25,100,250) ){
  outC1  <- NULL
  #numClusterTest  <-  c(25,100,250)
  minrs <- c(25)
  mincs <- c(2,3,4,6)
  pvalueCutoff  <- 1e-3
  for(numCluster in numClusterTest){
    print(paste("starting the processing for number of cluster", numCluster, sep=" "))
    for(minr in minrs){
      for(minc in mincs){
	intMClust <-biclust(intM.binary, method=BCBimax(), minr=minr, minc=minc, number=numCluster)
	out  <- goSignificantClusterMulticore(intMClust, intM, entrezName, pvalueCutoff=pvalueCutoff )
	#out  <- goSignificantCluster(intMClust, intM, entrezName, pvalueCutoff=pvalueCutoff )
	out$runParam  <- list("bimax",pvalueCutoff, numCluster )
	outC1  <- cbind(outC1 , out)
	save(file="outC1.RData", outC1)
      }
    }
    print(paste("finshed the processing for number of cluster", numCluster, sep=" "))
  }
  return(outC1)
}
