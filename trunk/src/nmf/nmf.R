if(FALSE){
#intM <- read.table(file="../data/eurexpress/matrix.csv", header=TRUE, row.names=2, check.names=FALSE, sep=",")
#intM<-intM[,-1:0]
#Mdim <- dim(intM)i
#perform 10% instead of max 2) nmf ns version
library(NMF)
library(doMC)
library(foreach)
library(bigmemory)
library(synchronicity)
require(biclust)
getEntrezname <- function()
{
  library(org.Mm.eg.db)
  entrezName  <- as.list(org.Mm.egALIAS2EG)
  detach(package:org.Mm.eg.db)
  #detach(package:Category)
  #detach(package:Category)
  detach(package:AnnotationDbi, force=T)
  return(entrezName)
}
findConsecutiveMax <- function(b, num)
{
  col = b[,num]
  b.order=order(col, decreasing=T)
  a=b[b.order,]
 out = rep(FALSE, dim(a)[1])
 colmean = mean(col) 
  for(i in 1:length(out)){
    if((max(a[i,])==a[i,num]) && (a[i,num]> colmean) ){
      #out[i]= TRUE
      imax = i
    }
  }
 if(imax > 0) out[1:imax] = TRUE
 out = out[order(b.order, decreasing=T)]
 out
}

findMax <- function(b, num)
{
  col = b[,num]
  #a=b[order(col, decreasing=T),]
a = b
 out = rep(FALSE, dim(a)[1])
 colmean = mean(col) 
  for(i in 1:length(out)){
    if((max(a[i,])==a[i,num]) && (a[i,num]> colmean) ){
      out[i]= TRUE
      #imax = i
    }
  }
 #if(imax > 0) out[1:imax] = TRUE
 out
}
number=50
meth ="nsNMF"
cmd = sprintf("res = nmf(as.matrix(intM), rank=%d, meth=meth, nrun=100, .pbackend=25, .opt='vP')", number) 
#if(TRUE){
  #eval(parse(text=cmd))
  #save(file=paste("nmf/res.",meth ,".RData", sep=""), res)
#}else{
  #load("nmf/res.RData")
#}
print("here")
w = basis(res)
h = coef(res)
clust = new("Biclust")
Parameters = list()
Parameters$Call = cmd 
Parameters$Method = "NMF"   
clust@Parameters = Parameters 
clust@RowxNumber = matrix(FALSE, Mdim[1], number) 
clust@NumberxCol = matrix(FALSE, number, Mdim[2]) 
clust@Number = number
  for(num in 1:number){
    clust@RowxNumber[ ,num]= findMax(b=w, num=num)
    clust@NumberxCol[num, ]= t(findMax(t(h), num ))
  }
clust@RowxNumber[1:10,1:10]
#load(file="../data/tissuetab.RData")
#tissueType = colnames(intM)[clust@NumberxCol[clust, ]] 
#tissueTypeEmapTerm =tissuetab$emap_term[tissuetab$emap_id%in%tissueType]
source("myFunc.R")
#debug(goSignificantClusterMulticore)
print("here1")
entrezName = getEntrezname()
out = goSignificantClusterMulticore(clust, intM, entrezName, pvalueCutoff = 0.5)






}
writeClusterFile  <-  function(intMClust, intM, novartis=FALSE, genetab.common = genetab.common, biclusterDir="biclusters/", sepSymbol=":", numCluster = 25, threshold =0){
  library(biclust)
  if (!file.exists(biclusterDir)) dir.create(biclusterDir)
  #discretize
  intM.binary <- as.matrix(intM)
  intM.binary[,] <- 0
  intM.binary[intM > threshold] <- 1
  Mdim <- dim(intM)
  #intMClust = biclust(intM.binary, method=BCBimax(), minr = 25, minc=5 , number=numCluster)
  #out  <- goSignificantCluster(intMBCB, intM, entrezName, pvalueCutoff=5e-6 )

  geneUniverse = rownames(intM)
  if(novartis){
    geneIdNames = "mgi_symbol"
    ext <- ".Nova"
  }else{
    geneIdNames = "tmp_gene_symbol"
    ext  <- ""
  }

  background.table =  genetab.common[genetab.common[,geneIdNames] %in% geneUniverse,]
  background.table = background.table[!duplicated(background.table$ensembl_gene_id)] 
  background = data.frame(regionId = paste("chr",background.table$chromosome_name,sepSymbol, background.table$start_position,sepSymbol,background.table$end_position,sep=""), 
			  ensemblId = background.table$ensembl_gene_id)
  filename  <- paste(biclusterDir,"/","background",ext,".txt",sep="")
  write.table(file=filename,x=background, quote=F, row.names=F)
  for(clust in seq(1,intMClust@Number)){
    geneSample  <- rownames(intM)[intMClust@RowxNumber[,clust]]
    bimax.table =  genetab.common[genetab.common[,geneIdNames] %in% geneSample,] 
  bimax.table = bimax.table[!duplicated(bimax.table$ensembl_gene_id)] 
    bimax = data.frame(regionId = paste("chr",bimax.table$chromosome_name,sepSymbol,bimax.table$start_position,sepSymbol,bimax.table$end_position,sep=""), ensemblId = bimax.table$ensembl_gene_id)
    filename  <- paste(biclusterDir,"/","bimax",numCluster,"C.",clust, ext,".txt",sep="")
    write.table(file=filename,x=bimax, quote=F, row.names=F)
  }
  #return(intMClust)
}

  date = Sys.Date()
outputDir = "../result/" 
if (!file.exists(outputDir)) dir.create(outputDir)
outputDir = paste(outputDir,date,"manual_enrichment_analyses",sep="") 
if (!file.exists(outputDir)) dir.create(outputDir)
biclusterDir = paste(outputDir,"bicluster.nsNMF",sep="/") 
if (!file.exists(biclusterDir)) dir.create(biclusterDir)
writeClusterFile(intMClust=clust, intM=intM, genetab.common=genetab.common, biclusterDir=biclusterDir)
