goSignificantCluster  <- function(cluster, geneMat, entrezName, pvalueCutoff = 0.5 ){ 
  library(biclust)
  library(GOstats)
  library(org.Mm.eg.db)
  geneUniverse  <- entrezName[rownames(intM)]
  geneUniverse.entrez  <- unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1])))

  sigBicluter <- 0
  Allcluster  <- NULL 
  sigGO  <- NULL
  GOcategories  <- c("BP","MF", "CC")
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
goClust <- function(i, geneMat=geneMat, cluster=cluster, entrezName=entrezName, pvalueCutoff)
{
  GOcategories  <- c("BP","MF", "CC")
  geneSample  <- entrezName[rownames(geneMat)[cluster@RowxNumber[,i]]]
  geneSample.entrez  <- 	unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1])))
  goCatSigCnt <- 0
  Allcluster <- NULL
  for(goCategory in GOcategories){
    params <- new("GOHyperGParams", geneIds = as.vector(geneSample.entrez) , universeGeneIds = as.vector(geneUniverse.entrez), annotation="org.Mm.eg.db", ontology = goCategory, pvalueCutoff = pvalueCutoff, conditional = TRUE, testDirection = "over")
    hgOver <- hyperGTest(params)
    Allcluster  <- c(Allcluster, hgOver)	
    df  <-  summary(hgOver)
    if (dim(df)[1] > 0 )
      goCatSigCnt  <- goCatSigCnt + 1

    #print(sigBicluter)
    #sigGO  <- c(sigGO, sigCategories(hgOver))
  }
  #if (goCatSigCnt > 0) 
  #sigBicluter <- sigBicluter + 1
  return(Allcluster)

}

goSignificantClusterMulticore  <- function(cluster, geneMat, entrezName, pvalueCutoff = 0.5 ){ 
  library(biclust)
  library(GOstats)
  library(org.Mm.eg.db)
  library(multicore)
  geneUniverse  <- entrezName[rownames(intM)]
  geneUniverse.entrez  <- unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1])))

  sigBicluter <- 0
  Allcluster  <- NULL 
  sigGO  <- NULL
  #for(i in seq(1,cluster@Number))
  clustseq <- seq(1,cluster@Number)
  Allcluster  <- mclapply(clustseq, goClust, cluster=cluster, entrezName=entrezName, pvalueCutoff ) 
  sigBicluter  <- (sigBicluter/(cluster@Number))
  #sigGO  <- unique(sigGO)
  out  <- NULL
  #out$sigBicluter  <- sigBicluter
  #out$sigGO  <- length(sigGO)
  out$Allcluster  <-  Allcluster
  return(out)
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
  numClusterTest  <-  c(25,100,250)
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
biClustSearchMulticore  <- function(intM, entrezName){
  outC1  <- NULL
  numClusterTest  <-  c(25,100,250)
  minrs <- c(5, 10, 25, 30 ,40)
  mincs <- c(5, 10, 25, 30 ,40)
  pvalueCutoff  <- 1e-3
  for(numCluster in numClusterTest){
    for(minr in minrs){
      for(minc in mincs){
	intMClust <-biclust(intM.binary, method=BCBimax(), minr=minr, minc=minc, number=numCluster)
	out  <- goSignificantCluster(intMClust, intM, entrezName, pvalueCutoff=pvalueCutoff )
	out$runParam  <- list(biclustmethod,pvalueCutoff, numCluster )
	outC1  <- cbind(outC1 , out)
      }
    }
  }
  return(outC1)
}
