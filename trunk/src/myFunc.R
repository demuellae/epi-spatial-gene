goSignificantCluster  <- function(cluster, geneMat, entrezName, pvalueCutoff = 0.5 ){  
	geneUniverse  <- entrezName[rownames(intM)]
	geneUniverse.entrez  <- unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1])))

	sigBicluter <- 0
	Allcluster  <- NULL 
	sigGO  <- NULL
	#for(i in seq(1,1)){
	for(i in seq(1,cluster@Number)){
		print(i)
		geneSample  <- entrezName[rownames(geneMat)[intMBCB@RowxNumber[,i]]]
		geneSample.entrez  <- 	unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1])))
		params <- new("GOHyperGParams", geneIds = as.vector(geneSample.entrez) , universeGeneIds = as.vector(geneUniverse.entrez), annotation="org.Mm.eg.db", ontology = "MF", pvalueCutoff = pvalueCutoff, conditional = TRUE, testDirection = "over")
		hgOver <- hyperGTest(params)
		Allcluster  <- c(Allcluster, hgOver)	
		df  <-  summary(hgOver)
		if (dim(df)[1] > 0 ){
			sigBicluter <- sigBicluter + 1
		}
		print(sigBicluter)
		##pVal  <- pvalues(hgOver)
		#g  <- names(pVal)[pvalues(hgOver)  > pvalueCutoff] 
		#g > 
		sigGO  <- c(sigGO, sigCategories(hgOver))

	}
	sigBicluter  <- (sigBicluter/(cluster@Number))
	sigGO  <- unique(sigGO)
	out  <- NULL
	out$sigBicluter  <- sigBicluter
	out$sigGO  <- length(sigGO)
	out$Allcluster  <-  Allcluster
	return(out)
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

biClustSearch  <- function(intM, entrezName){
	#creates an regression over biclustering methods, number of clusters and parameter of bicluster.
	numClusterTest  <-  c(25,50,100,250)
	pvalueCutoff  <- 1e-3
	#biclustMethods  <- ("BCCC", "BCXmotifs", "BCPlaid", "BCSpectral", "BCBimax", "BCQuest")
	biclustMethods  <- ("BCCC", "BCBimax")
	for(numCluster in numClusterTest){
		for(biclustmethod in biclustMethods){
			if (biclustmethod == "BCBimax"){	
				intMClust <-biclust(intM.binary, method=BCBimax(), minr=ceiling(5000/(1.5 * numCluster)), minc=ceiling(800/(1.5 * numCluster)) , number=numCluster)
			} else if(biclustmethod == "BCCC"){
				intMClust <-biclust(as.matrix(intM), method=BCCC(), delta=1.5,  alpha=1, number=numCluster)
			}
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
