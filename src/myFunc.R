goSignificantCluster  <- function(cluster, geneMat, entrezName ){  
	geneUniverse  <- entrezName[rownames(intM)]
	geneUniverse.entrez  <- unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1])))

	frac <- 0
	for(i in seq(1,cluster@Number)){

		geneSample  <- entrezName[rownames(geneMat)[intMBCB@RowxNumber[,i]]]
		geneSample.entrez  <- 	unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1])))
		params <- new("GOHyperGParams", geneIds = as.vector(geneSample.entrez) , universeGeneIds = as.vector(geneUniverse.entrez), annotation="org.Mm.eg.db", ontology = "MF", pvalueCutoff = 0.005, conditional = TRUE, testDirection = "over") 
		hgOver <- hyperGTest(params) 
		df  <-  summary(hgOver)
		if (dim(df)[1] > 0 ){
			frac <- frac + 1
		}
	}
	return(frac/(cluster@Number))
}

