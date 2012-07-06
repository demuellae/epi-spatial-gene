## Avinash Das, 12 Jun 2012 
# this code process the interaction matrix present http://www.eurexpress.org. The interaction matrix was obtained on request to the author of paper http://www.plosbiology.org/article/info%3Adoi%2F10.1371%2Fjournal.pbio.1000582. The interaction data is at location ./data/eurexpress/matrix_1224667147767.txt. The interaction matrix contain data 
#10: strong, 4: medium, 1: weak, 0: no expression
#reading the file matrix.csv is comma separated with 3 lines removed from  matrix_1224667147767.txt        
intM <- read.table(file="../data/eurexpress/matrix.csv", header=TRUE, row.names=2, check.names=FALSE, sep=",")
# removing the first column as it have assay information. intM is now numeric 

intM<-intM[,-1:0]
Mdim <- dim(intM)
#dim(intM) (gene x tissue)
#[1] 5510  811

source("numAnno.R")
library(MASS)

#number of tissue where gene are strongly expressed
MnumAnno  <- apply(intM, 1, numAnno, val=10)
#histogram of number of tissues in which genes are strongly expressed")
pdf("../doc/stronglyExpressed.pdf")
truehist(MnumAnno, xlab="number of tissue with strong expression in gene")

#null distribution is B(n,p), with p is average number of tissue in which a gene is expressed
nullDist <- dbinom(seq(0, max(MnumAnno)),size=Mdim[2], prob=sum(MnumAnno)/(Mdim[1] * Mdim[2]))
lines(x=seq(0,max(MnumAnno)), y=nullDist, col=3)
dev.off()      

#number of tissue where gene are strongly or medium expressed
pdf("../doc/mediumExpressed.pdf")
MnumAnno  <- MnumAnno + apply(intM, 1, numAnno, val=4)
truehist(MnumAnno,xlab="number of tissue with strong expression in gene")
nullDist <- dbinom(seq(0, max(MnumAnno)),size=Mdim[2], prob=sum(MnumAnno)/(Mdim[1] * Mdim[2]))
lines(x=seq(0,max(MnumAnno)), y=nullDist, col=3)
dev.off()

#number of tissue where gene are strongly, medium or weakly expressed
pdf("../doc/weakExpressed.pdf")
MnumAnno  <- MnumAnno + apply(intM, 1, numAnno, val=1)
truehist(MnumAnno,xlab="number of tissue with strong expression in gene")
nullDist <- dbinom(seq(0, max(MnumAnno)),size=Mdim[2], prob=sum(MnumAnno)/(Mdim[1] * Mdim[2]))
lines(x=seq(0,max(MnumAnno)), y=nullDist, col=3)
dev.off()


#trying to visualize the matrix if there are some patterns. 
col<-gray(seq(256,0)/256)

pdf("../doc/interactionMatrix.pdf")
image(z=z<-as.matrix(intM), col=col, xlab="tissue type", ylab="gene expression")

image(z=z<-as.matrix(intM), xlab="tissue type", col=col,ylab="gene expression")
dev.off()

# PCA and sorting according to first principal component along tissue
pdf("../doc/interactionMatrixtissuePC1.pdf")
Mpca <- prcomp(intM)

Mpc1 <- intM[,order(Mpca$rotation[1,])]
image(z=z<-as.matrix(Mpc1), col=col, xlab="tissue type", ylab="gene exp")
dev.off()

# PCA and sorting according to first principal component along gene
pdf("../doc/interactionMatrixGenePC.pdf")
MpcaG <- prcomp(t(intM))
Mpc1G <- intM[order(MpcaG$rotation[1,]),]
image(z=z<-as.matrix(Mpc1G), col=col, xlab="tissue type", ylab="gene exp")
dev.off()
# PCA and sorting according to first principal component along gene and gene
pdf("../doc/interactionMatrixTissueGenePC.pdf")
Mpc1TG <- Mpc1[order(MpcaTG$rotation[1,]),]
image(z=z<-as.matrix(Mpc1G), col=col, xlab="tissue type", ylab="gene exp")
dev.off()


# hierarchical clustering along both rows and column of intM
pdf("../doc/interactionMatrixHirarchichal.pdf")
MclustG   <- hclust(as.dist((1 - cor(intM))/2))
MclustA   <- hclust(as.dist((1 - cor(t(intM)))/2))
hclustM <- intM[MclustA$order, MclustG$order]
image(z=z<-as.matrix(hclustM), col=col, xlab="anatomy", ylab="gene exp")
dev.off()






#finding GO annotations for a given gene set




#Convert MGI gene identifiers to entrez gene id
#source("http://bioconductor.org/biocLite.R")
#biocLite("org.Mm.eg.db") 
library(org.Mm.eg.db)
#xx  <- xx[!is.na(xx)]
# NOTE there are some duplicate entries and NA 
entrezName  <- as.list(org.Mm.egALIAS2EG)

#GO annontation example
library(GOstats); library(GO.db);  library(annotate) # Loads the required libraries.
library(RColorBrewer); library(xtable); library(Rgraphviz)
#library(ALL); library(genefilter); 
#goann <- as.list(GOTERM) # Retrieves full set of GO annotations.
#zz <- eapply(GOTERM, function(x) x@Ontology); table(unlist(zz)) # Calculates the number of annotations for each ontology category.

#library(ath1121501.db); affySample <- c("266592_at", "266703_at", "266199_at", "246949_at", "267370_at", "267115_s_at", "266489_at", "259845_at", "266295_at", "262632_at"); geneSample <- as.vector(unlist(mget(affySample, ath1121501ACCNUM, ifnotfound=NA))); library(ath1121501cdf); affyUniverse <- ls(ath1121501cdf); geneUniverse <- as.vector(unlist(mget(affyUniverse, ath1121501ACCNUM, ifnotfound=NA))); params <- new("GOHyperGParams", geneIds = geneSample, universeGeneIds = geneUniverse, annotation="ath1121501", ontology = "MF", pvalueCutoff = 0.5, conditional = FALSE, testDirection = "over"); hgOver <- hyperGTest(params); summary(hgOver); htmlReport(hgOver, file = "MyhyperGresult.html")

#geneSample <- as.vector(unlist(mget(entrezName, org.Mm.egACCNUM, ifnotfound=NA)));
#geneSample  <-  unlist(lapply(names(geneSample), function(x) (geneSample[[x]][1])))
#geneUniverse <- as.vector(unlist(mget(entrezName, org.Mm.egACCNUM, ifnotfound=NA)));
#geneUniverse  <-  unlist(lapply(names(geneUniverse), function(x) (geneUniverse[[x]][1])))
params <- new("GOHyperGParams", geneIds = as.vector(entrezName[1:1000]) , universeGeneIds = as.vector(entrezName), annotation="org.Mm.eg.db", ontology = "MF", pvalueCutoff = 0.5, conditional = FALSE, testDirection = "over"); 
hgOver <- hyperGTest(params); summary(hgOver); htmlReport(hgOver, file = "MyhyperGresult.html")

# alternative go annotation using GOHyperGAll function 

source("GOHyperGAll.txt")
readGOorg(myfile="../data/gene_association.mgi", colno=c(5,11,9), org="Mus muscus")
gene2GOlist(rootUK=T)

#some initial biclustering

library(biclust)

#creating a binary interaction matrix

#intM.binary <- as.matrix(intM)
#intM.binary[,] <- 0
#intM.binary[intM > 0] <- 1

intM.binary <- discretize(intM)

# BCXmotif
pdf("../doc/intMBCB.pdf")
intMBCB <-biclust(intM.binary, method=BCBimax(), minr = 100, minc=16, number=100)
#intMBCB <-biclust(intM.binary, method=BCXmotifs(), alpha=0.5, number=100) 
parallelCoordinates( x=intM.binary, bicResult=intMBCB, number=4)  
bubbleplot(intM.binary, intMBCB)
dev.off()


source("./myFunc.R")
outC  <- NULL
x  <- c(1e-6)
#choosing best threshold to obtain best clusters.
for(i in seq(1,1)){

	out  <- goSignificantCluster(intMBCB, intM, entrezName, pvalueCutoff=x[i] )
	outC  <- cbind(outC , out)
}


# temperary run on all clusters 
numCluster  <- 100
intMBCB <-biclust(intM.binary, method=BCBimax(), minr = ceiling(5000/(1.5 * numCluster)), minc=ceiling(800/(1.5 * numCluster)) , number=numCluster)
out  <- goSignificantCluster(intMBCB, intM, entrezName, pvalueCutoff=5e-6 )

pThershold  <- 1e-2 ;
sigBicluster  <- 0
for(i in seq(1, numCluster)){
	pVal  <- pvalues(out$Allcluster[[i]])
       #print(pVal[1])	
	print(sum(pVal  < pThershold))
	if (sum(pVal  < pThershold) > 0)
		sigBicluster  <- sigBicluster + 1
	print(sigBicluster)
}









#choosing the optimal number of clusters
numClusterTest  <-  c(25,50,100,150, 250)
outC  <- NULL
for(i in seq(1,length(numClusterTest))){
	numCluster  <-   numClusterTest[i]

	intMBCB <-biclust(intM.binary, method=BCBimax(), minr = ceiling(5000/(1.5 * numCluster)), minc=ceiling(800/(1.5 * numCluster)) , number=numCluster)
	out  <- goSignificantCluster(intMBCB, intM, entrezName, pvalueCutoff=5e-6 )
	outC  <- cbind(outC , out)
	print(i)
}

#numClusterTest  <-  c(25,50,100,150, 250)
numClusterTest  <-  c(1)
outC1  <- NULL
for(i in seq(1,length(numClusterTest))){
	numCluster  <-   numClusterTest[i]

	intMBCCC <-biclust(as.matrix(intM), method=BCCC(), delta=1.5,  alpha=1, number=numCluster)
	out  <- goSignificantCluster(intMBCCC, intM, entrezName, pvalueCutoff=5e-6 )
	outC1  <- cbind(outC1 , out)
	print(i)
}



#plotting beta function for the doc
colors <- c("red", "blue", "black")
x <- seq(0,1,by=0.01)
jpeg("../doc/beta.jpg")
plot(x, dbeta(x, 2, 40), type="l", ylab="PDF", lty=2)
lines(x, dbeta(x, 2, 5), lwd=2, col=colors[1])
lines(x, dbeta(x, 20, 2), lwd=2, col=colors[2])
labels <- c("active", "variable", "inactive")
legend("topright", inset=.05, title="Beta priors",
       labels, lwd=2, lty=c(1, 1, 2), col=colors)
dev.off()



# xml tree parser
library(XML)
tissue.tree <- xmlTreeParse("../data/EMAPtree.xml", 
			    handlers=list(anatomy=function(x,attr) {x},
					  component=function(x,attr) {x}, 
					  childrenId=function(x,attr) {x}, 
					  startElement=function(x,attr){ NULL}),
			    asTree=T)
v <- treeApply(x$children, function(x) cat(class(x),"\n"))
