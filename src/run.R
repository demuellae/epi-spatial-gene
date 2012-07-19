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

# alternative go annotation using GOHyperGAll function 



library(biclust)

#creating a binary interaction matrix

intM.binary <- as.matrix(intM)
intM.binary[,] <- 0
intM.binary[intM > 0] <- 1
#intM.binary <- discretize(intM)

#pdf("../doc/intMBCB.pdf")
#intMBCB <-biclust(intM.binary, method=BCBimax(), minr = 100, minc=16, number=100)
##intMBCB <-biclust(intM.binary, method=BCXmotifs(), alpha=0.5, number=100) 
#parallelCoordinates( x=intM.binary, bicResult=intMBCB, number=4)  
#bubbleplot(intM.binary, intMBCB)
#dev.off()


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
#numClusterTest  <-  c(25,50,100,150, 250)
numClusterTest  <-  c(25)
outC  <- NULL
for(i in seq(1,length(numClusterTest))){
	numCluster  <-   numClusterTest[i]

	#intMBCB <-biclust(intM.binary, method=BCBimax(), minr = ceiling(5000/(1.5 * numCluster)), minc=ceiling(800/(1.5 * numCluster)) , number=numCluster)
	intMBCB <-biclust(intM.binary, method=BCBimax(), minr = 10,
			  minc=10,
			  number=numCluster)
	out  <- goSignificantCluster(intMBCB, intM, entrezName, pvalueCutoff=5e-6 )
	outC  <- cbind(outC , out)
	print(i)
}

#numClusterTest  <-  c(25,50,100,150, 250)

outC1  <- biClustSearch(intM, entrezName)
outC3  <- biClustPSQM(intM, entrezName)

bimax25C1.table <- genetabNAN[genetabNAN$entrezgene %in% geneIds(temp[[1]][[2]]),]
bimax25C1 <- data.frame(regionId = paste("chr",bimax25C1.table$chromosome_name,"_",bimax25C1.table$start_position,"_",bimax25C1.table$end_position,sep=""), ensemblId = bimax25C1.table$ensembl_gene_id)

write.table(file="../data/biclusters/bimax25C1.txt",x=bimax25C1, quote=F, row.names=F)

background.table <- genetabNAN[genetabNAN$entrezgene %in% geneIds(temp[[1]][[2]]),]
background <- data.frame(regionId = paste("chr",background.table$chromosome_name,"_",background.table$start_position,"_",background.table$end_position,sep=""), ensemblId = background.table$ensembl_gene_id)

write.table(file="../data/biclusters/background.txt",x=background, quote=F, row.names=F)


#analysis of novartis dataset
#save("file=genetab.novartis.RData",genetab.novartis)
#load("file=genetab.novartis.RData")
load("./dataset/novartis/downloaded/GSE1133_GPL1073_GNF1M_mouse.RData")
library(Biobase)
library(BiocGenerics)
intNovartis = exprs(mouse)
#GNF1M = parse.table("./dataset/processed/gnf1m.annot2007.tsv")
#save(file="./GNF1M_annotMapping.RData",GNF1M)
load(file="./GNF1M_annotMapping.RData")
intNovartis.annotated  <- intNovartis[rownames(intNovartis) %in% rownames(GNF1M),]
intNovartis.symbol  <- intNovartis.annotated[!is.na(GNF1M[match(rownames(intNovartis.annotated), rownames(GNF1M)) , "Symbol"]),] 
symbol = GNF1M[match(rownames(intNovartis.symbol), rownames(GNF1M)) , "Symbol"] 
rownames(intNovartis.symbol)  <- symbol
thresholdQuantile = sum(intM.binary ==1)/4468610
intNovartis.sort = sort(unlist(as.list(intNovartis.symbol)), decreasing=T)
novartisExprThereshold = intNovartis.sort[ceiling(length(intNovartis.sort)*thresholdQuantile)]
intNovartis.binary <- as.matrix(intNovartis.symbol)
intNovartis.binary[,] <- 0
intNovartis.binary[intNovartis.symbol > novartisExprThereshold] = 1
intMClust = biclust(intNovartis.binary, method=BCBimax(), minr = 10, minc=10 , number=25)
symbol2entrezName  <- as.list(org.Mm.egSYMBOL2EG)
out  <- goSignificantCluster(intMClust, intNovartis.symbol, entrezName=symbol2entrezName, pvalueCutoff=5e-6 )



pThershold  <- c(1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-4, 5e-5, 1e-5, 5e-6,1e-6 )
numClust=100
intMClust = biclust(intNovartis.binary, method=BCBimax(), minr = 10, minc=10 , number=numClust)
out.ish  <- goSignificantCluster(intMClust, intNovartis.symbol, entrezName=symbol2entrezName, pvalueCutoff=5e-6 )
intMBCB <-biclust(intM.binary, method=BCBimax(), minr = 10, minc=10, number=numClust)
out.nova  <- goSignificantCluster(intMBCB, intM, entrezName, pvalueCutoff=5e-6 )
numSig.ish  <- NULL
numSig.nova  <- NULL
#novaClustNum <- intMClust@Number
novaClustNum <- 20
ishClustNum <- 20
for(i in seq(1, length(pThershold))){
numSig <- numSigCluster(out=out.ish, Number=ishClustNum, pThershold=pThershold[i])
numSig.ish  <- c(numSig.ish, numSig)
numSig <- numSigCluster(out=out.nova, Number=novaClustNum, pThershold=pThershold[i])
numSig.nova  <- c(numSig.nova, numSig)
}
jpeg(paste("ISHvsNOVA",numClust,".jpg",sep="")) 
counts <- rbind(numSig.ish/numClust, numSig.nova/numClust)
colnames  <- pThershold
barplot(counts, main=paste("% of significant bicluster for ", numClust," clusters", sep=""),
  xlab="p-value", col=c("darkblue","red"),
  legend = c("% of ISH significant cluster","% of Nova significant cluster"), beside=TRUE)
dev.off()



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
