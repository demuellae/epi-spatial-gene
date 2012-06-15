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
pdf("../doc/interactionMatrix.pdf")
image(z=z<-as.matrix(intM), col=heat.colors(11), xlab="tissue type", ylab="gene expression")
dev.off()

# PCA and sorting according to first principal component along tissue
pdf("../doc/interactionMatrixtissuePC1.pdf")
Mpca <- prcomp(intM)

Mpc1 <- intM[,order(Mpca$rotation[1,])]
image(z=z<-as.matrix(Mpc1), col=heat.colors(11), xlab="tissue type", ylab="gene exp")
dev.off()
col<-heat.colors(11)
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


