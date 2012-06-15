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
truehist(MnumAnno)
#null distribution is B(n,p), with p is average number of tissue in which a gene is expressed
nullDist <- dbinom(seq(0, max(MnumAnno)),size=Mdim[2], prob=sum(MnumAnno)/(Mdim[1] * Mdim[2]))
lines(x=seq(0,max(MnumAnno)), y=nullDist, col=3)


#number of tissue where gene are strongly or medium expressed
MnumAnno  <- MnumAnno + apply(intM, 1, numAnno, val=4)
truehist(MnumAnno)
nullDist <- dbinom(seq(0, max(MnumAnno)),size=Mdim[2], prob=sum(MnumAnno)/(Mdim[1] * Mdim[2]))
lines(x=seq(0,max(MnumAnno)), y=nullDist, col=3)

#number of tissue where gene are strongly, medium or weakly expressed
MnumAnno  <- MnumAnno + apply(intM, 1, numAnno, val=1)
truehist(MnumAnno)
nullDist <- dbinom(seq(0, max(MnumAnno)),size=Mdim[2], prob=sum(MnumAnno)/(Mdim[1] * Mdim[2]))
lines(x=seq(0,max(MnumAnno)), y=nullDist, col=3)

#trying to visualize the matrix if there are some patterns. 
image(z=z<-as.matrix(intM), col=col, xlab="anatomy", ylab="gene exp")

# PCA and sorting according to first principal component 
Mpca <- prcomp(intM)

Mpc1 <- intM[,order(Mpca$rotation[1,])]
image(z=z<-as.matrix(Mpc1), col=col, xlab="anatomy", ylab="gene exp")


# PCA and sorting according to first principal component 
MpcaG <- prcomp(t(intM))

Mpc1G <- intM[order(MpcaG$rotation[1,]),]


# hierarchical clustering along both rows and column of intM
MclustG   <- hclust(as.dist((1 - cor(intM))/2))
MclustA   <- hclust(as.dist((1 - cor(t(intM)))/2))
hclustM <- intM[MclustA$order, MclustG$order]
image(z=z<-as.matrix(hclustM), col=col, xlab="anatomy", ylab="gene exp")
