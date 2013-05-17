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


# multicore environment
library(org.Mm.eg.db)
library(GOstats); library(GO.db);  library(annotate)
library(biclust)
options(cores = 25)
source("../myFunc.R")
biClustSearchMulticore(intM, entrezName, numClusterTest=c(25,100,250) )


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
#library(XML)
#tissue.tree <- xmlTreeParse("../data/EMAPtree.xml", 
			    #handlers=list(anatomy=function(x,attr) {x},
					  #component=function(x,attr) {x}, 
					  #childrenId=function(x,attr) {x}, 
					  #startElement=function(x,attr){ NULL}),
			    #asTree=T)
#v <- treeApply(x$children, function(x) cat(class(x),"\n"))

library(XML)

URL <- "http://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/Stage_xml/TS23.xml"
root <- xmlTreeParse(URL, useInternalNodes = TRUE)

fn <- function(node) {
   id <- xmlAttrs(node)["id"]
   parent.id <- xmlAttrs(xmlParent(node))["id"]
   setNames(head(c(id, parent.id, NA), 2), c("id", "parent"))
}
parents.all <- t(xpathSApply(root, "//component", n))

parents = as.data.frame(parents.all[1:1411,])
findLeaf <- function(xx) which(!(xx$id %in% xx$parent))
temp =1:nrow(parents)
inx = list()
while(TRUE){
  inxT = findLeaf(parents[temp,])
  inx = c(inx,inxT)
  inx = c

}
findLeaf(parents)

motifWordcloud <- function(fileDir=".")
{
  numClust = 50 
  require(RColorBrewer)
  require(wordcloud)
  for(i in seq(1,numClust)){	
    jpeg(paste(fileDir, "/", "promoternmf50C.", i,".jpg", sep="")) 
    fishertab <- read.table(file=paste(fileDir, "/", "promoternmf50C.", i, ".bed.fisher", sep=""),header=FALSE, sep="\t")

    fishertab = fishertab[ fishertab[,7] < 0.05, ]
    fishertab = fishertab[!duplicated(fishertab[,1]), ] 
    fishertab$logq=-log(fishertab[,7] + 1e-1000)
    pal2 <- brewer.pal(8,"Dark2")
    wordcloud(fishertab[,1] , freq=fishertab$logq, scale=c(2,.01),
	      colors=pal2, max=50)
    dev.off()
  }
}



#this are function for dendogram analysis
P = read.table(file="../result/2013-05-08/P", header = F, sep=",", comment.char = "", stringsAsFactors =F,  strip.white =T)
T1 = read.table(file="../result/2013-05-08/Time", header = F, sep=",", comment.char = "", stringsAsFactors =F,  strip.white =T)
XX = read.table(file="../result/2013-05-08/XXcomplete", header = F, sep=",", comment.char = "", stringsAsFactors =F,  strip.white =T)
a <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
a$merge = matrix(0, 49,2)
for(ii in 1:49){
  temp = which(P==ii+50) - 50
  temp[temp <1] = temp[temp <1] -1
  a$merge[ii,] = temp
}
a$height = abs(T1[51:99])
#a$merge <- matrix(c(-1, -2,
                    #-3, -4,
                     #1,  2), nc=2, byrow=TRUE ) 
#a$height <- c(1, 1.5, 3)    # define merge heights
a$order <- 50:1              # order of leaves(trivial if hand-entered)
a$labels <- paste("a", 50:1, sep="")   # labels of leaves
class(a) <- "hclust"        # make it an hclust object
plot(a)                     # look at the result   

#convert to a dendrogram object if needed
ad <- as.dendrogram(a)


#creating hmm  object
environment(BaumWelch.dthmm.Tree) = asNamespace("HiddenMarkov")
environment(Estep.Tree) = asNamespace("HiddenMarkov")
environment(forwardback.Tree) = asNamespace("HiddenMarkov")

low = 1e-5 
high = 1 - low
medium = 0.1
object = NULL
P = unlist(P)
XXcur = unlist(XX[4780,])
object$x = XXcur
Pi = matrix(1/3,9,3) 
numLeaf = (1+length(XXcur))/2
delta = matrix(0, numLeaf, 3)
delta[XXcur[1: numLeaf] ==0,1] = 1
delta[XXcur[1: numLeaf] ==1,3] = 1
#delta = XXcur[1: numLeaf]
XXcur.approx = XXcur
XXcur.approx[XXcur==1] = .95
XXcur.approx[XXcur==0] = .05
tree.object <- dthmm(XXcur.approx, Pi, delta, "beta", list(shape1=c(1, 2, 3), shape2=c(3,5, 1)))
tree.object$P = P
control = bwcontrol(tol = 1e-6, posdiff = F, converge = expression(abs(diff) < tol))
tree.out = BaumWelch.dthmm.Tree(tree.object, control)
curve(dbeta(x,shape1=tree.out$pm$shape1[1], shape2=tree.out$pm$shape2[1]), from=0, to=1)
curve(dbeta(x,shape1=tree.out$pm$shape1[2], shape2=tree.out$pm$shape2[2]), from=0, to=1, add=T, col="red")
curve(dbeta(x,shape1=tree.out$pm$shape1[3], shape2=tree.out$pm$shape2[3]), from=0, to=1, add=T, col="green")
#toy example
x.toy = c(1,1,0,1,0.5)
P.toy = c(4,4,5,5,0)
numLeaf = (1+length(x.toy))/2
delta = matrix(0, numLeaf, 3)
delta[x.toy[1: numLeaf] ==0,1] = 1
delta[x.toy[1: numLeaf] ==1,3] = 1
#delta = x.toy[1: numLeaf]
x.toy.approx = x.toy
x.toy.approx[x.toy==1] = .95
x.toy.approx[x.toy==0] = .05
tree.object <- dthmm(x.toy.approx, Pi, delta, "beta", list(shape1=c(1, 0.5, 3), shape2=c(3,0.5, 1)))
tree.object$P = P.toy
tree.out = BaumWelch.dthmm.Tree(tree.object)



