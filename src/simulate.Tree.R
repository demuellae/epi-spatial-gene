# generate tree
#t <- function(p)
  #list(if(runif(1) < p) t(p) else 0, if(runif(1) < p) t(p) else 0)
numLeaf <- 100
n = 2*numLeaf -1 
h <- ceiling(log(numLeaf,base=2) +1)
P  <-  rep(0,n)
Pi = matrix(1/3,9,3) 

Pi[1,] <- c(.9996, .0002, .0002)
Pi[2,] <- c(.0002, .9996, .0002)
Pi[3,] <- c(.0002, .9996, .0002)
Pi[4,] <- c(.0002, .9996, .0002)
Pi[5,] <- c(.0002, .9996, .0002)
Pi[6,] <- c(.0002, .9996, .0002)
Pi[7,] <- c(.0002, .9996, .0002)
Pi[8,] <- c(.0002, .9996, .0002)
Pi[9,] <- c(.0002, .0002, .9996)

curr.nodes <-  seq(numLeaf); curr.parent <- numLeaf 
for(jj in seq(numLeaf - 1)){
  curr.child  <- sample(curr.nodes, replace=F, size=2)
  curr.parent  <- curr.parent + 1
  P[curr.child] <- curr.parent 
  curr.nodes <- curr.nodes[!(curr.nodes %in% curr.child)]
  curr.nodes <- c(curr.nodes, curr.parent)
}

#B  <- findBro(P)
#generate z

z <- rep(0, n)
z[n] <- 0 # initial state
Y <- rep(0, n)
shape <- matrix(0, 3,2)
shape[1,]= c(1, 5)
shape[2,] = c(5, 5)
shape[3,] = c(5, 1)

for(child in seq(1, numLeaf)) {
  #Randomly select z for the leaves
  z[child] <- sample(c(-1, 1), size=1, prob=c(.5,.5))
  Y[child] <- rbeta(1,shape1=shape[z[child]+2,1], shape2= shape[z[child]+2,2])
}
z[1:numLeaf] <- 1
z[1:numLeaf][Y[1:numLeaf] < .5] <- -1

for(parent in seq(numLeaf+1, n)){
  curr.child <- which(P==parent)
  z[parent] <- sample(-1:1, size=1, prob=Pi[z[curr.child[1]]+2 + 3*(z[curr.child[2]]+1),], replace=T)
}

for(node in seq(numLeaf+1, n)){
  Y[node]  <- rbeta(1,shape1=shape[z[node]+2,1], shape2= shape[z[node]+2,2])
}

source("myR/BaumWelch.dthmm.Tree.R")
source("myR/Estep.Tree.R")
source("myR/forwardback.Tree.R")
source("myR/makedensity.R")
source("myR/getj.R")
source("myR/dthmm.R")
delta <- matrix(0,n,3); delta[Y > .5, 3] <- 1; delta[Y <= .5,1] <- 1
#tree.object <- dthmm(Y, Pi, delta, "beta", list(shape1=c(1, 2, .5), shape2=c(3,5, .7)))
tree.object <- dthmm(Y, Pi, delta, "beta", list(shape1=shape[,1], shape2=shape[,2]))
tree.object$P = P
control = bwcontrol(tol = 1e-6, posdiff = F, converge = expression(abs(diff) < tol))
tree.out = BaumWelch.dthmm.Tree(tree.object, control)

write.table(file="P",x = c(P, 0) , row.names = F, col.names =F,  sep=",", quote=F )
write.table(file="z",x = z, row.names = F, col.names =F,  sep=",", quote=F )
write.table(file="Y",x = Y, row.names = F, col.names =F,  sep=",", quote=F )


##### binomial distribution
numLeaf <- 4
n = 2*numLeaf -1 
h <- ceiling(log(numLeaf,base=2) +1)
P  <-  rep(0,n)
Pi = matrix(1/3,9,3) 

Pi[1,] <- c(.9996, .0002, .0002)
Pi[2,] <- c(.0002, .9996, .0002)
Pi[3,] <- c(.0002, .9996, .0002)
Pi[4,] <- c(.0002, .9996, .0002)
Pi[5,] <- c(.0002, .9996, .0002)
Pi[6,] <- c(.0002, .9996, .0002)
Pi[7,] <- c(.0002, .9996, .0002)
Pi[8,] <- c(.0002, .9996, .0002)
Pi[9,] <- c(.0002, .0002, .9996)

curr.nodes <-  seq(numLeaf); curr.parent <- numLeaf 
for(jj in seq(numLeaf - 1)){
  curr.child  <- sample(curr.nodes, replace=F, size=2)
  curr.parent  <- curr.parent + 1
  P[curr.child] <- curr.parent 
  curr.nodes <- curr.nodes[!(curr.nodes %in% curr.child)]
  curr.nodes <- c(curr.nodes, curr.parent)
}

#B  <- findBro(P)
#generate z

z <- rep(0, n)
z[n] <- 0 # initial state
Y <- rep(0, n)
shape <- c(.02, .3, .7) 
for(child in seq(1, numLeaf)) {
  #Randomly select z for the leaves
  z[child] <- sample(c(-1, 1), size=1, prob=c(.5,.5))
  Y[child] <- rbinom(1,1, prob=shape[z[child]+2]) 
}
z[1:numLeaf] <- 1
z[1:numLeaf][Y[1:numLeaf] < .5] <- -1

for(parent in seq(numLeaf+1, n)){
  curr.child <- which(P==parent)
  z[parent] <- sample(-1:1, size=1, prob=Pi[z[curr.child[1]]+2 + 3*(z[curr.child[2]]+1),], replace=T)
}

for(node in seq(numLeaf+1, n)){
  Y[node]  <- rbinom(1,1,prob=shape[z[node]+2] ) 
}

source("myR/BaumWelch.dthmm.Tree.R")
source("myR/Estep.Tree.R")
source("myR/forwardback.Tree.R")
source("myR/makedensity.R")
source("myR/getj.R")
source("myR/dthmm.R")
source("myR/bwcontrol.R")
source("myR/Mstep.binom.R")

Pi = matrix(1/3,9,3) 

delta <- matrix(0,n,3); delta[Y > .5, 3] <- 1; delta[Y <= .5,1] <- 1
#tree.object <- dthmm(Y, Pi, delta, "beta", list(shape1=c(1, 2, .5), shape2=c(3,5, .7)))
tree.object <- dthmm(Y, Pi, delta, "binom", list(prob=shape), list(size=matrix(1, nrow=1, ncol=n)))
tree.object$P = P
control = bwcontrol(tol = 1e-6, posdiff = F, converge = expression(abs(diff) < tol))

write.table(file="P",x = c(P, 0) , row.names = F, col.names =F,  sep=",", quote=F )
write.table(file="z",x = z, row.names = F, col.names =F,  sep=",", quote=F )
write.table(file="Y",x = Y, row.names = F, col.names =F,  sep=",", quote=F )

tree.out = BaumWelch.dthmm.Tree(tree.object, control)




# tissue-tree infomation
#library(XML)

#URL <- "http://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/Stage_xml/TS23.xml"
#root <- xmlTreeParse(URL, useInternalNodes = TRUE)

fn <- function(node) {
   id <- xmlAttrs(node)["name"]
   parent.id <- xmlAttrs(xmlParent(node))["name"]
   setNames(head(c(id, parent.id, NA), 2), c("id", "parent"))
}
parents.all <- t(xpathSApply(root, "//component", fn))

parents = as.data.frame(parents.all[1:1411,])
#findLeaf <- function(xx) which(!(xx$id %in% xx$parent))
#temp =1:nrow(parents)
#nodes <- NULL
#P <- NULL
#while(TRUE){
#  inxT = findLeaf(parents[temp,])
#  nodes <- c(nodes,as.character(parents[inxT, "id"]))
#  P <- c(P ,as.character(parents[inxT, "parent"]))
#  temp <- temp[-inxT]
#  if(length(temp) == 0) break;
#}

#P <- c(10,10,11,11,12,12,12,12,13,11,14,14,14,0)
# binary tree simulation
binaryP <- P
binaryY <- Y
numLeaf <- 9

"binaryTree" <- function(P, Y, numLeaf) {
  binaryP <- P
  binaryY <- Y
  offset <- 0
  childOffset <- mat.or.vec(length(P), 1)
  for (parent in seq(numLeaf+1, length(P))) {
      children <- which(P == parent)
      binaryP[children[1] + childOffset[children[1]]] <- parent + offset
      binaryY[children[1] + childOffset[children[1]]] <- Y[children[1]]
      if (length(children) > 1) {
      	 binaryP[children[2] + childOffset[children[2]]] <- parent + offset
	 binaryY[children[2] + childOffset[children[2]]] <- Y[children[2]]
      }
      childOffset[parent] <- offset
      if (length(children) > 2) {
      	 length(binaryP) <- length(binaryP) + length(children) - 2
	 for (idx in seq(3, length(children))) {
	     offset <- offset + 1
	     binaryP[children[idx] + childOffset[children[idx]]] <- parent + offset
	     binaryP[parent + offset - 1] <- parent + offset
	     binaryY[children[idx] + childOffset[children[idx]]] <- Y[children[idx]]
	     binaryY[parent + offset - 1] <- Y[parent - 1]
	     childOffset[parent] <- offset
         }
      }   
  }
  Y[length(binaryY)] <- Y[length(Y)]
  binaryP[length(binaryP)] <- 0
  list(P=binaryP, Y=binaryY)
}
binaryP[length(binaryP)] <- 0
