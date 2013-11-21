# generate tree
#t <- function(p)
  #list(if(runif(1) < p) t(p) else 0, if(runif(1) < p) t(p) else 0)
numLeaf <- 100
n = 2*numLeaf -1 
h <- ceiling(log(numLeaf,base=2) +1)
P  <-  rep(0,n-1)
Pi[1,] = c(0,1e-3, 1e-6)
Pi[2,]=c(.1, 0, .05)
Pi[3,] = c(1e-6,1e-3,0)
diag(Pi) = 1 - rowSums(Pi)
Pi = matrix(1/3,9,3) 

curr.nodes <-  seq(numLeaf); curr.parent <- numLeaf 
for(jj in seq(numLeaf - 1)){
  curr.child  <- sample(curr.nodes, replace=F, size=2)
  curr.parent  <- curr.parent + 1
  P[curr.child] <- curr.parent 
  curr.nodes <- curr.nodes[!(curr.nodes %in% curr.child)]
  curr.nodes <- c(curr.nodes, curr.parent)
}
B  <- findBro(P)
#generate z
z <- rep(0, n)
z[n] <- 0 # initial state 
for(parent in seq(n,numLeaf+1, -1)){
  curr.child <- which(P==parent)
  z[curr.child] <- sample(-1:1, size=length(curr.child), prob=Pi[z[parent]+2,])
}
# generate Y
Y <- rep(0, n)
shape <- matrix(0, 3,2)
shape[1,]= c(.2, .9)
shape[2,] = c(.2, .5)
shape[3,] = c(.9, .5)
for(node in seq(n)){
  Y[node]  <- rbeta(1,shape1=shape[z[node]+2,1], shape2= shape[z[node]+2,2])
}

write.table(file="P",x = c(P, 0) , row.names = F, col.names =F,  sep=",", quote=F )
write.table(file="z",x = z, row.names = F, col.names =F,  sep=",", quote=F )
write.table(file="Y",x = Y, row.names = F, col.names =F,  sep=",", quote=F )


# tissue-tree infomation
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

