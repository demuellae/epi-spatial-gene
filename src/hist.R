histAnnotation <- function(mat) {
out <- apply(mat, 1, numMedium)
return(out)  
}

numStrong <- function(mat){
out <- sum(mat[mat==3])/3
return(out)
}
numMedium <- function(mat){
out <- length(mat[mat > 1])
return(out)
}
numWeak <- function(mat){
out <- length(mat[mat > 0])
return(out)
}
