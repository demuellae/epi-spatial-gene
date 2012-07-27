BicBin <- function(G, alpha, beta, p, proc_genes=TRUE) {
  M <- nrow(G) # TFs
  N <- ncol(G) # Genes

  max_scores <- c();

  # Selected TFs
  x <- rep.int(0, M);
  x[rnorm(M) < 0] <- 1

  # Selected Genes
  y <- rep.int(0, N); 
  y[rnorm(N) < 0] <- 1

  while (length(max_scores) < 2 || diff(max_scores)[length(max_scores)-1] != 0) {
    if (proc_genes) {
      n <- sum(y)

      s <- G %*% y # Rowsums over selected Genes (i.e., y == 1)
      ss <- sort(s, decr=TRUE, index.return=TRUE);

      C <- score(k <- cumsum(ss$x), m <- 1:M, n, p, alpha, beta)

      hits <- 1:which.max(C)
      x <- rep.int(0, M); x[ss$ix[hits]] <- 1 

      proc_genes <- FALSE;
    } else {
      m <- sum(x)

      s <- as.matrix(t(x) %*% G) # Colsums over selected TFs (i.e., x == 1)
      ss <- sort(s, decr=TRUE, index.return=TRUE);

      C <- score(k <- cumsum(ss$x), m, n <- 1:N, p, alpha, beta)

      hits <- 1:which.max(C)
      y <- rep.int(0, N); y[ss$ix[hits]] <- 1 

      proc_genes <- TRUE;
    }

    if(max(C) >= max(max_scores)){
      outx <- x
      outy  <- y
    } 
    max_scores <- c(max_scores, max(C))
  }

  max_score  <- max(max_scores)
  out <- list()
  out$score  <- max_score
  out$x <- outx
  out$y  <- outy
  return(out)
}

score <- function(k, m, n, p, alpha, beta) {
  mnp <- m*n*p
  res1 <- (k - mnp)^2 / (m^alpha * n^beta * (k+mnp))
  res2 <- (k - mnp)^2 / (3 * m^(1+alpha) * n^(1+beta) * p)

  hits <- k-(2*mnp) > 0
  res <- c(res1[hits], res2[!hits])

  # Same as: if (k - mnp > mnp) { res1 } else { res2 }
}

###################### Example #########################

BCBIN <- function(mat, numRep=300, alpha=0.5, beta=0.5, number=25)
{

  require(SparseM)
  require(biclust)
  mat <- as.matrix.csr(mat)
  M <- nrow(mat) # TFs
  N <- ncol(mat) # Genes
  clust = new("Biclust")
  Parameters = list()
  Parameters$Call = sprintf( "BCBIN(mat, numRep=%d, alpha=%f, beta=%f, number=%d)", numRep,
			    alpha, beta, number) 
  Parameters$Method = "BICBIN"   
  clust@Parameters = Parameters 
  clust@RowxNumber = matrix(FALSE, M, number) 
  clust@NumberxCol = matrix(FALSE, number, N) 
  clust@Number = number
  for(num in 1:number){

    p <- sum(attr(mat, "ra")) / (M*N)

    max <- list();
    max$score <- 0;

    for (i in 1:numRep) {
      print(i)
      res <- BicBin(mat, alpha, beta, p, proc_genes=TRUE)
      if (res$score > max$score) { max <- res }
      res <- BicBin(mat, alpha, beta, p, proc_genes=FALSE)
      if (res$score > max$score) { max <- res }
    }
    clust@RowxNumber[,num]= res$x==1
    clust@NumberxCol[num, ]= res$y==1

    mat[clust@RowxNumber[,num],clust@NumberxCol[num, ] ] <- 0
  }
  return(clust)
}
