BaumWelch.dthmm.Tree <- function (object, control = bwcontrol(), ...){
  makeSymmetric <- function(aa)
  {
    mm = ncol(aa)
    temp = expand.grid(1:mm, 1:mm)
    naInx = which(is.nan(aa[,1]))
    equiC =unlist(apply(temp ,1, function(x) which((x[1] == temp[,2]) & (x[2] ==temp[,1])) ))
    for(inx in naInx){
      if(!(equiC[inx] %in% naInx) ) aa[inx,] = aa[equiC[inx], ]
      else aa[inx, ] = 1/mm
    }
   (aa + aa[equiC,])/2 
  }
    x <- object$x
    Pi <- object$Pi
    delta <- object$delta
    distn <- object$distn
    pm <- object$pm
    P <- object$P
    tol <- control$tol
    numLeaf <- (length(P) + 1)/2 
   findBro <- function(P)
   {
     #Bro = list()
     Bro = P
     for(ii in seq(numLeaf+1, length(P))){
       inx = which(P == ii)
       stopifnot(length(inx) ==2)
       Bro[inx[1]] = inx[2]
       Bro[inx[2]] = inx[1]
       #for(tt in seq(1,length(inx))){
	   #Bro[[inx[tt]]] = inx[-tt] 
       #}
     }
     Bro
   }
   Bro = findBro(P)
    if (distn[1]!="glm"){
        Mstep <- parse(text=paste("Mstep.", distn,
                       "(x, cond, pm, object$pn)", sep=""))
    } else{
        Mstep <- parse(text=paste("Mstep.glm",
                  "(x, cond, pm, object$pn, distn[2], distn[3])", sep=""))
    }
    m <- ncol(Pi)
    n <- length(x)
    oldLL <- -Inf
    for (iter in 1:control$maxiter) {
	cond  <- Estep.Tree(x, P, Bro, numLeaf, Pi, delta, distn, pm, object$pn) 
        #cond <- Estep.Tree(x, Pi, delta, distn, pm, object$pn)
        diff <- cond$LL - oldLL
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(cond$LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        #if (diff < 0 & control$posdiff) stop("Worse log-likelihood on last iteration")
        if (eval(control$converge)) break
        #----  Mstep  ----
        Pi <- diag(1/apply(cond$v, MARGIN = 2, FUN = sum)) %*% 
            apply(cond$v, MARGIN = c(2, 3), FUN = sum)
	Pi = makeSymmetric(Pi)
	if(any(is.na(Pi))) browser()
        if (object$nonstat) delta <- cond$u[1:numLeaf, ]
        else delta <- compdelta(Pi)
        pm <- eval(Mstep)
	print(pm)
        oldLL <- cond$LL
    }
    object$delta <- delta
    object$Pi <- Pi
    object$u <- cond$u    
    object$v <- cond$v
    object$pm <- pm
    object$LL <- cond$LL
    object$iter <- iter
    object$diff <- diff
    return(object)
}
environment(Estep.Tree) = asNamespace("HiddenMarkov")


