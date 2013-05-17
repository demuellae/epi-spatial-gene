"forwardback.Tree" <-
function(x, PiF, PiB, P, delta, distn, pm, pn = NULL, fortran=TRUE){
    #parent: parent array of the trees
    #child : child array of the trees
  myOuter <- function(aa, bb)
  {
    bbLen = 1:length(bb)
    if(is.na(aa[bbLen])) aa[bbLen] = bb
    else aa = unlist(as.list(outer(aa[bbLen], bb)))
    aa
  }

    m <- nrow(Pi)
    n <- length(x)
    dfunc <- makedensity(distn)
    prob <- matrix(as.double(0), nrow=n, ncol=m)
    for (k in 1:m)
        prob[,k] <- dfunc(x=x, getj(pm, k), pn, log=FALSE)
    #   forward probabilities alpha_ij
    phi <- matrix(as.double(NA), m, m^2) # initialize to na. 
    logalpha <- matrix(as.double(rep(0, m*n)), nrow=n)
    lscale <- as.double(rep(0, length(P)))
    if (fortran!=TRUE){
        #  loop1 using R code
        for (i in (1:length(P))){
            if (i > numLeaf)  phiT <- phi[i,] %*% PiF
	    else phiT  <-  delta 
            phiT <- phiT*prob[i,]
            sumphiT <- sum(phiT)
	    phiT <-  phiT/sumphiT
            lscaleT <- lscale[i] + log(sumphi) 
            logalpha[i,] <- log(phiT) + lscaleT
	    if(i < length(P)) {
	      phi[P[i],] <- myOuter(phi[P[i],], phiT)
	      lscale[P[i]] <- lscale[,P[i]] + lscaleT 
	    }
        }
        LL <- lscaleT
    } else{
        if (!is.double(PiF)) stop("PiF is not double precision")
        if (!is.double(prob)) stop("prob is not double precision")
        memory0 <- rep(as.double(0), m)
        loop1 <- .Fortran("loop1", m, n, phi, prob, PiF, logalpha,
                          lscale, memory0, PACKAGE="HiddenMarkov")
        logalpha <- loop1[[6]]
        LL <- loop1[[7]]
    }
    #   backward probabilities beta_ij
    logbeta <- matrix(as.double(rep(0, m*n)), nrow=n)
    phi <- as.double(rep(1/m, m))
    lscale <- as.double(log(m))
    if (fortran!=TRUE){
        #  loop2 using R code
        for (i in seq(n-1, 1, -1)){
            phi <- PiB %*% (prob[i+1,]*phi)
            logbeta[i,] <- log(phi) + lscale
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
        }
    } else{
        memory0 <- rep(as.double(0), m)
        loop2 <- .Fortran("loop2", m, n, phi, prob, PiB, logbeta,
                          lscale, memory0, PACKAGE="HiddenMarkov")
        logbeta <- loop2[[6]]
    }
    return(list(logalpha=logalpha, logbeta=logbeta, LL=LL))
}

