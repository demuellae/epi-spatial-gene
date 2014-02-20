"forwardback.Tree" <-
    function(x, numLeaf, PiF, PiB, P, Bro, delta, distn, pm, pn = NULL, fortran=F){
	#parent: parent array of the trees
	#child : child array of the trees
	myOuter <- function(aa, bb)
	{
	    bbLen = 1:length(bb)
	    if(any(is.na(aa[bbLen]))) aa[bbLen] = bb
	    else aa = unlist(as.list(outer(aa[bbLen], bb)))
	    aa
	}
browser()

	m <- ncol(Pi)
	n <- length(x)
	dfunc <- makedensity(distn)
	prob <- matrix(as.double(0), nrow=n, ncol=m)
	for (k in 1:m)
	    prob[,k] <- dfunc(x=x, getj(pm, k), pn, log=FALSE)
	#   forward probabilities alpha_ij
	# contains (phi(l) outer phi(r))
	phi2 <- matrix(as.double(NA), nrow=n, ncol=m^2) # initialize to na, contains phi(l)*phi(r) 
	phi <- matrix(as.double(rep(0, m*n)), nrow=n) # contains simple phi(t)
	logalpha <- matrix(as.double(rep(0, m*n)), nrow=n)
	logalpha2<- matrix(as.double(NA), nrow=n, ncol=m^2) # contains (alpha(l) out alpha(r))
	#print(prob)
	lscale2 <- as.double(rep(0, length(P))) # contains log(w_r* w_l)
	lscale <- as.double(rep(0, length(P))) # contains log(w_t)
	if (fortran!=TRUE){
	    #  loop1 using R code
	    for (i in (1:n)){
		#if(i==1) browser()
		if (i > numLeaf)  phiT <- phi2[i,] %*% PiF
		else phiT  <-  delta[i,] 
		phiT <- phiT*prob[i,] # prob[i,] is 1xM vector
		sumphi <- sum(phiT) # phi(l) out phi(r) B 1'
		phi[i,] <-  phiT/sumphi # this is phi(t)
		lscale[i] <- lscale2[i] + log(sumphi)
		#print(lscaleT) 
		logalpha[i,] <- log(unlist(phi[i,])) + lscale[i]
		if(i < n) {
		    phi2[P[i],] <- myOuter(phi2[P[i],], phi[i,])
		    lscale2[P[i]] <- lscale2[P[i]] + lscale[i] 
		    logalpha2[P[i],] <- log(phi2[P[i],]) + lscale2[P[i]]
		}
	    }
	    LL <- lscale[n]
	} else{
	    if (!is.double(PiF)) stop("PiF is not double precision")
	    if (!is.double(prob)) stop("prob is not double precision")
	    memory0 <- rep(as.double(0), m)
	    loop1 <- .Fortran("loop1", m, n, phi, prob, PiF, logalpha,
			      lscale, memory0, PACKAGE="HiddenMarkov")
	    logalpha <- loop1[[6]]
	    LL <- loop1[[7]]
	}
	#print(logalpha)
	#   backward probabilities beta_ij
	#browser()
	logbeta <- matrix(as.double(rep(0, m*n)), nrow=n)
	theta  <-  matrix(as.double(NA), nrow=n, ncol=m)
	lscaleB <- as.double(rep(0, length(P)))
	theta[n,] <- as.double(rep(1/m, m))
	lscaleB[n] <- as.double(log(m))
	if (fortran!=TRUE){
	    #  loop2 using R code
	    for (i in seq(n-1, 1, -1)){
		thetaT <- PiF %*% (prob[P[i],]*theta[P[i],])
		dim(thetaT)  <- c(m,m)
		thetaT = phi[Bro[i],] %*% thetaT
		sumtheta <- sum(thetaT)
		theta[i,] <- thetaT/sumtheta
		lscaleB[i] <- lscaleB[P[i]] + lscale[Bro[i]] + log(sumtheta)
		logbeta[i,] <- log(theta[i,]) + lscaleB[i] 
	    }
	} else{
	    memory0 <- rep(as.double(0), m)
	    loop2 <- .Fortran("loop2", m, n, phi, prob, PiB, logbeta,
			      lscale, memory0, PACKAGE="HiddenMarkov")
	    logbeta <- loop2[[6]]
	}
	browser()
	#stopifnot(all.equal(log(rowSums(exp(logbeta + logalpha))) ,rep(LL, dim(logbeta)[1])))
	return(list(logalpha=logalpha, logbeta=logbeta, LL=LL, logalpha2=logalpha2))
    }
environment(forwardback.Tree) = asNamespace("HiddenMarkov")
