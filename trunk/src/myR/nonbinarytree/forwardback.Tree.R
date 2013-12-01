"forwardback.Tree" <-
    function(x, numLeaf, Pi, P, delta, distn, pm, pn = NULL, fortran=F){
	#parent: parent array of the trees
	#child : child array of the trees
	myOuter <- function(aa, bb)
	{
	    bbLen = 1:length(bb)
	    if(any(is.na(aa[bbLen]))) aa[bbLen] = bb
	    else aa = unlist(as.list(outer(aa[bbLen], bb)))
	    aa
	}

	m <- ncol(Pi)
	n <- length(x)
	dfunc <- makedensity(distn)
	prob <- matrix(as.double(0), nrow=n, ncol=m)
	for (k in 1:m)
	    prob[,k] <- dfunc(x=x, getj(pm, k), pn, log=FALSE)
	#   forward probabilities alpha_ij
	# contains (phi(l) outer phi(r))
	phiT <- numeric(m)
	phi2 <- matrix(1, nrow=n, ncol=m) # initialize to na, contains phi(l)*phi(r) 
	phi <- matrix(as.double(rep(0, m*n)), nrow=n) # contains simple phi(t)
	logalpha <- matrix(as.double(rep(0, m*n)), nrow=n)
	# Dan logalpha2 is no more required. Safely ignore it
	logalpha2<- matrix(as.double(NA), nrow=n, ncol=m^2) # contains (alpha(l) out alpha(r))
	#print(prob)
	lscale2 <- double(n) # contains log(w_r* w_l)
	lscale <- double(n) # contains log(w_t)
	if (fortran!=TRUE){
	    #  loop1 using R code
	    for (i in (1:n)){
		#if(i>0) browser()
		if (i <= numLeaf) 
		  phiT[]  <-  delta[i,] * prob[i,]
		else phiT[] <- phi2[i,] * prob[i,]
		#phiT <- phiT*prob[i,] # prob[i,] is 1xM vector
		sumphi <- sum(phiT) # phi(l) out phi(r) B 1'
		phi[i,] <-  phiT/sumphi # this is phi(t)
		lscale[i] <- lscale2[i] + log(sumphi)
		#print(lscaleT) 
		#if(lscale[i] < -1000) browser()
		logalpha[i,] <- log(phi[i,]) + lscale[i]
		if(i < n) {
		    # contains prod (phi gamma)
		    phi2[P[i],] <- phi2[P[i],] *  (phi[i,,drop=F] %*% Pi)
		    lscale2[P[i]] <- lscale2[P[i]] + lscale[i] 
		    logalpha2[P[i],] <- phi2[P[i],] + lscale2[P[i]]
		}
	    }
	    LL <- lscale[n]
	} else{
	    if (!is.double(Pi)) stop("Pi is not double precision")
	    if (!is.double(prob)) stop("prob is not double precision")
	    memory0 <- rep(as.double(0), m)
	    loop1 <- .Fortran("loop1", m, n, phi, prob, Pi, logalpha,
			      lscale, memory0, PACKAGE="HiddenMarkov")
	    logalpha <- loop1[[6]]
	    LL <- loop1[[7]]
	}
	#print(logalpha)
	#   backward probabilities beta_ij
	#browser()
	logbeta <- matrix(as.double(rep(0, m*n)), nrow=n)
	theta  <-  matrix(as.double(NA), nrow=n, ncol=m)
	lscaleB <- as.double(rep(0, n))
	theta[n,] <- as.double(rep(1/m, m))
	lscaleB[n] <- as.double(log(m))
	temp <- numeric(3)
	thetaT <- matrix(0, nrow=m, ncol=1)
	if (fortran!=TRUE){
	    #  loop2 using R code
	    for (i in seq(n-1, 1, -1)){
		temp[] <- phi2[P[i],] * prob[P[i],]*theta[P[i],]/ (phi[i,] %*% Pi) # all ops component wise
		thetaT[] <- Pi %*% temp
		#dim(thetaT)  <- c(m,m)
		sumtheta <- sum(thetaT)
		theta[i,] <- thetaT/sumtheta
		lscaleB[i] <- lscaleB[P[i]] + (lscale2[P[i]] - lscale[i])+ log(sumtheta)
		logbeta[i,] <- log(theta[i,]) + lscaleB[i] 
	    }
	} else{
	    memory0 <- rep(as.double(0), m)
	    loop2 <- .Fortran("loop2", m, n, phi, prob, Pi, logbeta,
			      lscale, memory0, PACKAGE="HiddenMarkov")
	    logbeta <- loop2[[6]]
	}
	stopifnot(all.equal(log(rowSums(exp(logbeta + logalpha))) ,rep(LL, dim(logbeta)[1])))
	#browser()
	return(list(logalpha=logalpha, logbeta=logbeta, LL=LL, logalpha2=logalpha2))
    }
environment(forwardback.Tree) = asNamespace("HiddenMarkov")
