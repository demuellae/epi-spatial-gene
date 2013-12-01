"Estep.Tree" <-
    function (x, P, Bro, numLeaf, Pi, delta, distn, pm, pn = NULL) 
    {
	threeD2twoD <- function(aa)
	{
	    mm = dim(aa)[2] 
	    aa = as.data.frame(aa)
	    colN =  paste("t", 1:mm, sep="")
	    colnames(aa)= colN
	    temp = expand.grid(1:mm, 1:mm) 
	    #aa$Var1 = temp$Var1
	    #aa$Var2 = temp$Var2
	    temp1 = (aggregate(aa, list(Var1 = temp$Var1), FUN=sum)[,colN]   + aggregate(aa, by= list(Var2 = temp$Var2), FUN=sum)[,colN])/6
	    as.matrix(temp1) 
	}
	dfunc <- makedensity(distn)
	m <- ncol(Pi)
	n <- length(x)
	#PiF  <- Pi
	#PiB  <- threeD2twoD(Pi)
	y <- forwardback.Tree(x, numLeaf, Pi, P, delta, distn, pm, pn, fortran=F)
	logbeta <- y$logbeta
	logalpha2 <- y$logalpha2
	logalpha  <- y$logalpha
	LL <- y$LL
	u <- exp(logalpha + logbeta - LL)
	#browser()
	if(!isTRUE(all.equal(rowSums(u), rep(1, dim(u)[1]))))
	  browser()

	stopifnot(all.equal(rowSums(u), rep(1, dim(u)[1])))
	v <- array(NA, dim = c(n - 1, m, m)) #absent for root , now defined for child or each edge of tree
	for (k in 1:m) {
	    logprob <- dfunc(x=x[-1:-numLeaf], getj(pm, k),
			     getj(pn, -1), log=TRUE)
	    logPi <- matrix(log(Pi[, k]), byrow = TRUE, nrow = n - 
			    1, ncol = m)
	    logPbeta <- matrix(logprob + logbeta[-1:-numLeaf, k],
			       byrow = FALSE, nrow = n - numLeaf, ncol = m) # there is only nonLeaf, therefore -numLeaf for index
	    v[, , k] <- logPi + logalpha[ -n, ] + logPbeta[(P[-n] - numLeaf), ] - LL
	}
	v <- exp(v)
	if(!isTRUE(all.equal(apply(v,1, sum), rep(1, dim(v)[1]))))
	{
	  cat("all v not equal to 1 \n")
	  browser()
	}
	return(list(u = u, v = v, LL = LL))
    }
environment(Estep.Tree) = asNamespace("HiddenMarkov")
