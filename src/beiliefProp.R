# Copy of R-code given in the Appendix 
# A.3 HMM with bivariate normal state-dependent distributions

bivnorm.HMM.pn2pw =  function(mu,sigma,gamma,m)
{
  tsigma = log(sigma)
  #tcorr=log((1+corr)/(1-corr))
  tgamma = NULL
  if(m>1)
  {
    foo    = log(gamma/diag(gamma))
    tgamma = as.vector(foo[!diag(m)])
  }
  parvect = c(as.vector(mu),as.vector(tsigma),
	      tgamma)
  parvect
}

bivnorm.HMM.pw2pn = function(parvect,m)
{
  mu    = matrix(parvect[1:m],m,1)
  sigma = matrix(exp(parvect[(m+1):(2*m)]),m,1)
  #temp=exp(parvect[(4*m+1):(5*m)])
  #corr=(temp-1)/(temp+1)
  gamma = diag(m)
  if(m>1)
  {
    gamma[!gamma]=exp(parvect[(2*m+1):length(parvect)])
    gamma=gamma/apply(gamma,1,sum)
  }
  print(gamma)
 # this was initialization with stationary distribution 
  #delta=solve(t(diag(m)-gamma+1),rep(1,m))
# initializing state1 is initial state.
  delta = c(1,rep(0, m-1))
  list(mu=mu,sigma=sigma,gamma=gamma,
       delta=delta)
}


bivnorm.HMM.mllk=function(parvect,x,m,...)
{
  require(mvtnorm)
  n        = dim(x)[1]
  p        = bivnorm.HMM.pw2pn(parvect,m)
  foo      = p$delta
  #covs     = array(NA,c(2,2,m))
  #for (j in 1:m)                             
  #{
    #covs[,,j] = diag(p$sigma[j,])%*%
    #matrix(c(1, p$corr[j],p$corr[j],1),2,2)%*%
    #diag(p$sigma[j,])
  #}
  covs = (p$sigma)^2
  P        = rep(NA,m)
  lscale   = 0                                
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      P[j] = dnorm(x=x[i],  
		     mean=p$mu[j], sd=covs[j])
      #print(P[j])
    }
    foo    = foo%*%p$gamma*P
    sumfoo = sum(foo)
    lscale = lscale+log(sumfoo)
    foo    = foo/sumfoo
  }
  mllk     = -lscale                    
  mllk
}


bivnorm.HMM.mle=function(x,m,mu0,sigma0,gamma0,...)
{
  n     = dim(x)[1]
  start = bivnorm.HMM.pn2pw(mu0,sigma0,gamma0,m)
  print(start)
  mod   = nlm(bivnorm.HMM.mllk,p=start,x=x,m=m,
	      steptol = 1e-4,iterlim = 10000)
  mllk  = mod$minimum
  code  = mod$code
  p     = bivnorm.HMM.pw2pn(mod$estimate,m)
  np    = m*(m+4)
  AIC   = 2*(mllk+np)
  BIC   = 2*mllk+np*log(n)
  list(mu=p$mu,sigma=p$sigma,
       gamma=p$gamma,delta=p$delta,code=code,
       mllk=mllk,AIC=AIC,BIC=BIC)
}

bedF = read.table(file="DGFHCMSig.chr1.bed", header = F, sep="\t", stringsAsFactors =F,  strip.white =T)
chrmLen = tail(bedF, 1)[3]
bedF$zerostrt = c(0, bedF$V3[1:length(bedF$V3) -1])
dgf = unlist(apply(bedF[,-1], 1, function(x)  {c(rep(0,x[1]-x[4]), rep(x[3], x[2] - x[1]))}))
names(dgf) = NULL
#Kernsmooth for smoothing the data
wind = wind_= 10
sigmaS = sigmaS_ = 1
sigmaT = sigma_T = sigmaN
distS = exp(-((-wind:wind)/sigmaS)^2 /2 )
#performs bilinear filtering of 1D data.
dgf.temp = temp
dgf.smooth = unlist(lapply(seq(wind+1, length(dgf.temp)-wind) , function(x) {
	temp = dgf.temp[(x-wind):(x+wind)]
	distI = exp(-((temp - dgf.temp[x])/sigmaT)^2 / 2)
	wtemp = distI * distS
	wtemp = sum(wtemp * temp)/sum(wtemp)
	}) )

#sample Rcpp code
library(inline)
library(Rcpp)
bilateralCode = "
  int wind = as<int>(wind_);
  NumericVector dgf(dgf_);
  // Automatic  initialization to zero
  int n=dgf.size();
  NumericVector out(n - (2*wind));
  NumericVector distS(2*wind +1);
  double temp;
  double sigmaS = as<double>(sigmaS_);
  double sigmaT = as<double>(sigmaT_);
  double sumWeight;

  for (int i = 0; i < 2*wind+1; i++) {
    distS[i] = exp(- pow((i-wind)/sigmaS,2) /2);
    //std::cout<<distS[i] <<std::endl;
  }
  for (int ii = wind; ii < n-wind; ii++) {
    out[ii-wind]=0;
   sumWeight = 0; 
    for (int jj = 0; jj <= 2*wind; jj++) {
      temp = distS[jj] * exp(- pow((dgf[jj+ ii - wind] - dgf[ii])/sigmaT, 2)/2);
      //std::cout << temp << std::endl;
      out[ii-wind] += dgf[jj+ ii - wind] * temp;
      sumWeight  += temp;  
    }
    out[ii-wind] =  out[ii-wind]/sumWeight;
  }
  
  return out;
"
bilateral = cxxfunction(signature(wind_="integer",dgf_="numeric", sigmaS_="numeric", sigmaT_="numeric"), bilateralCode, plugin="Rcpp", includes="#include <math.h>
			#include <iostream>")
dgf_ = dgf
wind_=10
sigmaS_=1
sigmaN = sd(dgf[dgf !=0])
sigmaT_=sigmaN
out = bilateral(dgf_=dgf_, wind_= wind_, sigmaS_=sigmaS_, sigmaT_=sigmaT_)
#testing
temp = dgf[20959546:20961533]
out = bilateral(dgf_=temp, wind_= wind_, sigmaS_=sigmaS_, sigmaT_=sigmaT_)
jpeg("doc/bilateral.jpg")
plot(temp[200:500], type="l", xlab="DNASe across region near plink gene chr1:20,959,746-20,960,046", ylab="DNase Signal", title="Bilateral filter with sigmaS = 1, sigmaT= SigmaN")
lines(out[190:490], type="l", col=2)
legend(1, 35, c("DNASE","DNASE SMOOTH"), cex=0.8,  col=c("blue","red"),lty=1);
dev.off()
dgf_ = dgf[10000000:30000000]
out = bilateral(dgf_=dgf_, wind_= wind_, sigmaS_=sigmaS_, sigmaT_=sigmaT_)

#subsetting=seq( 149240568, (149240568+2492405))
dgf.smooth = out[10000000:30000000]

dgf.diff = dgf.smooth[-1] - dgf.smooth[-length(dgf.smooth)]
dgf.diff.pos = dgf.diff[dgf.diff > 0]
dgf.diff.neg = dgf.diff[dgf.diff < 0]
mu0.pos = quantile(dgf.diff.neg, probs=0.25)
mu0.neg = quantile(dgf.diff.pos, probs=0.75)
mu0.pos.extrem = quantile(dgf.diff.neg, probs=0.95)
mu0.neg.extrem = quantile(dgf.diff.pos, probs=0.95)
mu0 = c(0, mu0.pos, mu0.neg, 0)
sd0.pos = sd( dgf.diff.pos[ dgf.diff.pos > mu0.pos] ) 
sd0.neg = sd( dgf.diff.neg[ dgf.diff.neg < mu0.neg] )
sd0.bind = sd(c(dgf.diff.pos[ dgf.diff.pos >  mu0.pos.extrem], dgf.diff.neg[ dgf.diff.neg <  mu0.neg.extrem]))
sigma0 = c(sd(dgf.diff), sd0.pos, sd0.neg, sd0.bind) 
sm=1e-4
bg=.05
gammaFrom1 = c(.005, sm, sm)
gammaFrom1 = c(1-sum(gammaFrom1), gammaFrom1)
gammaFrom2 = c(sm, .01, sm)
gammaFrom2 = c(gammaFrom2[1],  1-sum(gammaFrom2), gammaFrom2[c(2,3)])
gammaFrom3 = c(sm, sm,20/21)
gammaFrom3 = c(gammaFrom3,  1-sum(gammaFrom3))
gammaFrom4 = c(1/20,1/20,.25/20,17.75/20)
#gammaij: i is starting state and j is ending state
gamma0 = rbind(gammaFrom1, gammaFrom2, gammaFrom3, gammaFrom4)

hmmOut = bivnorm.HMM.mle(x=as.matrix(dgf.diff[1:1000]), m=4, mu0=mu0,sigma0=sigma0,gamma0=gamma0)


dgf.pos = dgf.smooth[dgf.smooth >= 1]
mu0.open = mean(dgf.pos)
sd0.open = sd(dgf.pos)
dgf.closed = c(dgf.smooth[dgf.smooth == 0], dgf.smooth[dgf.smooth < quantile(dgf.pos, probs=0.01)])
mu0.closed = mean(dgf.closed) 
sd0.closed = sd(dgf.closed)
temp = dgf.smooth[dgf.smooth < quantile(dgf.pos, probs=0.05)]
dgf.bound = c(rep(length(temp),1), temp )
mu0.bound = mean(dgf.bound) 
sd0.bound = sd(dgf.bound)
mu0 = c(mu0.closed, mu0.open, mu0.bound)
sigma0 = c(sd0.closed, sd0.open, sd0.bound)
gamma0 = rbind(c(99/100, 1/100, 1e-14),
	       c(.5/40, 39/40, .5/40),
	       c(1e-14, 1/21,20/21))
hmmOut = bivnorm.HMM.mle(x=as.matrix(dgf.smooth), m=3, mu0=mu0,sigma0=sigma0,gamma0=gamma0)
hmmOut = bivnorm.HMM.mle(x=as.matrix(dgf.smooth[10920000:109900000]), m=3, mu0=mu0,sigma0=sigma0,gamma0=gamma0)

