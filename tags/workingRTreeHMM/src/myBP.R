# M is number of states
# T is time (sequence length)
#uu(M \x T) 2 dimensional matrix = p(c(t) =j/X(1:T) ) (4.13)
#vv(M^2 \x T-1) = p(c(t-1)=j,c(t)=k/X(1:T) (4.14) 
#GammaC transition probability state(t-1) x state(t)
#emission(k \x x): emission probability state x observation
# alpha forward probability M X T
# beta backward probability

# calculate Alpha
# calculate Beta
beta.HMM.forwardBackward = function(XX, MM, TT, shape, GammaC, delta=NULL)
{
  if (is.null (delta)) delta = solve(t (diag(MM) - GammaC +1) , rep(1 , m ) )
  #shape = data.frame(rbind(shape1, shape2))
  allProbs = outer( XX, shape, function(s,t) dbeta(s , t[[1]][1], t[[1]][2]))
  Lalpha = Lbeta = martix(NA, MM, TT)
  temp =   delta*allProbs[1,]
  sumtemp = sum(temp)
  ll = log(sumtemp) 
  temp = temp/sumtemp
  Lalpha[,1] = log(temp) + ll
  for(tt in seq(2,TT)){
    temp = temp %*% GammaC * allProbs[tt,] 
    sumtemp = sum(temp)
    ll = ll + log(sumtemp)
    temp = temp/sumtemp
    Lalpha[,tt] = log(temp) + ll
  }
  Lbeta[,TT] = rep(0,MM)
  temp = rep(1/MM, MM)
  ll = log(MM)
  for(tt in seq(TT-1,1)){
    temp = GammaC %*% (allProbs[tt+1,]*temp) 
    Lbeta[,tt] = ll + log(temp)
    sumtemp = sum(temp)
    temp = temp/sumtemp
    ll = log(temp) + ll
  }
  list(Lalpha=Lalpha, Lbeta=Lbeta)
}


betaF <- function(lshape, XX, uuJ )
{
  shape = exp(lshape)
  print(shape)
  print(sum(dbeta( XX,shape1=shape[1], shape2=shape[2], log=T)))
  #-mean(log(dnorm(inp, mean=param[1], sd=param[2])))
  -mean(uuJ* dbeta( XX,shape1=shape[1], shape2=shape[2], log=T))
}

pois.HMM.EM <-function(XX,MM, TT, shape,GammaC,delta,maxiter=1000,tol=1e-6,...)
{
  #lambda.next = lambda
  shape.next = shape
  #shape2.next = shape$shape2
  GammaC.next = GammaC
  delta.next = delta
  for ( iter in 1: maxiter )
  {
    LallProbs = outer( XX, shape, function(s,t) dbeta(s , t[[1]][1], t[[1]][2], log=T))
    #luuT = LalphaT +  LbetaT - logProbT
    #lvvT = LalphaT* LGammaCC * t(emission[,XXT] * LbetaT)/ logProbT
    #lvvT = outer(LalphaT, LbetaT+LallprobT, FUN="+" ) + LGammaCC
    # This is possibility but it creates underflow problem
    #vv = (exp(Lalpha[,1:(TT-1)]) %*% exp(t(Lbeta[,2:TT] + t(Lallprobs[2:TT,])))) * GammaCC  
    lalphaBeta = beta.HMM.forwardBackward(XX, MM, TT, shape, GammaC, delta=NULL)
    Lalpha = lalphaBeta$Lalpha
    Lbeta = lalphaBeta$Lbeta
    c=max(Lalpha[,TT])
    llk=c+log(sum(exp(Lalpha[,TT]-c)))
    for(jj in 1:MM){
      for(kk in 1:MM){
	GammaC.next[jj,kk] = GammaC[jj,kk] * sum(exp(Lalpha[jj, 1:(TT-1)] + LallProbs[2:TT, k] + Lbeta[k, 2:TT]- llk))
      }
      #lambda.next[j]=sum(exp(Lalpha[j,]+Lbeta[j,]-llk)*x)/sum(exp(Lalpha[j,]+Lbeta[j,]-llk))
      shape0 = log(unlist(shape[,jj]))
      uuJ = exp(Lalpha[jj,]+Lbeta[jj,]-llk)
      temp = nlm(betaF, shape0, XX, uuJ, iterlim=1000)
      shape.next[,jj] = temp$uJJ
    }
    GammaC.next=GammaC.next/apply(GammaC.next,1,sum)
    delta.next=exp(Lalpha[,1]+Lbeta[,1]-llk)
    delta.next=delta.next/sum(delta.next)
    crit = sum(abs(shape-shape.next))+
    sum(abs(GammaC-GammaC.next))+
    sum(abs(delta-delta.next))
    if(crit<tol)
    {
      np=m*m+m-1
      AIC=-2*(llk-np)
      BIC=-2*llk+np*log(n)
      return(list(lambda=lambda,GammaC=GammaC,delta=delta,
		  mllk=-llk,AIC=AIC,BIC=BIC))
    }
    lambda=lambda.next
    GammaC=GammaC.next
    delta=delta.next
  }
  print(paste("Noconvergenceafter",maxiter,"iterations"))
  NA
}
f <- function(lshape, XX, uuJ )
{
  shape = exp(lshape)
  print(shape)
  print(sum(dbeta( XX,shape1=shape[1], shape2=shape[2], log=T)))
  #-mean(log(dnorm(inp, mean=param[1], sd=param[2])))
  -mean(uuJ* dbeta( XX,shape1=shape[1], shape2=shape[2], log=T))
}


func1 <- function(lshape, inp )
{
   shape = exp(lshape)
  -mean(dbeta(inp, shape1=shape[1], shape2=shape[2], log=T))
}
inp = rbeta(1000, .2, .9)
lshape = log(c(.2,.9))
out = nlm(func1, lshape, inp=inp, iterlim=1000)
exp(out$estimate)

func = function(param, inp)
{
  -mean(log(dnorm(inp, mean=param[1], sd=param[2])))
}

inp = rnorm(100, mean=10, sd=1)
param0 = c(9,10)
out = nlm(func, param0, inp=inp, iterlim=1000)
