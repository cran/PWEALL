hxbeta<-function(x=c(0.5,1),y=seq(.1,1,by=0.01),d=rep(1,length(y)),tfix=2,K=20,eps=1.0e-06){
  ###x: output time points
  ###y: observed times
  ###d: non-censoring indicators
  ###tfix: maximum timepoint where hazard is estimated
  ###K: tuning parameter, the smaller the smoother of the estimate
  nx=length(x)
  ny=length(y)
  lambda<-rep(0,nx)
  xs=x/tfix
  te<-sort(unique(y[d==1&y<tfix],fromLast=T))
  nte<-length(te)
  de<-re<-rep(0,nte)
  for(i in 1:nte){
    de[i]<-sum(d==1&abs(y-te[i])<eps)
    re[i]<-sum(y>=te[i])
  }
  de=de/ny;re=re/ny
  ys=te/tfix  
  atemp=floor(xs*K) ###do we have to use integer?
  #atemp=xs*K
  ax=1+atemp
  bx=K-atemp
  btemp=de/re
  for (i in 1:nx){
    ctemp=dbeta(ys,shape1=ax[i],shape2=bx[i])*btemp
    lambda[i]<-sum(ctemp)
  }
  lambda=lambda/tfix
  list(lambda=lambda)
}