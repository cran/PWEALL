wlrcom<-function(y=seq(1,5,by=1),d=rep(1,length(y)),z=c(1,0,1,0,1),tfix=max(y),p=c(1),q=c(1),eps=1.0e-08){
  n<-length(y)
  ymax<-max(y)
  tte<-sort(unique(y[d==1],fromLast=T))
  te<-tte[tte<=tfix]
  abc<-wlrutil(y=y,d=d,z=z,te=te,eps=eps)$mfunc
  ns<-length(p)
  cnames<-c("log-rank","Gehan","Tarone-Ware","Peto-Peto","mPeto-Peto")
  for (i in 1:ns){
      add<-paste0("FH","p=",p[i],",","q=",q[i])
      cnames<-c(cnames,add)
  }
  
  weights<-matrix(0,nrow=length(te),ncol=5+ns) 
  weights[,1]<-1
  weights[,2]<-abc[,6]
  weights[,3]<-sqrt(abc[,6])
  weights[,4]<-abc[,18]
  weights[,5]<-abc[,18]*abc[,6]/(abc[,6]+1/n)
  for (i in 1:ns){
    weights[,5+i]<-abc[,9]^p[i]*(1-abc[,9])^q[i]
    #weights[,7+i]<-abc[,9]^p[i]*(1-abc[,9])^q[i]
    #weights[,7+i]<-abc[,21]^p[i]*(1-abc[,21])^q[i]
  }
  def<-wlrcal(n=n,te=te,tfix=tfix,dd1=abc[,1],dd0=abc[,2],r1=abc[,4],r0=abc[,5],weights=weights)
  xpvalue<-matrix(2*(1-pnorm(abs(def$wlr))),nrow=1)
  xtest<-matrix(def$test,nrow=1)
  xvar<-matrix(def$var,nrow=1)
  xwlr<-matrix(def$wlr,nrow=1)
  colnames(xpvalue)<-colnames(xtest)<-colnames(xvar)<-colnames(xwlr)<-cnames
  list(n=n,test=xtest,var=xvar,wlr=xwlr,pvalue=xpvalue)
}
