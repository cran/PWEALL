############################################################################################################################
#   A function to simulate the power of LR and Cox
#     account for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
pwesim<-function(t=seq(1,2,by=0.1),taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
                     n=1000,rn=200,testtype=c(1,2,3,4)){
  ##t: the time point where the overall log hazard ratio beta is calculated
  ##pi1: proportion of treatment group
  ##rate11: hazard before dilution for the treatment group
  ##rate21: hazard after dilution for the treatment group
  ##rate31: hazard for treatment discontinuation for the treatment group
  ##ratec1: hazard for loss to follow-up for the treatment group
  ##rate10: hazard before dilution for the control group
  ##rate20: hazard after dilution for the control group
  ##rate30: hazard for treatment discontinuation for the control group
  ##ratec0: hazard for loss to follow-up for the control group
  ##tchange: points at which hazard changes
  ##rp21, rp20 re-randomization prob for the tx and contl groups

  nt<-length(t)
  ntest<-length(testtype)
  outr<-array(0,dim=c(rn,nt,ntest))
  for (r in 1:rn){
    E<-T<-C<-Z<-delta<-rep(0,n)
    E<-rpwu(nr=n,u=u,ut=ut)$r
    Z<-rbinom(n,1,pi1)
    n1<-sum(Z)
    n0<-sum(1-Z)
    C[Z==1]<-rpwe(nr=n1,rate=ratec1,tchange=tchange)$r
    C[Z==0]<-rpwe(nr=n0,rate=ratec0,tchange=tchange)$r
    T[Z==1]<-rpwecx(nr=n1,rate1=rate11,rate2=rate21,rate3=rate31,
                    rate4=rate41,rate5=rate51,tchange=tchange,type=type1,rp2=rp21)$r
    T[Z==0]<-rpwecx(nr=n0,rate1=rate10,rate2=rate20,rate3=rate30,
                    rate4=rate40,rate5=rate50,tchange=tchange,type=type0,rp2=rp20)$r
    for (i in 1:nt){
      y<-pmin(pmin(T,C),t[i]-E)
      y1<-pmin(C,t[i]-E)
      delta[T<=y1]<-1
      ###log-rank
      xlr<-survdiff(Surv(y, delta) ~ Z)
      outr[r,i,1]<-(xlr$obs[2]-xlr$exp[2])/sqrt(xlr$var[2,2])
      xcox<-coxph(Surv(y, delta) ~ Z,method="breslow")
      outr[r,i,2]<-xcox$coef/sqrt(xcox$var)
      xcox0<-coxph(Surv(y, delta) ~ Z,method="breslow",iter.max =0)
      ares<-residuals(xcox0,type='score')
      amean<-mean(ares)
      outr[r,i,3]<-(xlr$obs[2]-xlr$exp[2])/sum((ares-amean)^2)
      outr[r,i,4]<-xcox$coef
    }
  }
  list(outr=outr)
}

