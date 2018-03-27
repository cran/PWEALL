############################################################################################################################
#   A function to simulate the power of RMST difference
#     account for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
#rmstsim<-function(tc=c(1,2),t=tc+0.2,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
rmstsim<-function(tcut=c(1,2),tstudy=tcut+0.2,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                 rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                 rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                 rate10=rate11,rate20=rate10,rate30=rate31,
                 rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                 tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
                 n=1000,rn=200,eps=1.0E-08){
  ##tc: time at which RMST is calcualted
  ##t: the time points (from FPI) where data is cut.
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
  nt<-length(tcut)
  outr<-rmst1<-rmst0<-se1<-se0<-drmst<-sed<-matrix(0,nrow=rn,ncol=nt)
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
      y<-pmin(pmin(T,C),tstudy[i]-E)
      y1<-pmin(C,tstudy[i]-E)
      delta<-rep(0,n)
      delta[T<=y1]<-1
      ###RMST
      yy1<-y[Z==1];dd1<-delta[Z==1]
      a1<-rmsth(y=yy1,d=dd1,tcut=tcut[i],eps=eps)
      yy0<-y[Z==0];dd0<-delta[Z==0]
      a0<-rmsth(y=yy0,d=dd0,tcut=tcut[i],eps=eps)
      outr[r,i]<-(a1$rmst-a0$rmst)/sqrt(a1$var/n1+a0$var/n0)
      rmst1[r,i]<-a1$rmst
      rmst0[r,i]<-a0$rmst
      se1[r,i]<-sqrt(a1$var/n1)
      se0[r,i]<-sqrt(a0$var/n0)
      drmst[r,i]<-a1$rmst-a0$rmst
      sed<-sqrt(a1$var/n1+a0$var/n0)
    }
  }
  list(outr=outr,rmst1=rmst1,se1=se1,rmst0=rmst0,se0=se0,drmst=drmst,sed=sed)
}



