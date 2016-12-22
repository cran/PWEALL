############################################################################################################################
#   A utility function to calculate inner integration in the double integration
#     that is needed for variance calculation (this is another version of function "INNERVAR")
#   version 1.0 (08/17/2016)
#   this version is currently used (8/17/2016)
############################################################################################################################

######Parameters "taur" to "ut" are not needed, we keep them as backward compliant ####8/17/2016
innervar<-function(t=seq(0,10,by=0.5),taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                       rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                       rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                       rate10=rate11,rate20=rate10,rate30=rate31,
                       rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                       tchange=c(0,1),type1=1,type0=1,
                       eps=1.0e-2,veps=1.0e-2,beta=0){
  ##t: the time points where inner part of the double integration is calculated
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

  nt<-length(t)
  ratemax<-max(c(max(c(rate21,rate20)),max(c(rate11+rate31,rate10+rate30)),max(c(ratec1,ratec0)),max(abs(c(rate11+rate31-rate21,rate10+rate30-rate20)))))
  err<-veps/ratemax
  tmax<-max(c(t,tchange,taur))+err

  nr<-length(rate11)
  tplus<-rep(0,nr)
  tplus[nr]<-tmax
  if (nr>1) tplus[-nr]<-tchange[-1]

  nn<-rep(1,nr)
  nn[1]<-ceiling((tplus[1]-tchange[1])/err)
  atchange<-rep(0,nn[1])
  atchange<-seq(tchange[1],tplus[1],by=(tplus[1]-tchange[1])/nn[1])[1:nn[1]]
  if (nr>=2){
    for (i in 2:nr){
      nn[i]<-ceiling((tplus[i]-tchange[i])/err)
      atchange<-c(atchange,seq(tchange[i],tplus[i],by=(tplus[i]-tchange[i])/nn[i])[1:nn[i]])
    }
  }

  anr<-length(atchange)+1
  atplus<-rep(0,anr)
  atplus[anr]<-tmax
  atplus[-anr]<-atchange

  atchange1<-sort(unique(c(atchange,t),fromLast=T))
  anr<-length(atchange1)+1
  atplus<-rep(0,anr)
  atplus[anr]<-tmax
  atplus[-anr]<-atchange1
  nplus<-length(atplus)

  t41<-pwefvplus(t=atplus,rate1=rate11,rate2=rate21,rate3=rate31,
                 rate4=rate41,rate5=rate51,rate6=ratec1,
                 tchange=tchange,type=type1,eps=eps)
  t40<-pwefvplus(t=atplus,rate1=rate10,rate2=rate20,rate3=rate30,
                 rate4=rate40,rate5=rate50,rate6=ratec0,
                 tchange=tchange,type=type0,eps=eps)
  t21<-pwefv2(t=atplus,rate1=rate11,rate2=rate11+rate31+ratec1,tchange=tchange,eps=eps)
  t20<-pwefv2(t=atplus,rate1=rate10,rate2=rate10+rate30+ratec0,tchange=tchange,eps=eps)

  dk1<-(t41$f0[-1]+t21$f0[-1]-t41$f0[-nplus]-t21$f0[-nplus])
  dk0<-(t40$f0[-1]+t20$f0[-1]-t40$f0[-nplus]-t20$f0[-nplus])
  tk1<-(t41$f1[-1]+t21$f1[-1]-t41$f1[-nplus]-t21$f1[-nplus])/dk1
  tk0<-(t40$f1[-1]+t20$f1[-1]-t40$f1[-nplus]-t20$f1[-nplus])/dk0

  ST11<-pwecx(t=tk1,rate1=rate11,rate2=rate21,rate3=rate31,
              rate4=rate41,rate5=rate51,
              tchange=tchange,type=type1,eps=eps)$surv
  ST10<-pwecx(t=tk1,rate1=rate10,rate2=rate20,rate3=rate30,
              rate4=rate40,rate5=rate50,
              tchange=tchange,type=type0,eps=eps)$surv
  SC11<-pwe(t=tk1,rate=ratec1,tchange=tchange)$surv
  SC10<-pwe(t=tk1,rate=ratec0,tchange=tchange)$surv

  ST01<-pwecx(t=tk0,rate1=rate11,rate2=rate21,rate3=rate31,
              rate4=rate41,rate5=rate51,
              tchange=tchange,type=type1,eps=eps)$surv
  ST00<-pwecx(t=tk0,rate1=rate10,rate2=rate20,rate3=rate30,
              rate4=rate40,rate5=rate50,
              tchange=tchange,type=type0,eps=eps)$surv
  SC01<-pwe(t=tk0,rate=ratec1,tchange=tchange)$surv
  SC00<-pwe(t=tk0,rate=ratec0,tchange=tchange)$surv

  bb1<-(1-pi1)*ST10*SC10
  bb0<-(1-pi1)*ST00*SC00
  aa1<-pi1*exp(beta)*ST11*SC11
  aa0<-pi1*exp(beta)*ST01*SC01
  r1bs<-(aa1+bb1)
  r0bs<-(aa0+bb0)
  q1bs<-aa1/r1bs
  q0bs<-aa0/r0bs

  qf1<-qf2<-rep(0,nt)
  for (i in 1:nt){
    uppi<-sum(atplus<t[i])
    if (uppi>0){
      qf1[i]<-pi1*sum(q1bs[1:uppi]/r1bs[1:uppi]*dk1[1:uppi])+(1-pi1)*sum(q0bs[1:uppi]/r0bs[1:uppi]*dk0[1:uppi])
      qf2[i]<-pi1*sum((1-q1bs[1:uppi])/r1bs[1:uppi]*dk1[1:uppi])+(1-pi1)*sum((1-q0bs[1:uppi])/r0bs[1:uppi]*dk0[1:uppi])
    }
  }
  list(qf1=qf1,qf2=qf2)
}
