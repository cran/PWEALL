############################################################################################################################
#   A function to calculate variance of LR and Cox
#     account for delayed treatment, discontinued treatment and non-uniform entry
#   version 1.0 (07/20/2016)
#   version 2.0 (08/17/2016)
#   version 3.0 (10/26/2016) corrected the cut-points definition, see Part A
############################################################################################################################
overallvar<-function(tfix=2.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
                     eps=1.0e-2,veps=1.0e-2,beta=0){
  ##tfix: the time point where the overall log hazard ratio beta is calculated
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
  
  ######Part A################################################################################
  #ratemax<-max(c(max(c(rate21,rate20)),max(c(rate11+rate31,rate10+rate30)),max(c(ratec1,ratec0)),max(abs(c(rate11+rate31-rate21,rate10+rate30-rate20)))))
  #Changed 3/2/2017 so that ratemax won't be too big or too small
  ratemax<-max(abs(rate11-rate10))+max(abs(rate21-rate20))+max(abs(rate31-rate30))+max(abs(rate41-rate40))+max(abs(rate51-rate50))+max(abs(ratec1-ratec0)) #Changed 3/1/2017 so that ratemax won't be too big
  rateb<-max(0.01,min(ratemax,1))
  err<-veps/rateb
  tmax<-max(c(tfix,tchange,taur))+err

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

  ###############################################################################################

  #On 8/17/2016 we use this to replace the following
  atchange1<-sort(unique(c(atchange,tfix),fromLast=T))
  anr<-length(atchange1)+1
  atplus<-rep(0,anr)
  atplus[anr]<-tmax
  atplus[-anr]<-atchange1

  ##Replaced############
  #anr<-length(atchange)+1
  #atplus<-rep(0,anr)
  #atplus[anr]<-tmax
  #atplus[-anr]<-atchange

  ats<-atplus[atplus<(tfix-0.1*err)]
  atsupp<-c(ats,tfix)
  nsupp<-length(atsupp)

  BigK1<-pwecxpwuforvar(tfix=tfix,t=atsupp,taur=taur,u=u,ut=ut,
            rate1=rate11,rate2=rate21,rate3=rate31,
            rate4=rate41,rate5=rate51,ratec=ratec1,
            tchange=tchange,type=type1,rp2=rp21,eps=eps)
  BigK0<-pwecxpwuforvar(tfix=tfix,t=atsupp,taur=taur,u=u,ut=ut,
            rate1=rate10,rate2=rate20,rate3=rate30,
            rate4=rate40,rate5=rate50,ratec=ratec0,
            tchange=tchange,type=type0,rp2=rp20,eps=eps)
  dk1<-BigK1$f0[-1]-BigK1$f0[-nsupp]
  dk0<-BigK0$f0[-1]-BigK0$f0[-nsupp]
  ###12/22/2016
  #tk1<-(BigK1$f1[-1]-BigK1$f1[-nsupp])/dk1
  #tk0<-(BigK0$f1[-1]-BigK0$f1[-nsupp])/dk0
  ###3/1/2017
  adk1<-(dk1>1.0e-08)
  adk0<-(dk0>1.0e-08)
  tk1<-tk0<-atsupp[-nsupp]
  tk1[adk1==1]<-(BigK1$f1[-1]-BigK1$f1[-nsupp])[adk1==1]/dk1[adk1==1]
  tk0[adk0==1]<-(BigK0$f1[-1]-BigK0$f1[-nsupp])[adk0==1]/dk0[adk0==1]
  
  ST11<-pwecx(t=tk1,rate1=rate11,rate2=rate21,rate3=rate31,rate4=rate41,rate5=rate51,tchange=tchange,type=type1,rp2=rp21,eps=eps)$surv
  ST10<-pwecx(t=tk1,rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,tchange=tchange,type=type0,rp2=rp20,eps=eps)$surv
  SC11<-pwe(t=tk1,rate=ratec1,tchange=tchange)$surv
  SC10<-pwe(t=tk1,rate=ratec0,tchange=tchange)$surv

  ST01<-pwecx(t=tk0,rate1=rate11,rate2=rate21,rate3=rate31,rate4=rate41,rate5=rate51,tchange=tchange,type=type1,rp2=rp21,eps=eps)$surv
  ST00<-pwecx(t=tk0,rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,tchange=tchange,type=type0,rp2=rp20,eps=eps)$surv
  SC01<-pwe(t=tk0,rate=ratec1,tchange=tchange)$surv
  SC00<-pwe(t=tk0,rate=ratec0,tchange=tchange)$surv

  bb1<-(1-pi1)*ST10*SC10
  bb0<-(1-pi1)*ST00*SC00

  aa1<-pi1*exp(beta)*ST11*SC11
  aa0<-pi1*exp(beta)*ST01*SC01
  q1bs<-aa1/(aa1+bb1)
  q0bs<-aa0/(aa0+bb0)

  EA<-pi1*sum((1-q1bs)*dk1)-(1-pi1)*sum(q0bs*dk0)
  EA2<-pi1*sum((1-q1bs)^2*dk1)+(1-pi1)*sum(q0bs^2*dk0)
  xdenom<-pi1*sum(q1bs*(1-q1bs)*dk1)+(1-pi1)*sum(q0bs*(1-q0bs)*dk0)

  #ak1<-INNERVAR(t=tk1,taur=taur,u=u,ut=ut,pi1=pi1,rate11=rate11,rate21=rate21,rate31=rate31,ratec1=ratec1,rate10=rate10,rate20=rate20,rate30=rate30,ratec0=ratec0,tchange=tchange,
  #              eps=eps,veps=veps,beta=beta)
  #ak0<-INNERVAR(t=tk0,taur=taur,u=u,ut=ut,pi1=pi1,rate11=rate11,rate21=rate21,rate31=rate31,ratec1=ratec1,rate10=rate10,rate20=rate20,rate30=rate30,ratec0=ratec0,tchange=tchange,
  #              eps=eps,veps=veps,beta=beta)
  ak1<-innervar(t=tk1,taur=taur,u=u,ut=ut,pi1=pi1,rate11=rate11,
                rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,
                ratec1=ratec1,rate10=rate10,rate20=rate20,rate30=rate30,
                rate40=rate40,rate50=rate50,ratec0=ratec0,
                tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
                eps=eps,veps=veps,beta=beta)
  ak0<-innervar(t=tk0,taur=taur,u=u,ut=ut,pi1=pi1,rate11=rate11,
                rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,
                ratec1=ratec1,rate10=rate10,rate20=rate20,rate30=rate30,
                rate40=rate40,rate50=rate50,ratec0=ratec0,rp21=rp21,rp20=rp20,
                tchange=tchange,type1=type1,type0=type0,
                eps=eps,veps=veps,beta=beta)

  AB<-pi1*sum(q1bs*(1-q1bs)*ak1$qf1*dk1)-(1-pi1)*sum(q0bs^2*ak0$qf1*dk0)
  AB<-AB-exp(beta)*pi1*sum((1-q1bs)^2*ak1$qf2*dk1)+exp(beta)*(1-pi1)*sum(q0bs*(1-q0bs)*ak0$qf2*dk0)
  xnum<-EA2-EA^2+2*AB
  vbeta<-xnum/xdenom^2
  vs<-xnum
  list(vbeta=vbeta,vs=vs,xdenom=xdenom,EA=EA, EA2=EA2,AB=AB)
}
