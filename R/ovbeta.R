
############################################################################################################################
#   A function to calculate the overall log hazard ratio (beta) at time point tfix
#   accounting for delayed treatment effect, treatment crossover effect, non-uniform study entry
#   and differential loss to follow-up
#   version 1.0 (07/20/2016)
#   checked 8/17/2016
############################################################################################################################
ovbeta<-function(tfix=2.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                         rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                         rate10=rate11,rate20=rate10,rate30=rate31,rate40=rate20,rate50=rate20,ratec0=c(0.4,0.3),
                         tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,beta0=0,epsbeta=1.0e-4,iterbeta=25){
  ##tfix: the time point where the overall log hazard ratio beta is calculated
  ##pi1: proportion of treatment group
  ##rate11: hazard before crossover for the treatment group
  ##rate21: hazard after crossover for the treatment group
  ##rate31: hazard for crossover for the treatment group
  ##rate41: additional hazard after crossover for the treatment group under more complex settings
  ##rate51: additional hazard after crossover for the treatment group under more complex settings
  ##ratec1: hazard for loss to follow-up for the treatment group
  ##type1: type of crossover for the treatment group

  ##rate10: hazard before crossover for the control group
  ##rate20: hazard after crossover for the control group
  ##rate30: hazard for crossover for the control group
  ##rate40: additional hazard after crossover for the control group under more complex settings
  ##rate50: additional hazard after crossover for the control group under more complex settings
  ##ratec0: hazard for loss to follow-up for the control group
  ##type0: type of crossover for the control group
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

  anr<-length(atchange)+1
  atplus<-rep(0,anr)
  atplus[anr]<-tmax+0.1*err
  atplus[-anr]<-atchange

  ats<-atplus[atplus<(tfix-0.1*err)]
  #atsupp<-c(ats,tfix)
  atsupp<-c(ats)
  nsupp<-length(atsupp)
  ###############################################################################################


  BigK1<-pwecxpwuforvar(tfix=tfix,t=atsupp,taur=taur,u=u,ut=ut,rate1=rate11,rate2=rate21,rate3=rate31,rate4=rate41,rate5=rate51,ratec=ratec1,tchange=tchange,type=type1,rp2=rp21,eps=eps)
  BigK0<-pwecxpwuforvar(tfix=tfix,t=atsupp,taur=taur,u=u,ut=ut,rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,ratec=ratec0,tchange=tchange,type=type0,rp2=rp20,eps=eps)
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

  aa1<-bb1<-(1-pi1)*ST10*SC10
  aa0<-bb0<-(1-pi1)*ST00*SC00

  b0<-b1<-beta0
  ik<-0;kerr<-1.0;bhist<-xnumhist<-xdenomhist<-rep(0,iterbeta)
  while (kerr>epsbeta & ik<=iterbeta){
    ik<-ik+1
    bhist[ik]<-b0
    aa1<-pi1*exp(b0)*ST11*SC11
    aa0<-pi1*exp(b0)*ST01*SC01
    q1bs<-aa1/(aa1+bb1)
    q0bs<-aa0/(aa0+bb0)

    xnum<-pi1*sum((1-q1bs)*dk1)-(1-pi1)*sum(q0bs*dk0)
    xdenom<-pi1*sum(q1bs*(1-q1bs)*dk1)+(1-pi1)*sum(q0bs*(1-q0bs)*dk0)
    xnumhist[ik]<-xnum;xdenomhist[ik]<-xdenom


    b1<-b0+xnum/xdenom
    kerr<-abs(b1-b0)
    b0<-b1
  }
  list(b1=b1,hr=exp(b1),err=kerr,iter=ik,bhist=bhist[1:ik],xnum=xnumhist[1:ik],xdenom=xdenomhist[1:ik],atsupp=atsupp)
}








