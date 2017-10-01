############################################################################################################################
#   A function to calculate covariance of LR and Cox
#     account for delayed treatment, discontinued treatment and non-uniform entry
#   version 1.0 (08/17/2016)
############################################################################################################################

overallcovp1<-function(tfix=2.0,tfix0=1.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                    rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                    rate41=rate21,rate51=rate51,ratec1=c(0.5,0.6),
                    rate10=rate11,rate20=rate10,rate30=rate31,
                    rate40=rate20,rate50=rate20,ratec0=ratec1,
                    tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
                    eps=1.0e-2,veps=1.0e-2,beta=0,beta0=0){
  ##tfix: the time point where the overall log hazard ratio beta is calculated
  ##tfix0: the time point where the overall log hazard ratio beta0 is calculated
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

  atchange1<-sort(unique(c(atchange,tfix,tfix-tfix0),fromLast=T))
  aind<-(atchange1[-1]-atchange1[-length(atchange1)]>err/10)
  ats<-c(0,atchange1[-1][aind==1])
  anr<-length(ats)+1
  atplus<-rep(0,anr)
  #atplus[anr]<-tmax
  atplus[anr]<-tmax+0.1*err
  atplus[-anr]<-ats
  nplus<-length(atplus)

  #anr<-length(atchange1)+1
  #atplus<-rep(0,anr)
  #atplus[anr]<-tmax
  #atplus[anr]<-tmax+0.1*err
  #atplus[-anr]<-atchange1

  covbeta1<-EA1<-0.0
  if (sum((tfix-tfix0)<=atplus & atplus<(tfix-err/10))>0){
    ats<-atplus[(tfix-tfix0)<=atplus & atplus<(tfix-err/10)]
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
    


    ST11<-pwecx(t=tk1,rate1=rate11,rate2=rate21,rate3=rate31,
          rate4=rate41,rate5=rate51,tchange=tchange,type=type1,rp2=rp21,eps=eps)$surv
    ST10<-pwecx(t=tk1,rate1=rate10,rate2=rate20,rate3=rate30,
          rate4=rate40,rate5=rate50,tchange=tchange,type=type0,rp2=rp20,eps=eps)$surv
    SC11<-pwe(t=tk1,rate=ratec1,tchange=tchange)$surv
    SC10<-pwe(t=tk1,rate=ratec0,tchange=tchange)$surv

    ST01<-pwecx(t=tk0,rate1=rate11,rate2=rate21,rate3=rate31,
          rate4=rate41,rate5=rate51,tchange=tchange,type=type1,rp2=rp21,eps=eps)$surv
    ST00<-pwecx(t=tk0,rate1=rate10,rate2=rate20,rate3=rate30,
          rate4=rate40,rate5=rate50,tchange=tchange,type=type0,rp2=rp20,eps=eps)$surv
    SC01<-pwe(t=tk0,rate=ratec1,tchange=tchange)$surv
    SC00<-pwe(t=tk0,rate=ratec0,tchange=tchange)$surv

    bb1<-(1-pi1)*ST10*SC10
    bb0<-(1-pi1)*ST00*SC00

    aa1<-pi1*exp(beta)*ST11*SC11
    aa0<-pi1*exp(beta)*ST01*SC01
    q1bs<-aa1/(aa1+bb1)
    q0bs<-aa0/(aa0+bb0)

    ak1<-innervar(t=tk1-(tfix-tfix0),taur=taur,u=u,ut=ut,pi1=pi1,
        rate11=rate11,rate21=rate21,rate31=rate31,
        rate41=rate41,rate51=rate51,ratec1=ratec1,
        rate10=rate10,rate20=rate20,rate30=rate30,
        rate40=rate40,rate50=rate50,ratec0=ratec0,
        tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
        eps=eps,veps=veps,beta=beta0)
    ak0<-innervar(t=tk0-(tfix-tfix0),taur=taur,u=u,ut=ut,pi1=pi1,
        rate11=rate11,rate21=rate21,rate31=rate31,
        rate41=rate41,rate51=rate51,ratec1=ratec1,
        rate10=rate10,rate20=rate20,rate30=rate30,
        rate40=rate40,rate50=rate50,ratec0=ratec0,
        tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
        eps=eps,veps=veps,beta=beta0)

    covbeta1<-pi1*sum(q1bs*(1-q1bs)*ak1$qf1*dk1)-(1-pi1)*sum(q0bs^2*ak0$qf1*dk0)
    covbeta1<-covbeta1-exp(beta0)*pi1*sum((1-q1bs)^2*ak1$qf2*dk1)+exp(beta0)*(1-pi1)*sum(q0bs*(1-q0bs)*ak0$qf2*dk0)
    EA1<-pi1*sum((1-q1bs)*dk1)-(1-pi1)*sum(q0bs*dk0)
  }
  list(covbeta1=covbeta1, EA1=EA1)
}

