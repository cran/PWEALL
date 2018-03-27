############################################################################################################################
#   A utility function to calculate the true restricted mean survival time (RMST) and its variance
#     account for delayed treatment, discontinued treatment and non-uniform entry
#   version 1.0 (07/13/2017)
#   version 2.0 (12/20/2017) #"tfix" is changed to "tcut", "tcut" is changed to "tstudy"   
#   
#   
############################################################################################################################
#rmstutil<-function(tfix=2.0,tcut=5.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
rmstutil<-function(tcut=2.0,tstudy=5.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
                     rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
                     rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
                     tchange=c(0,1),type=1,rp2=0.5,
                     eps=1.0e-2,veps=1.0e-2){
  ##tcut: the cut-point where the restricted mean survival times (RMST) is calculated
  ##tstudy: the study time point from the first patient in, it must be larger than tcut. Combining with the recruitment curve, 
  ##      this impacts the censoring distribution and the variance of the RMST
  ##taur: recruitment period
  ##u,ut: recruitment curves
  ##rate1: hazard before crossover
  ##rate2: hazard after crossover
  ##rate3: hazard for crossover
  ##rate4: additonal hazard after crossover
  ##rate5: additonal hazard for crossover
  ##ratec: hazard for loss to follow-up
  ##tchange: points at which hazard rate changes
  ##type: type of crossover 1(default): Markov, 2: Semi-Markov, 3: Hybrid
  ##eps,veps: error tolerance for numerical integration
  ##rp2: re-randomization prob
  
  tin<-c(tcut)
  r1<-rate1;r2<-rate2;r3<-rate3;r4<-rate4;r5<-rate5;rc<-rep(0,length(tchange));
  ### The integration
  a2<-pwefv2(t=tin,rate1=r1,rate2=(r1+r3+rc),tchange=tchange,eps=eps)
  a4<-pwefvplus(t=tin,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,rate6=rc,tchange=tchange,type=type,rp2=rp2,eps=eps)
  Lrmst<-tcut-tcut*(a2$f0[1]+a4$f0[1])+(a2$f1[1]+a4$f1[1])

  
  ######Part A################################################################################
  ratemax<-max(rate1)+max(rate2)+max(rate3)+max(rate4)+max(rate5)+max(ratec)
  rateb<-max(0.01,min(ratemax,1))
  err<-veps/rateb
  tmax<-max(c(tstudy,tcut,tchange,taur))+err

  nr<-length(rate1)
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
  atchange1<-sort(unique(c(atchange,tstudy),fromLast=T))
  anr<-length(atchange1)+1
  atplus<-rep(0,anr)
  atplus[anr]<-tmax
  atplus[-anr]<-atchange1

  ats<-atplus[atplus<(tstudy-0.1*err)]
  atsupp<-c(ats,tstudy)
  nsupp<-length(atsupp)

  ### The integration
  ka2<-pwefv2(t=atsupp,rate1=r1,rate2=(r1+r3+rc),tchange=tchange,eps=eps)
  ka4<-pwefvplus(t=atsupp,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,rate6=rc,tchange=tchange,type=type,rp2=rp2,eps=eps)
  
  dk<-(ka2$f0[-1]+ka4$f0[-1])-(ka2$f0[-nsupp]+ka4$f0[-nsupp])
  
  adk<-(dk>1.0e-08)
  tk<-atsupp[-nsupp]
  tk[adk==1]<-(ka2$f1[-1]+ka4$f1[-1]-ka2$f1[-nsupp]-ka4$f1[-nsupp])[adk==1]/dk[adk==1]

  ST<-pwecx(t=tk,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,tchange=tchange,type=type,rp2=rp2,eps=eps)$surv
  SC<-pwe(t=tk,rate=ratec,tchange=tchange)$surv
  HT<-pwu(t=tstudy-tk,u=u,ut=ut)$dist
  HTcut<-pwu(t=tstudy,u=u,ut=ut)$dist
  HT[tstudy<=tk]<-(-1.0)

  ka2<-pwefv2(t=tk,rate1=r1,rate2=(r1+r3+rc),tchange=tchange,eps=eps)
  ka4<-pwefvplus(t=tk,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,rate6=rc,tchange=tchange,type=type,rp2=rp2,eps=eps)
  
  Lr<-tk-tk*(ka2$f0+ka4$f0)+(ka2$f1+ka4$f1)

  temp1<-(Lrmst-Lr)/(ST^2*SC*HT)*dk
  temp2<-(Lrmst-Lr)^2/(ST^2*SC*HT)*dk
  
  vfix<-HTcut*sum(temp2[tk<=tcut])
  vadd<-HTcut*sum(temp1[tk<=tcut])
  
  list(tcut=tcut,tstudy=tstudy,rmst=Lrmst,var=vfix,vadd=vadd)
}
