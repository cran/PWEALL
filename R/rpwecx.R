##########################################################################################
#  random number generator for PieceWise Exp dist with CROSSOVER effect
#  version 1.0 (11/04/2016)
##########################################################################################
#rPWECROSSFRx<-function(nr=1,rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,tchange=c(0,1),type=1,ISEEDS=c(13,17,19,23,29)){
rpwecx<-function(nr=1,rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,tchange=c(0,1),
                 type=1,rp2=0.5){
  # rate1: hazard before crossover
  # rate2: hazard after crossover
  # rate3: hazard for crossover
  # rate4: hazard after crossover that will be combined with rate2
  # rate5: hazard after crossover that will be combined with rate2
  # tchange: the time points at which either rate1-rate5 changes
  # type: type of crossover, i.e. markov, semi-markov,hybird
  # rp2: re-randomization probability to receive the rescue treatment when semi-markov crossover occurs.
  #      When it happens, the overall hazard will be rp2*rate2(t-s)+(1-rp2)*rate4(t), where rate2 is the hazard 
  #           for the (semi-markov) rescue treatment and rate4 is hazard for the (markov) rescue treatment

  #set.seed(ISEEDS[3])
  if (min(rate3)<1.0e-08)rate3<-rate3+1.0e-08
  t3<-rpwe(nr=nr, rate=rate3, tchange=tchange)$r
  if (type==1){
    s1t3<-pwe(t=t3,rate=rate1,tchange=tchange)$surv
    s2t3<-pwe(t=t3,rate=rate2,tchange=tchange)$surv
    #set.seed(ISEEDS[1])
    x<-runif(nr)
    x1<-x*s2t3/s1t3
    y<-rep(0,nr)
    y[x>s1t3]<-qpwe(p=1-x[x>s1t3],rate=rate1,tchange=tchange)$q
    y[x<=s1t3]<-qpwe(p=1-x1[x<=s1t3],rate=rate2,tchange=tchange)$q
  }
  else if (type==2){
    s1t3<-pwe(t=t3,rate=rate1,tchange=tchange)$surv
    #set.seed(ISEEDS[1])
    x<-runif(nr)
    x1<-x/s1t3
    y<-rep(0,nr)
    y[x>s1t3]<-qpwe(p=1-x[x>s1t3],rate=rate1,tchange=tchange)$q
    y[x<=s1t3]<-t3[x<=s1t3]+qpwe(p=1-x1[x<=s1t3],rate=rate2,tchange=tchange)$q
  }
  ###Xiaodong double-checked type=3 on 10/10/2017
  else if (type==3){
    s1t3<-pwe(t=t3,rate=rate1,tchange=tchange)$surv
    s4t3<-pwe(t=t3,rate=(1-rp2)*rate4,tchange=tchange)$surv
    #set.seed(ISEEDS[1])
    x<-runif(nr);z<-runif(nr)
    x1<-x/s1t3;z1<-z*s4t3
    y<-rep(0,nr)
    y[x>s1t3]<-qpwe(p=1-x[x>s1t3],rate=rate1,tchange=tchange)$q
    if (sum(x<=s1t3)>0){
      atemp<-t3[x<=s1t3]+qpwe(p=1-x1[x<=s1t3],rate=rp2*rate2,tchange=tchange)$q
      btemp<-qpwe(p=1-z1[x<=s1t3],rate=(1-rp2)*rate4,tchange=tchange)$q
      y[x<=s1t3]<-pmin(atemp,btemp)
    }
  }
  list(r=y)
}


