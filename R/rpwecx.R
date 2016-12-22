##########################################################################################
#  random number generator for PieceWise Exp dist with CROSSOVER effect
#  version 1.0 (11/04/2016)
##########################################################################################
#rPWECROSSFRx<-function(nr=1,rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,tchange=c(0,1),type=1,ISEEDS=c(13,17,19,23,29)){
rpwecx<-function(nr=1,rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,tchange=c(0,1),type=1){
  # rate1: hazard before crossover
  # rate2: hazard after crossover
  # rate3: hazard for crossover
  # rate4: hazard after crossover that will be combined with rate2
  # rate5: hazard after crossover that will be combined with rate2
  # tchange: the time points at which either rate1-rate5 changes
  # type: type of crossover, i.e. markov, semi-markov,

  #set.seed(ISEEDS[3])
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
  list(r=y)
}


