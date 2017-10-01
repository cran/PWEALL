
##########################################################################################
#  The integration \int_0^t s^k f(s)G(s)ds$, $k=0,1,2$ where $f(s)$ is the density of PieceWise Exp dist 
#  with CROSSOVER effect and $G(s)$ is the survival function of the censoring time
#  version 1.0 (1/19/2017)
##########################################################################################
pwecxcens<-function(t=seq(0,10,by=0.5),rate1=c(1,0.5),rate2=rate1,
                rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,ratec=c(0.2,0.3),
                tchange=c(0,1),type=1,rp2=0.5,eps=1.0e-2){
  ##t: the time points where the function values are calculated
  ##tchange: points at which hazard changes
  ##rate1: hazard before crossover
  ##rate2: hazard after crossover
  ##rate3: hazard for crossover
  ##rate4: additional hazard after crossover
  ##rate5: additional hazard after crossover
  ##ratec: hazard for censoring
  ##type: type of crossover
  ##rp2: re-randomization prob for the semi-markov crossover (taking rate2)
  ##     rp2 is only used when type=3
  ##eps: error rate
  
  
  tin<-t
  r1<-rate1;r2<-rate2;r3<-rate3;r4<-rate4;r5<-rate5;rc<-ratec;
  ### The integration
  a2<-pwefv2(t=tin,rate1=r1,rate2=(r1+r3+rc),tchange=tchange,eps=eps)
  a4<-pwefvplus(t=tin,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,rate6=rc,tchange=tchange,type=type,rp2=rp2,eps=eps)
  
  du<-a2$f0+a4$f0
  
  ### The derivative of the integration
  b2<-pwecx(t=tin,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,tchange=tchange,type=type,rp2=rp2,eps=eps)
  b4<-pwe(t=tin,rate=rc,tchange=tchange)
  duprime<-b2$density*b4$surv
  
  list(du=du,duprime=duprime,s=b2$surv,sc=b4$surv)
}


