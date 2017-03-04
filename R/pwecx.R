
##########################################################################################
#  Distribution for PieceWise Exp dist with CROSSOVER effect
#  version 1.0 (11/04/2016)
##########################################################################################
pwecx<-function(t=seq(0,10,by=0.5),rate1=c(1,0.5),rate2=rate1,
                rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,
                tchange=c(0,1),type=1,eps=1.0e-2){
  ##t: the time points where the function values are calculated
  ##tchange: points at which hazard changes
  ##rate1: hazard before crossover
  ##rate2: hazard after crossover
  ##rate3: hazard for crossover
  ##rate4: additional hazard after crossover
  ##rate5: additional hazard after crossover
  ##type: type of crossover
  
  nt<-length(t)
  r1<-rate1;r2<-rate2;r3<-rate3;r4<-rate4;r5<-rate5
  temp<-hazard<-cumhazard<-density<-dist<-surv<-rep(0,nt)
  if (type==1){
    temp<-pwefv2(t=t,rate1=r3,rate2=(r1+r3-r2),tchange=tchange,eps=eps)$f0
    a1<-pwe(t=t,rate=r1,tchange=tchange)
    a2<-pwe(t=t,rate=r2,tchange=tchange)
    a13<-pwe(t=t,rate=(r1+r3),tchange=tchange)

    surv<-a2$surv*temp+a13$surv
    dist<-1-surv
    density<-a2$hazard*a2$surv*temp+a1$hazard*a13$surv
    cumhazard<--log(surv)
    hazard<-density/surv
  }
  else if (type==2){
    txone<-rep(1,length(tchange))
    temp<-fourhr(t=t,rate1=r3,rate2=(r1+r3),rate3=r2,rate4=r2,tchange=tchange,eps=eps)$fx
    temp1<-fourhr(t=t,rate1=r3,rate2=(r1+r3),rate3=txone,rate4=r2,tchange=tchange,eps=eps)$fx
    a1<-pwe(t=t,rate=r1,tchange=tchange)
    a13<-pwe(t=t,rate=(r1+r3),tchange=tchange)

    surv<-a13$surv+temp1
    dist<-1-surv
    density<-a1$hazard*a13$surv+temp
    cumhazard<--log(surv)
    hazard<-density/surv
  }
  list(hazard=hazard, cumhazard=cumhazard,density=density,dist=dist,surv=surv)
}


