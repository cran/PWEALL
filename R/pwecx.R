
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
  temp<-hazard<-cumhazard<-density<-dist<-surv<-rep(0,nt)
  if (type==1){
    temp<-pwefv2(t=t,rate1=rate3,rate2=(rate1+rate3-rate2),tchange=tchange,eps=eps)$f0
    a1<-pwe(t=t,rate=rate1,tchange=tchange)
    a2<-pwe(t=t,rate=rate2,tchange=tchange)
    a13<-pwe(t=t,rate=(rate1+rate3),tchange=tchange)

    surv<-a2$surv*temp+a13$surv
    dist<-1-surv
    density<-a2$hazard*a2$surv*temp+a1$hazard*a13$surv
    cumhazard<--log(surv)
    hazard<-density/surv
  }
  else if (type==2){
    r6<-rep(0,length(rate1))+0.00000001;r4<-r5<-rate2
    tem<-tem1<-temp
    tem<-pwefvplus(t=t,rate1=rate1,rate2=rate2,rate3=rate3,rate4=r4,rate5=r5,rate6=r6,tchange=tchange,type=type,eps=eps)$f0
    temp<-pwefv2(t=t,rate1=rate3,rate2=(rate1+rate3),tchange=tchange,eps=eps)$f0
    tem1<-temp-tem
    a1<-pwe(t=t,rate=rate1,tchange=tchange)
    a2<-pwe(t=t,rate=rate2,tchange=tchange)
    a13<-pwe(t=t,rate=(rate1+rate3),tchange=tchange)

    surv<-tem1+a13$surv
    dist<-1-surv
    density<-a1$hazard*a13$surv ##this is wrong
    cumhazard<--log(surv)
    hazard<-density/surv ##this is wrong
  }
  list(hazard=hazard, cumhazard=cumhazard,density=density,dist=dist,surv=surv)
}


