############################################################################################################################
#   A function to calculate the variance and covariance of estimated restricted mean survival time 
#     using data from different cut-off points
#     accounting for delayed treatment, discontinued treatment and non-uniform entry
#   version 1.0 (07/13/2017)
#   version 2.0 (12/20/2017) #"tfix" is changed to "tcut", "tcut" is changed to "tstudy"   
#   
############################################################################################################################
#rmstcov<-function(tfix=2.0,tcut=2.5,tfix1=3.0,tcut1=3.5,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
rmstcov<-function(t1cut=2.0,t1study=2.5,t2cut=3.0,t2study=3.5,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
                     rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
                     rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
                     tchange=c(0,1),type=1,rp2=0.5,
                     eps=1.0e-2,veps=1.0e-2){
  ##t1cut<=t2cut: the cut-points where the restricted mean survival times (RMST) are calculated
  ##t1study<=t2study: the study time points from first patient in, it must be larger than tcut. This will be used for study monitoring
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
  
  r1<-rate1;r2<-rate2;r3<-rate3;r4<-rate4;r5<-rate5;rc<-ratec;
  #aa<-rmstutil(tfix=tfix,tcut=tcut,taur=taur,u=u,ut=ut,
  aa<-rmstutil(tcut=t1cut,tstudy=t1study,taur=taur,u=u,ut=ut,
                     rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,ratec=rc,
                     tchange=tchange,type=type,rp2=rp2,eps=eps,veps=veps)
  bb<-rmstutil(tcut=t2cut,tstudy=t2study,taur=taur,u=u,ut=ut,
               rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,ratec=rc,
               tchange=tchange,type=type,rp2=rp2,eps=eps,veps=veps)
  cc<-rmstutil(tcut=t1cut,tstudy=t2study,taur=taur,u=u,ut=ut,
               rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,ratec=rc,
               tchange=tchange,type=type,rp2=rp2,eps=eps,veps=veps)
  xcov<-cc$var+(cc$vadd)*(bb$rmst-aa$rmst)
  xcov1<-cc$var+(cc$vadd)*(bb$rmst-cc$rmst)
  #list(tfix=tfix,tcut=tcut,tfix1=tfix1,tcut1=tcut1,rmst=aa$rmst,rmst1=bb$rmst,rmstx=cc$rmst,v=aa$var,v1=bb$var,cov=xcov,cov1=xcov1)
  list(t1cut=t1cut,t1study=t1study,t2cut=t2cut,t2study=t2study,rmst=aa$rmst,rmst1=bb$rmst,rmstx=cc$rmst,v=aa$var,v1=bb$var,cov=xcov,cov1=xcov1)
}
