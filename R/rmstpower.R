############################################################################################################################
#   A function to calculate powers at different cut-points based on difference of restricted mean survival times (RMST)
#   account for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
#rmstpower<-function(tfix=2,tcut=seq(tfix,tfix+2,by=0.5),alpha=0.05,twosided=1,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
rmstpower<-function(tcut=2,tstudy=seq(tcut,tcut+2,by=0.5),alpha=0.05,twosided=1,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,
                     n=1000) {
  ##tcut: time at which RMST is calculated
  ##tstudy: time points (from FPI) where powers are calculated 
  ##alpha: alpha level
  ##twodided: =1 two-sided test;=0 one-sided test
  ##taur: recruitment time
  ##u: recruitment rate in each interval
  ##ut: recruitment intervals
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
  ##type1: type of crossover for the treatment group
  ##type0: type of crossover for the control group
  ##n: total number of subject to be recruited
  ##rp21, rp20 re-randomization prob for the tx and contl groups
  

  nt<-length(tstudy)
  power<-rmst1<-rmst0<-drmst<-vrmst1<-vrmst0<-vdrmst<-rep(0,nt)
  n1<-n*pi1;n0<-n-n1
  alpha1<-alpha
  if (twosided==1) alpha1<-alpha/2
  for (i in 1:nt){
    t1<-tstudy[i]
    r1<-rate11;r2<-rate21;r3<-rate31;r4<-rate41;r5<-rate51;rc<-ratec1;
    bb<-rmstutil(tcut=tcut,tstudy=t1,taur=taur,u=u,ut=ut,
                 rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,ratec=rc,
                 tchange=tchange,type=type1,rp2=rp21,eps=eps,veps=veps)
    
    r1<-rate10;r2<-rate20;r3<-rate30;r4<-rate40;r5<-rate50;rc<-ratec0;
    aa<-rmstutil(tcut=tcut,tstudy=t1,taur=taur,u=u,ut=ut,
                 rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,ratec=rc,
                 tchange=tchange,type=type0,rp2=rp20,eps=eps,veps=veps)

    rmst1[i]<-bb$rmst;vrmst1[i]<-bb$var/n1
    rmst0[i]<-aa$rmst;vrmst0[i]<-aa$var/n0
    drmst[i]<-rmst1[i]-rmst0[i]
    vdrmst[i]<-vrmst1[i]+vrmst0[i]
    
    
    pb<-qnorm(1-alpha1)-drmst[i]/sqrt(vdrmst[i])
    power[i]<-1-pnorm(pb)
  }  

  list(power=power,rmst1=rmst1,se1=sqrt(vrmst1),rmst0=rmst0,se0=sqrt(vrmst0),drmst=drmst,sed=sqrt(vdrmst))
}
  



