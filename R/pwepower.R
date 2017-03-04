############################################################################################################################
#   A function to calculate powers at different timepoints
#   account for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
pwepower<-function(t=seq(0.1,3,by=0.5),alpha=0.05,twosided=1,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,eps=1.0e-2,veps=1.0e-2,epsbeta=1.0e-4,iterbeta=25,
                     n=1000) {
  ##t: time at which power is calculated
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
  ##testtype: type of statistic the power calculation is based =1, log-rank; =2, Cox model; =3, log-rank with robust var; =4, overall (log)HR

  nt<-length(t)
  power<-matrix(0,nrow=nt,ncol=10)
  for (i in 1:nt){
    tt<-t[i]
    a1x<-ovbeta(tfix=tt,taur=taur,u=u,ut=ut,pi1=pi1,
                rate11=rate11,rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,ratec1=ratec1,
                rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,rate50=rate50,ratec0=ratec0,
                tchange=tchange,type1=type1,type0=type0,eps=eps,veps=veps,epsbeta=epsbeta,iterbeta=iterbeta)

  ###variance under the alternative
  ###variance for overall hazard ratio
  avar1<-overallvar(tfix=tt,taur=taur,u=u,ut=ut,pi1=pi1,
                    rate11=rate11,rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,ratec1=ratec1,
                    rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,rate50=rate50,ratec0=ratec0,
                    tchange=tchange,type1=type1,type0=type0,eps=eps,veps=veps,beta=a1x$b1)
  ####variance for the log-rank test stat
  avar10<-overallvar(tfix=tt,taur=taur,u=u,ut=ut,pi1=pi1,
                     rate11=rate11,rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,ratec1=ratec1,
                     rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,rate50=rate50,ratec0=ratec0,
                     tchange=tchange,type1=type1,type0=type0,eps=eps,veps=veps,beta=0)

  alpha1<-alpha
  if (twosided==1) alpha1<-alpha/2
  nin<-n
  #bb<-qnorm(alpha1)*sqrt(avar0$vbeta/avar1$vbeta)-sqrt(nin)*a1x$b1/sqrt(avar1$vbeta)
  bb1<-qnorm(alpha1)*sqrt(avar1$xdenom*avar1$vbeta)-sqrt(nin)*a1x$b1/sqrt(avar1$vbeta)
  #aa<-qnorm(alpha1)*sqrt(avar0$vs/avar10$vs)-sqrt(nin)*avar10$EA/sqrt(avar10$vs)
  aa1<-qnorm(alpha1)*sqrt(avar10$xdenom/avar10$vs)-sqrt(nin)*avar10$EA/sqrt(avar10$vs)
  aa2<-qnorm(alpha1)-sqrt(nin)*avar10$EA/sqrt(avar10$xdenom)


  #coxpower<-pnorm(bb) ### beta devided by sqrt(avar0$vbeta)
  coxpower2<-pnorm(bb1) ### beta devided by sqrt(fisher info)
  #lrpower<-pnorm(aa)  ### un-standardized log-rank
  lrpower2<-pnorm(aa1) ### log-rank
  lrpower3<-pnorm(aa2) ### log-rank based on approximation
  #c(coxpower,coxpower2,lrpower,lrpower2,lrpower3)

  power[i,1]<-lrpower2
  power[i,2]<-coxpower2
  power[i,3]<-lrpower3
  }

  list(power=power)
}



