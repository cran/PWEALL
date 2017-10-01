############################################################################################################################
#   A function to calculate powers at different timepoints
#   accounting for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
pwepower<-function(t=seq(0.1,3,by=0.5),alpha=0.05,twosided=1,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,epsbeta=1.0e-4,iterbeta=25,
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
  ##rp21, rp20 re-randomization prob for the tx and contl groups
  
  nt<-length(t)
  power<-matrix(0,nrow=nt,ncol=20)
  xevent<-matrix(0,nrow=nt,ncol=3)
  xevent[,1]<-n*(1-pi1)*pwecxpwu(t=t,taur=taur,u=u,ut=ut,rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,rate5=rate50,ratec=ratec0,
                       tchange=tchange,type=type0,rp2=rp20,eps=eps)$du
  xevent[,2]<-n*pi1*pwecxpwu(t=t,taur=taur,u=u,ut=ut,rate1=rate11,rate2=rate21,rate3=rate31,rate4=rate41,rate5=rate51,ratec=ratec1,
                       tchange=tchange,type=type1,rp2=rp21,eps=eps)$du
  xevent[,3]<-xevent[,1]+xevent[,2]
  for (i in 1:nt){
    tt<-t[i]
    a1x<-ovbeta(tfix=tt,taur=taur,u=u,ut=ut,pi1=pi1,
                rate11=rate11,rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,ratec1=ratec1,
                rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,rate50=rate50,ratec0=ratec0,
                tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,eps=eps,veps=veps,epsbeta=epsbeta,iterbeta=iterbeta)

  ###variance under the alternative
  ###variance for overall hazard ratio
  avar1<-overallvar(tfix=tt,taur=taur,u=u,ut=ut,pi1=pi1,
                    rate11=rate11,rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,ratec1=ratec1,
                    rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,rate50=rate50,ratec0=ratec0,
                    tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,eps=eps,veps=veps,beta=a1x$b1)
  ####variance for the log-rank test stat
  avar10<-overallvar(tfix=tt,taur=taur,u=u,ut=ut,pi1=pi1,
                     rate11=rate11,rate21=rate21,rate31=rate31,rate41=rate41,rate51=rate51,ratec1=ratec1,
                     rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,rate50=rate50,ratec0=ratec0,
                     tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,eps=eps,veps=veps,beta=0)

  alpha1<-alpha
  if (twosided==1) alpha1<-alpha/2
  nin<-n
  #bb<-qnorm(alpha1)*sqrt(avar0$vbeta/avar1$vbeta)-sqrt(nin)*a1x$b1/sqrt(avar1$vbeta)
  bb2<-qnorm(alpha1)/sqrt(avar1$xdenom*avar1$vbeta)-sqrt(nin)*a1x$b1/sqrt(avar1$vbeta)
  bb3<-qnorm(alpha1)-sqrt(nin)*a1x$b1*sqrt(avar1$xdenom) ###if 1/vbeta approx. I(beta,t) 
  bb4<-qnorm(alpha1)-sqrt(nin)*a1x$b1*sqrt(pi1*(1-pi1)*xevent[i,3]/nin) ###if 1/vbeta approx. I(beta,t) approx. pi1*(1-pi1)*d(t)
  bb5<-qnorm(alpha1)-sqrt(nin)*a1x$b1/sqrt(nin/xevent[i,1]+nin/xevent[i,2]) ###if 1/vbeta approx. I(beta,t) approx. 1/[1/(pi1*d1(t))+1/(pi0*d0(t))]
  ###i.e. the variance of beta is estimated by assuming exp distributions.
  bb6<-qnorm(alpha1)-sqrt(nin)*a1x$b1*sqrt(avar10$xdenom) ###if 1/vbeta approx. I(0,t) 
  
  
  aa2<-qnorm(alpha1)*sqrt(avar10$xdenom/avar10$vs)-sqrt(nin)*avar10$EA/sqrt(avar10$vs)
  aa3<-qnorm(alpha1)-sqrt(nin)*avar10$EA/sqrt(avar10$xdenom) ###if V(0,t) approx. I(0,t), this is the method by Lakatos 
  aa4<-qnorm(alpha1)-sqrt(nin)*avar10$EA/sqrt(pi1*(1-pi1)*xevent[i,3]/nin) ###if V(0,t) approx. I(0,t) approx. pi1*(1-pi1)*d(t), 
  ## this is the method using number of observed events i.e. when pi1=0.5, var=d(t)/4  
  aa5<-qnorm(alpha1)-sqrt(nin)*avar10$EA*sqrt(nin/xevent[i,1]+nin/xevent[i,2]) ###if V(0,t) approx. I(0,t) approx. 1/[1/(pi1*d1(t))+1/(pi0*d0(t))] 
  aa6<-qnorm(alpha1)-sqrt(nin)*avar10$EA/sqrt(avar1$xdenom) ###if V(0,t) approx. I(beta,t) 
  
  #coxpower<-pnorm(bb) ### beta devided by sqrt(avar0$vbeta)
  coxpower2<-pnorm(bb2) ### Exact power based on Wald test
  coxpower3<-pnorm(bb3) ### Approx. var by Fisher information
  coxpower4<-pnorm(bb4) ### Approx. Fisher info by number of events i.e. 4/D(t)
  coxpower5<-pnorm(bb5) ### Approx. Fisher info by assuming exp dist. 1/D1(t)+1/D0(t)
  coxpower6<-pnorm(bb6) ### Approx. var by Fisher information at beta=0
  
  #lrpower<-pnorm(aa)  ### un-standardized log-rank
  lrpower2<-pnorm(aa2) ### Exact power based on log-rank/score test
  lrpower3<-pnorm(aa3) ### Approx. var by Fisher information, i.e. Lakatos's method
  lrpower4<-pnorm(aa4) ### Approx. Fisher info by number of events i.e. 4/D(t)
  lrpower5<-pnorm(aa5) ### Approx. Fisher info by assuming exp dist. 1/D1(t)+1/D0(t)
  lrpower6<-pnorm(aa6) ### Approx. var by Fisher information at beta=overall beta
  
  #c(coxpower,coxpower2,lrpower,lrpower2,lrpower3)
  
  power[i,1]<-lrpower2
  power[i,2]<-lrpower2
  power[i,3]<-lrpower3
  power[i,4]<-lrpower4
  power[i,5]<-lrpower5
  power[i,6]<-lrpower6
  
  
  power[i,11]<-coxpower2
  power[i,12]<-coxpower2
  power[i,13]<-coxpower3
  power[i,14]<-coxpower4
  power[i,15]<-coxpower5
  power[i,16]<-coxpower6
  }
  list(power=power)
}



