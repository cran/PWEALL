############################################################################################################################
#   A function to calculate powers at different timepoints for the non-inferiority trials
#   account for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
pwepowerni<-function(t=seq(0.1,3,by=0.5),nimargin=1.3,alpha=0.05,twosided=0,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,epsbeta=1.0e-4,iterbeta=25,
                     n=1000) {
  ##t: time at which power is calculated
  ##nimargin: non-inferiority margin
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
  power<-matrix(0,nrow=nt,ncol=10)
  xbeta<-rep(0,nt)
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
  alpha1<-alpha
  if (twosided==1) alpha1<-alpha/2

  bb<-qnorm(alpha1)/sqrt(avar1$xdenom*avar1$vbeta)+sqrt(n)*(log(nimargin)-a1x$b1)/sqrt(avar1$vbeta)
  bb1<-qnorm(alpha1)+sqrt(n)*(log(nimargin)-a1x$b1)/sqrt(avar1$vbeta) ###assuming Fisher information=variance of beta

  coxpower<-pnorm(bb) ### beta devided by sqrt(avar0$vbeta)
  coxpower1<-pnorm(bb1) ### beta devided by sqrt(fisher info)

  xbeta[i]<-a1x$b1
  power[i,1]<-coxpower
  power[i,2]<-coxpower1
  }
  list(beta=xbeta,power=power)
}



