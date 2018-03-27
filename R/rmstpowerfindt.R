############################################################################################################################
#   A function to find the exact time when a certain power is obtained using RMST
#     account for delayed treatment, discontinued treatment and non-uniform entry
############################################################################################################################
rmstpowerfindt<-function(power=0.9,alpha=0.05,twosided=1,tcut=2,tupp=5,tlow=3.0,taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                     rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                     rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                     rate10=rate11,rate20=rate10,rate30=rate31,
                     rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
                     tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,
                     n=1000,maxiter=20,itereps=0.001){
  ##power: the desired power
  ##alpha: alpha level
  ##twodided: =1 two-sided test;=0 one-sided test
  ##tcut: time point where RMST is calculated
  ##tupp: the upper time point where the power is supposed to be bigger than the  desired
  ##tlow: the lower time point where the power is supposed to be smaller than the  desired
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
  ierr<-1.0;k<-0
  tempupp<-tupp;templow<-tlow;
  while (ierr>itereps & k<maxiter) {
    ts<-c(templow,tempupp);tempd<-tempupp-templow
    aa<-rmstpower(tcut=tcut,tstudy=ts,alpha=alpha,twosided=twosided,taur=taur,u=u,ut=ut,pi1=pi1,
                 rate11=rate11,rate21=rate21,rate31=rate31,
                 rate41=rate41,rate51=rate51,ratec1=ratec1,
                 rate10=rate10,rate20=rate20,rate30=rate30,
                 rate40=rate40,rate50=rate50,ratec0=ratec0,
                 tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
                 eps=eps,veps=veps,n=n)
    pw<-aa$power
    abspw<-abs(pw-power);ierr<-min(abspw);k<-k+1
    at<-ts[abspw<=(ierr+0.000001)]
    if (pw[2]<power) {templow<-at;tempupp<-at+tempd}
    else if (pw[1]>power){templow<-at-tempd;tempupp<-at}
    else if (at>=(tempupp-0.000001)){templow<-at-tempd/2;tempupp<-at}
    else if (at<=(templow+0.000001)){templow<-at;tempupp<-at+tempd/2}
  }
  bb<-rmstpower(tcut=tcut,tstudy=at,alpha=alpha,twosided=twosided,taur=taur,u=u,ut=ut,pi1=pi1,
                rate11=rate11,rate21=rate21,rate31=rate31,
                rate41=rate41,rate51=rate51,ratec1=ratec1,
                rate10=rate10,rate20=rate20,rate30=rate30,
                rate40=rate40,rate50=rate50,ratec0=ratec0,
                tchange=tchange,type1=type1,type0=type0,rp21=rp21,rp20=rp20,
                eps=eps,veps=veps,n=n)
  pw1<-bb$power
  list(time=at,power=pw1,err=pw1-power,iter=k)
}

