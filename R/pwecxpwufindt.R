#################################################################################
##   CALCULATE TIMELINE
##   distribution function
##   version 1.0 (11/15/2016)
#################################################################################

pwecxpwufindt<-function(target=400,ntotal=1000,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                        rate11=c(1,0.5),rate21=c(0.8,0.9),rate31=c(0.7,0.4),rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
                        rate10=c(1,0.7),rate20=c(0.9,0.7),rate30=c(0.4,0.6),rate40=rate20,rate50=rate20,ratec0=c(0.3,0.3),
                        tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,init=taur,epsilon=0.000001,maxiter=100){
  ieps<-1
  iter<-0
  tau0<-init
  ptemp<-target/ntotal
  while (ieps>epsilon & iter<=maxiter){
    iter<-iter+1
    b<-pwecxpwu(t=c(tau0),taur=taur,u=u,ut=ut,rate1=rate11,rate2=rate21,rate3=rate31,
                   rate4=rate41,rate5=rate51,ratec=ratec1,tchange=tchange,type=type1,rp2=rp21,eps=eps)
    a<-pwecxpwu(t=c(tau0),taur=taur,u=u,ut=ut,rate1=rate10,rate2=rate20,rate3=rate30,
                   rate4=rate40,rate5=rate50,ratec=ratec0,tchange=tchange,type=type0,rp2=rp20,eps=eps)
    f<-(1-pi1)*a$du+pi1*b$du-ptemp
    fprime<-(1-pi1)*a$duprime+pi1*b$duprime
    tau1<-tau0-f/fprime
    ieps<-abs(tau1-tau0)
    tau0<-tau1
  }
  tvar<-target*(ntotal-target)/ntotal^2/fprime^2
  list(tau1=tau1,tvar=tvar,eps=ieps,iter=iter)
}
