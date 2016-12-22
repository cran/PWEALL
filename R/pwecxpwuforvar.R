############################################################################################################################
#   A utility function to get the neccessary elements in calculating the variance of LR and Cox model
#   version 1.0 (7/20/2016)
#   checked 8/16/2016 no change
############################################################################################################################
pwecxpwuforvar<-function(tfix=10,t=seq(0,10,by=0.5),taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
                         rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
                         tchange=c(0,1),type=1,eps=1.0e-2){
  nt<-length(t)
  nr<-length(rate1)
  nu<-length(u)
  f0<-f1<-rep(0,nt)
  abc2<-.Fortran("pwecxpwuforvar",as.double(tfix),as.integer(nt),as.integer(nr),as.integer(nu),as.double(t),
                 as.double(taur),as.double(u),as.double(ut),as.double(rate1),as.double(rate2),
                 as.double(rate3),as.double(rate4),as.double(rate5),as.double(ratec),
                 as.double(tchange),as.integer(type),as.double(eps),f01=as.double(f0),f11=as.double(f1))
  list(f0=abc2$f01,f1=abc2$f11)
}


