
################################################################################################################################
###A utility function calculating $int_0^t s^k\lambda_1(s)S_2(s)ds$ where k=0,1,2 and rate1=lambda_1 and S_2 has hazard rate2
##   version 1.0 (07/19/2016)
################################################################################################################################

pwefv2<-function(t=seq(0,5,by=0.5),rate1=c(0,5,0.8),rate2=rate1,tchange=c(0,3),eps=1.0e-2){
  nt<-length(t)
  nr<-length(rate1)
  fx<-matrix(0,nrow=nt,ncol=3)
  abc2<-.Fortran("pwefv2",as.integer(nt),as.integer(nr),as.double(t),as.double(rate1),as.double(rate2),as.double(tchange),
                 as.double(eps),fx1=as.double(fx))
  fx<-matrix(abc2$fx1,byrow=FALSE,ncol=3)
  list(f0=fx[,1],f1=fx[,2],f2=fx[,3])
}








