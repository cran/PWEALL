################################################################################################################################
###A utility function calculating $int_0^x s^le^{-s}ds/x^{l+1}$ where l=0,1,2
##   version 1.0 (07/19/2016)
################################################################################################################################

spf<-function(x=seq(-1,1,by=0.2),eps=1.0e-3){
  nx<-length(x)
  fx<-matrix(0,ncol=3,nrow=nx)
  abc2<-.Fortran("spf",as.integer(nx),as.double(x),as.double(eps),fx1=as.double(fx))
  fx<-matrix(abc2$fx1,byrow=FALSE,ncol=3)
  list(fx1=fx[,1],fx2=fx[,2],fx3=fx[,3])
}
