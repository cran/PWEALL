#################################################################################
##   PieceWise Exponential Distribution
##   quantile function
##   version 1.0 (11/04/2016)
#################################################################################

qpwe<-function(p=seq(0,1,by=0.1),rate=c(0,5,0.8),tchange=c(0,3)){
  np<-length(p)
  nr<-length(rate)
  outr<-rep(0,np)
  abc2<-.Fortran("qpwe",as.integer(np),as.integer(nr),as.double(p),as.double(rate),as.double(tchange),outr1=as.double(outr))
  outr<-abc2$outr1
  list(q=outr)
}

