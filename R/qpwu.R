#################################################################################
##   PieceWise Uniform Distribution
##   quantile function
##   version 1.0 (11/04/2016)
#################################################################################
"qpwu"<-function(p=seq(0,1,by=0.1),u=c(0,5,0.5),ut=c(1,2)){
  np<-length(p)
  nu<-length(u)
  ut0<-c(0,ut[-nu])
  rs<-sum(u*(ut-ut0))
  outr<-rep(0,np)
  if (abs(rs-1)<0.0000001) {
    abc2<-.Fortran("xqpwu",as.integer(np),as.integer(nu),as.double(p),as.double(u),as.double(ut),outr1=as.double(outr))
    outr<-abc2$outr1
  }
  list(q=outr)
}

