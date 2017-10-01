#################################################################################
##   PieceWise Uniform Distribution
##   distribution function
##   version 1.0 (11/04/2016)
#################################################################################
"pwu"<-function(t=seq(0,1,by=0.1),u=c(0,5,0.5),ut=c(1,2)){
  nt<-length(t)
  nu<-length(u)
  ut0<-c(0,ut[-nu])
  rs<-sum(u*(ut-ut0))
  outr<-rep(0,nt)
  if (abs(rs-1)<0.0000001) {
    abc2<-.Fortran("xpwu",as.integer(nt),as.integer(nu),as.double(t),as.double(u),as.double(ut),outr1=as.double(outr))
    outr<-abc2$outr1
  }
  list(dist=outr)
}
