############################################################################################################################
#   A utility function calculating $int_0^t f(s)G(s)GE(t-s)ds$ and its derivative
#   where f is the density of diluted PieceWise Exp
#         GE is the dist. function of PieceWise uniform from 0 to taur
#   version 1.0 (7/20/2016)
############################################################################################################################
pwecxpwu<-function(t=seq(0,10,by=0.5),taur=5,
                  u=c(1/taur,1/taur),ut=c(taur/2,taur),
                  rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
                  rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
                  tchange=c(0,1),type=1,eps=1.0e-2){
  nt<-length(t)
  nr<-length(rate1)
  nu<-length(u)
  du<-duprime<-rep(0,nt)
  abc2<-.Fortran("pwecxpwu",as.integer(nt),as.integer(nr),as.integer(nu),
                 as.double(t),as.double(taur),as.double(u),as.double(ut),
                 as.double(rate1),as.double(rate2),as.double(rate3),
                 as.double(rate4),as.double(rate5),as.double(ratec),
                 as.double(tchange),as.integer(type),as.double(eps),
                 du1=as.double(du),duprime1=as.double(duprime))
  list(du=abc2$du1,duprime=abc2$duprime1)
}
