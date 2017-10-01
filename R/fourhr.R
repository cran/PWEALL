################################################################################################################################
###A utility function calculating $int_0^t h_1(s)H_2(s)h_3(t-s)H4(t-s)ds$
##   version 1.0 (01/18/2017)
################################################################################################################################

fourhr<-function(t=seq(0,5,by=0.5),rate1=c(0,5,0.8),rate2=rate1,rate3=c(0.1,0.2),rate4=rate2,
                 tchange=c(0,3),eps=1.0e-2){
  nt<-length(t)
  nr<-length(tchange)
  fx<-rep(0,nt)
  abc2<-.Fortran("xfourhr",as.integer(nt),as.integer(nr),as.double(t),as.double(rate1),as.double(rate2),as.double(rate3),
                 as.double(rate4),as.double(tchange),as.double(eps),fx=as.double(fx))
  list(fx=abc2$fx)
}

