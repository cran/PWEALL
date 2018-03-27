############################################################################################################################
#   A function to estimate the restricted mean survival time (RMST) and its variance from data
#   version 1.0 (07/14/2017)
#   
#   
############################################################################################################################
rmsth<-function(y=c(1,2,3),d=c(1,1,0),tcut=2.0,eps=1.0e-08){
  ##y: observed times
  ##d: non-censoring indicator
  ##tcut: the time points where the restricted mean survival time (RMST) is calculated
  n<-length(y)
  te<-sort(unique(y[d==1],fromLast=T))
  nt<-length(te)
  rmst<-vrmst<-vadd<-1.0
  abc2<-.Fortran("xrmsth",as.integer(n),as.double(y),as.integer(d),as.double(tcut),as.integer(nt),as.double(te),
                         as.double(eps),rmst=as.double(rmst),vrmst=as.double(vrmst),vadd=as.double(vadd))
  list(tcut=tcut,rmst=abc2$rmst,var=abc2$vrmst,vadd=abc2$vadd)
}
