############################################################################################################################
#   A utility function to calculate some common functions in contructing weights
#   version 1.0 (07/14/2017)
#   
#   
############################################################################################################################
wlrutil<-function(y=c(1,2,3),d=c(1,0,1),z=c(1,0,0),te=c(1,3),eps=1.0e-08){
  ##y: observed times
  ##d: non-censoring indicator
  ##z: group indicator 1=tx and 0=contl
  ##tfix: the time points where weighted log-rank is calculated
  ##te: ordered event times
  ##weights: each column is different type of weight, each row is the weight at time points te
  ##te<-sort(unique(y[d==1],fromLast=T))
  nt<-length(te)
  mn<-21
  mfunc<-matrix(0,nrow=nt,ncol=mn)
  abc2<-.Fortran("xwlrutil",as.integer(length(y)),as.double(y),as.integer(d),as.integer(z),
                 as.integer(nt),as.double(te),as.double(eps),mfunc=as.double(mfunc))
  list(mfunc=matrix(abc2$mfunc,byrow=FALSE,ncol=mn))
}
