############################################################################################################################
#   A utility function to calculate weighted log-rank test and its variance given the weights
#   version 1.0 (07/14/2017)
#   
#   
############################################################################################################################
wlrcal<-function(n=10,te=c(1,2,3),tfix=2.0,dd1=c(1,0,1),dd0=c(0,1,0),r1=c(1,2,3),r0=c(1,2,3),weights=matrix(1,nrow=length(te),ncol=1),eps=1.0e-08){
  ##y: observed times
  ##d: non-censoring indicator
  ##z: group indicator 1=tx and 0=contl
  ##tfix: the time points where weighted log-rank is calculated
  ##te: ordered event times
  ##weights: each column is different type of weight, each row is the weight at time points te
  weights<-as.matrix(weights)
  nt<-length(te)
  nw<-ncol(weights)
  xtest<-xvtest<-xlr<-rep(0,nw)
  xlcor<-matrix(0,nrow=nw,ncol=nw)
  abc2<-.Fortran("xwlrcal",as.integer(n),as.integer(nt),as.double(te),as.double(tfix),
                 as.double(dd1),as.double(dd0),as.double(r1),as.double(r0),
                 as.integer(nw),as.double(weights),as.double(eps),
                 xtest=as.double(xtest),xvtest=as.double(xvtest),xlr=as.double(xlr),xlcor=as.double(xlcor))
  list(test=abc2$xtest,var=abc2$xvtest,wlr=abc2$xlr,wlcor=matrix(abc2$xlcor,ncol=nw,byrow=TRUE))
}


