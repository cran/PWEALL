
#################################################################################
##   PieceWise Exponential Distribution
##   hazard, cumulative hazard, density, distribution and survival functions
##   version 1.0 (7/19/2016)
#################################################################################

pwe<-function(t=seq(0,5,by=0.5),rate=c(0,5,0.8),tchange=c(0,3)){
       nt<-length(t)
       nr<-length(rate)
       outr<-matrix(0,ncol=5,nrow=nt)
       abc2<-.Fortran("xpwe",as.integer(nt),as.integer(nr), as.double(t),as.double(rate),as.double(tchange),outr1=as.double(outr))
       outr<-matrix(abc2$outr1,byrow=FALSE,ncol=5)
       list(hazard=outr[,1],cumhazard=outr[,2],density=outr[,3],dist=outr[,4],surv=outr[,5])
}



