#################################################################################
##   PieceWise Exponential Distribution
##   random number generation
##   version 1.0 (11/04/2016)
#################################################################################

rpwe<-function(nr=10,rate=c(0,5,0.8),tchange=c(0,3)){
  x<-1-runif(nr)
  y<-qpwe(p=x,rate=rate,tchange=tchange)$q
  list(r=y)
}

