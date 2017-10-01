#################################################################################
##   PieceWise Uniform Distribution
##   random number generation
##   version 1.0 (11/04/2016)
#################################################################################
rpwu<-function(nr=10,u=c(0,6,0.4),ut=c(1,2)){
  x<-1-runif(nr)
  y<-qpwu(p=x,u=u,ut=ut)$q
  list(r=y)
}
