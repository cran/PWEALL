\name{rmsth}
\alias{rmsth}
\title{Estimate the restricted mean survival time (RMST) and its variance from data}
\description{A function to estimate the restricted mean survival time (RMST) and its variance from data}
\usage{
rmsth(y=c(1,2,3),d=c(1,1,0),tcut=2.0,eps=1.0e-08)
}
\arguments{
  \item{y}{observed times}
  \item{d}{non-censoring indicators}
  \item{tcut}{time point at which rmst is calculated}
  \item{eps}{A small number representing the error tolerance when comparing the event times}
}
\details{More details}
\value{
  \item{tcut}{time point at which rmst is calculated}
  \item{rmst}{estimated RMST}
  \item{var}{estimated variance of \code{rmst}}
  \item{vadd}{estimated variance-covariance term of \code{rmst}}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This estimates the restricted mean survival time and its asymptotic variance}
\examples{
lamt<-0.8
lamc<-0.4
n<-3000
tcut<-2.0
truermst<-(1-exp(-lamt*tcut))/lamt
tt<-rexp(n)/lamt
cc<-rexp(n)/lamc
yy<-pmin(tt,cc)
dd<-rep(1,n)
dd[tt>cc]<-0
aest<-rmsth(y=yy,d=dd,tcut=tcut)
aest
}
\keyword{restricted mean survival time}
\keyword{variance}% __ONLY ONE__ keyword per line
