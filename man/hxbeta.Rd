\name{hxbeta}
\alias{hxbeta}
\title{A function to calculate the beta-smoothed hazard rate}
\description{A function to calculate the beta-smoothed hazard rate}
\usage{
hxbeta(x=c(0.5,1),y=seq(.1,1,by=0.01),d=rep(1,length(y)),
           tfix=2,K=20,eps=1.0e-06)
}
\arguments{
  \item{x}{time points where the estimated hazards are calculated}
  \item{y}{observed times}
  \item{d}{non-censoring indicators}
  \item{tfix}{maximum time point at which the hazard function is estimated}
  \item{K}{smooth parameter for the hazard estimate}
  \item{eps}{the error tolerance when comparing event times}
}
\details{V1:3/21/2018}
\value{
  \item{lambda}{estimated hazard at points \code{x}}
}
\author{Xiaodong Luo}
\examples{
n<-200
taur<-2.8
u<-c(1/taur,1/taur)
ut<-c(taur/2,taur)
tfix<-taur+2
tseq<-seq(0,tfix,by=0.1)
r11<-c(1,0.5)
r21<-c(0.5,0.8)
r31<-c(0.7,0.4)
r41<-r51<-r21
rc1<-c(0.5,0.6)
tchange<-c(0,1.873)

E<-T<-C<-d<-rep(0,n)
E<-rpwu(nr=n,u=u,ut=ut)$r
C<-rpwe(nr=n,rate=rc1,tchange=tchange)$r
T<-rpwecx(nr=n,rate1=r11,rate2=r21,rate3=r31,
               rate4=r41,rate5=r51,tchange=tchange,type=1)$r
y<-pmin(pmin(T,C),tfix-E)
y1<-pmin(C,tfix-E)
d[T<=y]<-1

lambda=hxbeta(x=tseq,y=y,d=d,tfix=tfix,K=20,eps=1.0e-06)$lambda
lambda
}
\keyword{smoothed estimate}
\keyword{hazard estimate}
