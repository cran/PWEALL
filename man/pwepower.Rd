\name{pwepower}
\alias{pwepower}
\title{Calculating the powers of various the test statistics for superiority trials}
\description{This will calculate the powers for the test statistics accouting for staggered entry, delayed treatment effect,  treatment crossover and loss to follow-up.}
\usage{
pwepower(t=seq(0.1,3,by=0.5),alpha=0.05,twosided=1,taur=1.2,
             u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
             rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
             rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
             rate10=rate11,rate20=rate10,rate30=rate31,
             rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
             tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
             eps=1.0e-2,veps=1.0e-2,epsbeta=1.0e-4,iterbeta=25,
             n=1000)
}
\arguments{
  \item{t}{a vector of time points at which power is calculated, \code{t} must be positive}
  \item{alpha}{type-1 error rate}
  \item{twosided}{twosided test or not}
  \item{taur}{Recruitment time}
  \item{u}{Piecewise constant recuitment rate}
  \item{ut}{Recruitment intervals}
  \item{pi1}{Allocation probability for the treatment group}
  \item{rate11}{Hazard before crossover for the treatment group}
  \item{rate21}{Hazard after crossover for the treatment group}
  \item{rate31}{Hazard for time to crossover for the treatment group}
  \item{rate41}{Hazard after crossover for the treatment group for complex case}
  \item{rate51}{Hazard after crossover for the treatment group for complex case}
  \item{ratec1}{Hazard for time to censoring for the treatment group}
  \item{rate10}{Hazard before crossover for the control group}
  \item{rate20}{Hazard after crossover for the control group}
  \item{rate30}{Hazard for time to crossover for the control group}
  \item{rate40}{Hazard after crossover for the control group for complex case}
  \item{rate50}{Hazard after crossover for the control group for complex case}
  \item{ratec0}{Hazard for time to censoring for the control group}
  \item{tchange}{A strictly increasing sequence of time points at which the event rates changes. The first element of tchange                 must be zero. It must have the same length as \code{rate11}, \code{rate21}, \code{rate31}, etc.}
  \item{type1}{Type of crossover in the treatment group}
  \item{type0}{Type of crossover in the control group}
  \item{rp21}{re-randomization prob for the treatment group}
  \item{rp20}{re-randomization prob for the control group}
  \item{eps}{error tolerence}
  \item{veps}{error tolenrence for calculating variance}
  \item{epsbeta}{error tolerance for calculating overall log HR}
  \item{iterbeta}{maximum number of iterations for calculating overall log HR}
  \item{n}{total number of subjects}
}
\details{
The hazard functions corresponding to \code{rate11},...,\code{rate51},\code{ratec1}, \code{rate10},...,\code{rate50},\code{ratec0} are all piecewise constant function taking the form  \eqn{\lambda(t)=\sum_{j=1}^m \lambda_j I(t_{j-1}\le t<t_j)}, where \eqn{\lambda_1,\ldots,\lambda_m} are the corresponding elements of the rates and \eqn{t_0,\ldots,t_{m-1}} are the corresponding elements of \cite{tchange}, \eqn{t_m=\infty}. Note that all the rates must have the same \code{tchange}.
}
\value{\item{power}{powers for various test statistics. Columns 2-6 are for log-rank and columns 12-16 are for cox model. Column 2 is the exact power based on log-rank/score test; column 3 uses variance approximated by Fisher information, i.e. Lakatos's method; column 4 uses approximated Fisher info by number of events i.e. 4/D(t); column 5 uses approximated Fisher info by assuming exp dist. 1/D1(t)+1/D0(t); column 6 uses Fisher information at beta. Column 12 is the exact power based on Wald test;  column 13 uses variance approximated by Fisher information; column 14 uses approximated Fisher info by number of events i.e. 4/D(t); column 15 uses approximated Fisher info by assuming exp dist. 1/D1(t)+1/D0(t); column 16 uses Fisher information at beta=0.}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{Version 1.0 (7/19/2016)}
\seealso{\code{\link{pwe}},\code{\link{rpwe}},\code{\link{qpwe}},\code{\link{ovbeta}},\code{\link{innervar}},
         \code{\link{pwepowerni}},\code{\link{pwepowereq}}}
\examples{
t<-seq(3,6,by=1)
taur<-1.2
u<-c(1/taur,1/taur)
ut<-c(taur/2,taur)
r11<-c(0.2,0.1)
r21<-r11
r31<-c(0.03,0.02)
r41<-r51<-r21
rc1<-c(0.01,0.02)
r10<-c(0.2,0.2)
r20<-r10
r30<-c(0.02,0.01)
r40<-r50<-r20
rc0<-c(0.02,0.01)
getpower<-pwepower(t=t,alpha=0.05,twosided=1,taur=taur,u=u,ut=ut,pi1=0.5,
                   rate11=r11,rate21=r21,rate31=r31,rate41=r41,rate51=r51,ratec1=rc1,
                   rate10=r10,rate20=r20,rate30=r30,rate40=r40,rate50=r50,ratec0=rc0,
                   tchange=c(0,1),type1=1,type0=1,n=1000)
###powers at each time point
cbind(t,getpower$power[,c(2:4,12:14)])
}
\keyword{piecewise exponential}
\keyword{power}
\keyword{piecewise uniform}
\keyword{delayed treatment effect}
\keyword{treatment crossover}


