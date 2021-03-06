\name{overallcovp2}
\alias{overallcovp2}
\title{calculate the other parts of the overall covariance}
\description{This will calculate the other parts of the overall covariance accouting for staggered entry, delayed treatment effect,  treatment crossover and loss to follow-up.}
\usage{
overallcovp2(tfix=2.0,tfix0=1.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                    rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
                    rate41=rate21,rate51=rate51,ratec1=c(0.5,0.6),
                    rate10=rate11,rate20=rate10,rate30=rate31,
                    rate40=rate20,rate50=rate20,ratec0=ratec1,
                    tchange=c(0,1),type1=1,type0=1,
                    rp21=0.5,rp20=0.5,
                    eps=1.0e-2,veps=1.0e-2,beta=0,beta0=0)
}
\arguments{
  \item{tfix}{The upper point where the overall covariance is computed.}
  \item{tfix0}{The lower point where the overall covariance is computed.}
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
  \item{rp21}{re-randomization prob in the treatment group}
  \item{rp20}{re-randomization prob in the control group}
  \item{eps}{A small number representing the error tolerance when calculating the utility function
              \deqn{\Phi_l(x)=\frac{\int_0^x s^l e^{-s}ds}{x^{l+1}}} with \eqn{l=0,1,2}.}
  \item{veps}{A small number representing the error tolerance when calculating the Fisher information.}
  \item{beta}{The value at which the covaraince is computed, upper bound}
  \item{beta0}{The value at which the covaraince is computed, lower bound}
}
\details{
The hazard functions corresponding to \code{rate11},...,\code{rate51},\code{ratec1}, \code{rate10},...,\code{rate50},\code{ratec0} are all piecewise constant function taking the form  \eqn{\lambda(t)=\sum_{j=1}^m \lambda_j I(t_{j-1}\le t<t_j)}, where \eqn{\lambda_1,\ldots,\lambda_m} are the corresponding elements of the rates and \eqn{t_0,\ldots,t_{m-1}} are the corresponding elements of \cite{tchange}, \eqn{t_m=\infty}. Note that all the rates must have the same \code{tchange}.
}
\value{\item{cov234}{The other part of the covariance}
       \item{covbeta2}{The second part of the covariance}
       \item{covbeta3}{The third part of the covariance}
       \item{covbeta4}{The fourth part of the covariance}
       \item{EA2}{The second score function}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{Version 1.0 (7/19/2016)}
\seealso{\code{\link{pwe}},\code{\link{rpwe}},\code{\link{qpwe}},\code{\link{ovbeta}},\code{\link{innervar}}}
\examples{
taur<-1.2
u<-c(1/taur,1/taur)
ut<-c(taur/2,taur)
r11<-c(1,0.5)
r21<-c(0.5,0.8)
r31<-c(0.7,0.4)
r41<-r51<-r21
rc1<-c(0.5,0.6)
r10<-c(1,0.7)
r20<-c(0.5,1)
r30<-c(0.3,0.4)
r40<-r50<-r20
rc0<-c(0.2,0.4)
getcov2<-overallcovp2(tfix=2.0,tfix0=1.0,taur=taur,u=u,ut=ut,pi1=0.5,
              rate11=r11,rate21=r21,rate31=r31,
              rate41=r41,rate51=r51,ratec1=rc1,
              rate10=r10,rate20=r20,rate30=r30,
              rate40=r40,rate50=r50,ratec0=rc0,
              tchange=c(0,1),type1=1,type0=1,
              eps=1.0e-2,veps=1.0e-2,beta=0,beta0=0)
getcov2
}
\keyword{piecewise exponential}
\keyword{piecewise uniform}
\keyword{delayed treatment effect}
\keyword{treatment crossover}

