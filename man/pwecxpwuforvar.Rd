\name{pwecxpwuforvar}
\alias{pwecxpwuforvar}
\title{calculate the utility function used for varaince calculation}
\description{This is a utility function to calculate the overall variance accouting for staggered entry, delayed treatment effect,  treatment crossover and loss to follow-up.}
\usage{
pwecxpwuforvar(tfix=10,t=seq(0,10,by=0.5),taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
    rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
         tchange=c(0,1),type=1,rp2=0.5,eps=1.0e-2)
}
\arguments{
  \item{tfix}{The upper point where the integral is computed.}
  \item{t}{A vector of lower bounds where the integral is computed.}
  \item{taur}{Recruitment time}
  \item{u}{Piecewise constant recuitment rate}
  \item{ut}{Recruitment intervals}
  \item{rate1}{Hazard before crossover}
  \item{rate2}{Hazard after crossover}
  \item{rate3}{Hazard for time to crossover}
  \item{rate4}{Hazard after crossover for complex case}
  \item{rate5}{Hazard after crossover for complex case}
  \item{ratec}{Hazard for time to censoring}
  \item{tchange}{A strictly increasing sequence of time points at which the event rates changes. The first element of tchange                 must be zero. It must have the same length as \code{rate1}, \code{rate2}, \code{rate3}, etc.}
  \item{type}{Type of crossover}
  \item{rp2}{re-randomization prob}
  \item{eps}{A small number representing the error tolerance when calculating the utility function
              \deqn{\Phi_l(x)=\frac{\int_0^x s^l e^{-s}ds}{x^{l+1}}} with \eqn{l=0,1,2}.}
}
\details{This is to calculate the function
\deqn{B_l(t,s)=\int_0^s x^l G_E(t-x)\widetilde{f}(x)S_C(x)dx,}
where \eqn{G_E} is the accrual function defined by \code{taur}, \code{u} and \code{ut}, \eqn{S_C} is the piecewise exponential survival function of the censoring time, defined by \code{tchange} and \code{ratec}, and \eqn{\widetilde{f}} is the density for the event distribution subject to crossover defined by \code{tchange}, \code{rate1} to \code{rate5} and \code{type}. This function is useful when calculating the overall varaince and covariance.
}
\value{\item{f0}{the integral when \eqn{l=0}}
       \item{f1}{the integral when \eqn{l=1}}
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
getf<-pwecxpwuforvar(tfix=3,t=seq(0,3,by=1),taur=taur,u=u,ut=ut,
                 rate1=r11,rate2=r21,rate3=r31,rate4=r41,rate5=r51,ratec=rc1,
                 tchange=c(0,1),type=1,eps=1.0e-2)
getf
}
\keyword{piecewise exponential}
\keyword{overall hazard ratio}
\keyword{piecewise uniform}
\keyword{delayed treatment effect}
\keyword{treatment crossover}
