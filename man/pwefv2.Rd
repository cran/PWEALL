\name{pwefv2}
\alias{pwefv2}
\title{A utility function}
\description{This will $int_0^t s^k lambda_1(s)S_2(s)ds$ where k=0,1,2
and rate1=lambda_1 and S_2 has hazard rate2}
\usage{
pwefv2(t=seq(0,5,by=0.5),rate1=c(0,5,0.8),
      rate2=rate1,tchange=c(0,3),eps=1.0e-2)
}
\arguments{
  \item{t}{A vector of time points}
  \item{rate1}{piecewise constant event rate}
  \item{rate2}{piecewise constant event rate}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates and tchange must have the same length.}
  \item{eps}{tolerance}
}
\details{Let \eqn{h_1,h_2} correspond to \code{rate1},\code{rate2}, and \eqn{H_1,H_2} be the corresponding survival functions. 
This function will calculate
\deqn{\int_0^t s^k h_1(s)H_2(s)ds,\hspace{1cm} k=0,1,2.}
}
\value{\item{f0}{values when \eqn{k=0}}
       \item{f1}{values when \eqn{k=1}}
       \item{f2}{values when \eqn{k=2}}
}
\references{
Luo et al. (2018) Design and monitoring of survival trials in complex scenarios, Statistics in Medicine <doi: https://doi.org/10.1002/sim.7975>.
}
\author{Xiaodong Luo}
\note{This will provide the number of events.}
\seealso{\code{\link{rpwe}}}
\examples{
r1<-c(0.6,0.3)
r2<-c(0.6,0.6)
tchange<-c(0,1.75)
pwefun<-pwefv2(t=seq(0,5,by=0.5),rate1=r1,rate2=r2,
              tchange=tchange,eps=1.0e-2)
pwefun
}
\keyword{piecewise exponential}

