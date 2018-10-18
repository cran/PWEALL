\name{fourhr}
\alias{fourhr}
\title{A utility functon}
\description{This will calculate the more complex integration}
\usage{
fourhr(t=seq(0,5,by=0.5),rate1=c(0,5,0.8),rate2=rate1,
                   rate3=c(0.1,0.2),rate4=rate2,tchange=c(0,3),eps=1.0e-2)
}
\arguments{
  \item{t}{A vector of time points}
  \item{rate1}{piecewise constant event rate}
  \item{rate2}{piecewise constant event rate}
  \item{rate3}{piecewise constant event rate}
  \item{rate4}{additional piecewise constant}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates and tchange must have the same length.}
  \item{eps}{tolerance}
}
\details{Let \eqn{h_1,\ldots,h_4} correspond to \code{rate1},...,\code{rate4}, and \eqn{H_1,\ldots,H_4} be the corresponding survival functions. 
We calculate
\deqn{\int_0^t h_1(s)H_2(s)h_3(t-s)H_4(t-s)ds.}
}
\value{\item{fx}{values}
}
\references{
Luo et al. (2018) Design and monitoring of survival trials in complex scenarios, Statistics in Medicine <doi: https://doi.org/10.1002/sim.7975>.
}
\author{Xiaodong Luo}
\note{This provides the result of the complex integration}
\seealso{\code{\link{rpwe}}}
\examples{
r1<-c(0.6,0.3)
r2<-c(0.6,0.6)
r3<-c(0.1,0.2)
r4<-c(0.5,0.4)
tchange<-c(0,1.75)
fourhrfun<-fourhr(t=seq(0,5,by=0.5),rate1=r1,rate2=r2,rate3=r3,
                 rate4=r4,tchange=c(0,3),eps=1.0e-2)
fourhrfun
}
\keyword{piecewise exponential}
\keyword{treatment crossover}

