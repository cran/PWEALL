\name{pwecx}
\alias{pwecx}
\title{Various function for piecewise exponential distribution with crossover effect}
\description{This will calculate the functions according to the piecewise exponential distribution with crossover}
\usage{
pwecx(t=seq(0,10,by=0.5),rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
      rate4=rate2,rate5=rate2,tchange=c(0,1),type=1,rp2=0.5,eps=1.0e-2)
}
\arguments{
  \item{t}{a vector of time points}
  \item{rate1}{piecewise constant event rate before crossover}
  \item{rate2}{piecewise constant event rate after crossover}
  \item{rate3}{piecewise constant event rate for crossover}
  \item{rate4}{additional piecewise constant event rate for more complex crossover}
  \item{rate5}{additional piecewise constant event rate for more complex crossover}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates \code{rate1} to \code{rate5} and tchange must have the same length.}
  \item{type}{type of crossover, i.e. 1: markov, 2: semi-markov, 3: hybrid case 1(as indicated in the reference), 4: hybrid case 2, 5: hybrid case 3.}
  \item{rp2}{re-randomization prob}
  \item{eps}{tolerance}
}
\details{More details}
\value{\item{hazard}{Hazard function}
       \item{cumhazard}{Cumulative hazard function}
       \item{density}{Density function}
       \item{dist}{Distribution function}
       \item{surv}{Survival function}
}
\references{
Luo et al. (2018) Design and monitoring of survival trials in complex scenarios, Statistics in Medicine <doi: https://doi.org/10.1002/sim.7975>.
}
\author{Xiaodong Luo}
\note{This provides a random number generator of the piecewise exponetial distribution with crossover}
\seealso{\code{\link{rpwe}}}
\examples{
r1<-c(0.6,0.3)
r2<-c(0.6,0.6)
r3<-c(0.1,0.2)
r4<-c(0.5,0.4)
r5<-c(0.4,0.5)
pwecxfun<-pwecx(t=seq(0,10,by=0.5),rate1=r1,rate2=r2,rate3=r3,rate4=r4,
                rate5=r5,tchange=c(0,1),type=1,eps=1.0e-2)
pwecxfun$surv
}
\keyword{piecewise exponential}
\keyword{crossover}
\keyword{random number generator}% __ONLY ONE__ keyword per line
