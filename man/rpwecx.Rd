\name{rpwecx}
\alias{rpwecx}
\title{Piecewise exponential distribution with crossover effect: random number generation}
\description{This will generate random numbers according to the piecewise exponential distribution with crossover}
\usage{
rpwecx(nr=1,rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
rate4=rate2,rate5=rate2,tchange=c(0,1),type=1,rp2=0.5)
}
\arguments{
  \item{nr}{number of random numbers to be generated}
  \item{rate1}{piecewise constant event rate before crossover}
  \item{rate2}{piecewise constant event rate after crossover}
  \item{rate3}{piecewise constant event rate for crossover}
  \item{rate4}{additional piecewise constant event rate for more complex crossover}
  \item{rate5}{additional piecewise constant event rate for more complex crossover}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates \code{rate1} to \code{rate6} and tchange must have the same length.}
  \item{type}{type of crossover, 1=markov, 2=semi-markov, 3=hybrid}
  \item{rp2}{re-randomization probability to receive the rescue treatment when semi-markov crossover occurs.
        When it happens, the overall hazard will be pi2*r2(t-s)+(1-pi2)*r4(t), where r2 is the hazard 
             for the semi-markov rescue treatment and r4 is hazard for the markov rescue treatment.} 
}
\details{More details}
\value{\item{r}{random numbers for the event time}
       \item{rx}{random numbers for the crossover time}
       \item{cxind}{indicators for the crossover, the first column indicates whether crossover occurs, i.e. \eqn{\code{rx}<\code{r}}. When \code{type}=3,4,5, the second column of \code{cxind} indicates whether it crosses to the arm with \code{rate2}}
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
pwecxr<-rpwecx(nr=10,rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,tchange=c(0,1),type=1)
pwecxr$r
}
\keyword{piecewise exponential}
\keyword{treatment crossover}
\keyword{random number generator}% __ONLY ONE__ keyword per line
