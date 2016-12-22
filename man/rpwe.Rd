\name{rpwe}
\alias{rpwe}
\title{Piecewise exponential distribution: random number generation}
\description{This will generate random numbers according to the specified piecewise exponential distribution}
\usage{
rpwe(nr=10,rate=c(0,5,0.8),tchange=c(0,3))
}
\arguments{
  \item{nr}{number of random numbers to be generated}
  \item{rate}{piecewise constant event rate}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. rate and tchange must have the same length.}
}
\details{More details}
\value{\item{r}{random numbers}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This provides a random number generator of the piecewise exponetial distribution}
\seealso{piecewise exponential}
\examples{
nr<-10
rate<-c(0.6,0.3)
tchange<-c(0,1.75)
pwer<-rpwe(nr=nr,rate=rate,tchange=tchange)
pwer
}
\keyword{piecewise exponential}
\keyword{random number generator}% __ONLY ONE__ keyword per line
