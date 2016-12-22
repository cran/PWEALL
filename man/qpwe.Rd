\name{qpwe}
\alias{qpwe}
\title{Piecewise exponential distribution: quantile function}
\description{This will provide the quantile function of the specified piecewise exponential distribution}
\usage{
qpwe(p=seq(0,1,by=0.1),rate=c(0,5,0.8),tchange=c(0,3))
}
\arguments{
  \item{p}{a vector of probabilities}
  \item{rate}{piecewise constant event rate}
  \item{tchange}{time points at which event rate changes. This must be an strictly increasing sequence starting from zero. rate and tchange must have the same length.}
}
\details{More details}
\value{\item{q}{quantiles}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This provides the quantile function related to the piecewise exponetial distribution}
\seealso{piecewise exponential}
\examples{
p<-seq(0,1,by=0.1)
rate<-c(0.6,0.3)
tchange<-c(0,1.75)
pweq<-qpwe(p=p,rate=rate,tchange=tchange)
pweq
}
\keyword{piecewise exponential}
\keyword{quantiles}% __ONLY ONE__ keyword per line
