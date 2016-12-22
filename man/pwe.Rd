\name{pwe}
\alias{pwe}
\title{Piecewise exponential distribution: hazard, cumulative hazard, density, distribution, survival}
\description{This will provide the related functions of the specified piecewise exponential distribution.}
\usage{
pwe(t=seq(0,5,by=0.5),rate=c(0,5,0.8),tchange=c(0,3))
}
\arguments{
  \item{t}{A vector of time points.}
  \item{rate}{A vector of event rates}
  \item{tchange}{A strictly increasing sequence of time points at which the event rate changes. The first element of tchange must be zero. It must have the same length as \code{rate}.}
}
\details{Let \eqn{\lambda(t)=\sum_{j=1}^m \lambda_j I(t_{j-1}\le t<t_j)} be the hazard function, where \eqn{\lambda_1,\ldots,\lambda_m} are the corresponding elements of \cite{rate} and \eqn{t_0,\ldots,t_{m-1}} are the corresponding elements of \cite{tchange}, \eqn{t_m=\infty}. The cumulative hazard function
\deqn{\Lambda(t)=\sum_{j=1}^m \lambda_j(t\wedge t_j-t\wedge t_{j-1}),}
the survival function \eqn{S(t)=\exp\{-\Lambda(t)\}}, the distribution function \eqn{F(t)=1-S(t)} and the density function
\eqn{f(t)=\lambda(t)S(t)}.
}
\value{\item{hazard}{Hazard function}
       \item{cumhazard}{Cumulative hazard function}
       \item{density}{Density function}
       \item{dist}{Distribution function}
       \item{surv}{Survival function}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{Version 1.0 (7/19/2016)}
\seealso{\code{\link{rpwe}},\code{\link{qpwe}}}
\examples{
t<-seq(0,3,by=0.1)
rate<-c(0.6,0.3)
tchange<-c(0,1.75)
pwefun<-pwe(t=t,rate=rate,tchange=tchange)
pwefun
}
\keyword{piecewise exponetial}
