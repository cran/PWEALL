\name{qpwu}
\alias{qpwu}
\title{Piecewise uniform distribution: quantile function}
\description{This will provide the quantile function of the specified piecewise uniform distribution}
\usage{
qpwu(p=seq(0,1,by=0.1),u=c(0,5,0.5),ut=c(1,2))
}
\arguments{
  \item{p}{a vector of probabilities}
  \item{u}{piecewise constant density}
  \item{ut}{time points at which event rate changes. This must be an strictly increasing sequence. \code{ut} and \code{u} must have the same length.}
}
\details{Let \eqn{f(t)=\sum_{j=1}^m u_j I(t_{j-1}<t\le t_j)} be the density function, where \eqn{u_1,\ldots,u_m} are the corresponding elements of \cite{u} and \eqn{t_1,\ldots,t_{m}} are the corresponding elements of \cite{ut} and \eqn{t_0=0}.
The distribution function \deqn{F(t)=\sum_{j=1}^m u_j(t\wedge t_j-t\wedge t_{j-1}). }
User must make sure that \eqn{\sum_{j=1}^m u_j (t_j-t_{j-1})=1} before using this function.}
\value{\item{q}{quantiles}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This provides the quantile function related to the piecewise uniform distribution}
\seealso{piecewise uniform}
\examples{
p<-seq(0,1,by=0.1)
u<-c(0.6,0.4)
ut<-c(1,2)
pwuq<-qpwu(p=p,u=u,ut=ut)
pwuq
}
\keyword{piecewise uniform}
\keyword{quantiles}% __ONLY ONE__ keyword per line
