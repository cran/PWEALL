\name{rpwu}
\alias{rpwu}
\title{Piecewise uniform distribution: random number generation}
\description{This will generate random numbers according to the specified piecewise uniform distribution}
\usage{
rpwu(nr=10,u=c(0,6,0.4),ut=c(1,2))
}
\arguments{
  \item{nr}{number of random numbers to be generated}
  \item{u}{piecewise constant density}
  \item{ut}{a strictly increasing sequence of time points defining the pieces. The first element must be strictly greater than zero. \code{u} and \code{ut} must have the same length.}
}
\details{Let \eqn{f(t)=\sum_{j=1}^m u_j I(t_{j-1}<t\le t_j)} be the density function, where \eqn{u_1,\ldots,u_m} are the corresponding elements of \cite{u} and \eqn{t_1,\ldots,t_{m}} are the corresponding elements of \cite{ut} and \eqn{t_0=0}.
The distribution function \deqn{F(t)=\sum_{j=1}^m u_j(t\wedge t_j-t\wedge t_{j-1}). }
User must make sure that \eqn{\sum_{j=1}^m u_j (t_j-t_{j-1})=1} before using this function.
}
\value{\item{r}{random numbers}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This provides a random number generator of the piecewise uniform distribution}
\seealso{\code{\link{rpwe}}}
\examples{
nr<-10
u<-c(0.6,0.4)
ut<-c(1,2)
pwur<-rpwu(nr=nr,u=u,ut=ut)
pwur
}
\keyword{piecewise uniform}
\keyword{random number generator}% __ONLY ONE__ keyword per line
