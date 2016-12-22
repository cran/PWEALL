\name{pwu}
\alias{pwu}
\title{Piecewise uniform distribution: distribution}
\description{This will calculate the distribution function of the piecewise uniform distribution}
\usage{
pwu(t=seq(0,1,by=0.1),u=c(0,5,0.5),ut=c(1,2))
}
\arguments{
  \item{t}{a vector of time points}
  \item{u}{piecewise constant density}
  \item{ut}{a strictly increasing sequence of time points defining the pieces. The first element must be strictly greater than zero. \code{u} and \code{ut} must have the same length.}
}
\details{Let \eqn{f(t)=\sum_{j=1}^m u_j I(t_{j-1}<t\le t_j)} be the density function, where \eqn{u_1,\ldots,u_m} are the corresponding elements of \cite{u} and \eqn{t_1,\ldots,t_{m}} are the corresponding elements of \cite{ut} and \eqn{t_0=0}.
The distribution function \deqn{F(t)=\sum_{j=1}^m u_j(t\wedge t_j-t\wedge t_{j-1}). }
User must make sure that \eqn{\sum_{j=1}^m u_j (t_j-t_{j-1})=1} before using this function.
}
\value{\item{dist}{distribution}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This provides distribution of the piecewise uniform distribution}
\seealso{\code{\link{pwe}}}
\examples{
t<-seq(-1,3,by=0.5)
u<-c(0.6,0.4)
ut<-c(1,2)
pwud<-pwu(t=t,u=u,ut=ut)
pwud
}
\keyword{piecewise uniform}
\keyword{distribution}

