\name{spf}
\alias{spf}
\title{A utility function}
\description{A utility function to calculate a ratio.}
\usage{
spf(x=seq(-1,1,by=0.2),eps=1.0e-3)
}
\arguments{
  \item{x}{A vector}
  \item{eps}{tolerance}
}
\details{This is to calculate
\deqn{\Phi_l(x)=\frac{\int_0^x s^le^{-s}ds}{x^{l+1}},\hspace{0.5cm}l=0,1,2.}
This function is well defined even when x=0. However, it is numerical chanllenging to calculate it when x is small. So when
\eqn{|x|\le \code{eps}} we approximate this function and the absolute error is \eqn{\code{eps}^5}.
}
\value{\item{fx1}{when \eqn{l=0};}
       \item{fx2}{when \eqn{l=1};}
       \item{fx3}{when \eqn{l=2}.}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{Version 1.0 (7/19/2016)}
\examples{
fun<-spf(x=seq(-1,1,by=0.2),eps=1.0e-3)
fun
}
\keyword{utility function}
