\name{pwecxpwu}
\alias{pwecxpwu}
\title{Integration of the density of piecewise exponential distribution with crossover effect, censoring and recruitment function}
\description{This will calculate the functions according to the piecewise exponential distribution with crossover}
\usage{
pwecxpwu(t=seq(0,10,by=0.5),taur=5,
        u=c(1/taur,1/taur),ut=c(taur/2,taur),
        rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
        rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
        tchange=c(0,1),type=1,rp2=0.5,eps=1.0e-2)
}
\arguments{
  \item{t}{a vector of time points}
  \item{taur}{recruitment time}
  \item{u}{recruitment rate}
  \item{ut}{recruitment interval, must have the same length as \code{u}}
  \item{rate1}{piecewise constant event rate before crossover}
  \item{rate2}{piecewise constant event rate after crossover}
  \item{rate3}{piecewise constant event rate for crossover}
  \item{rate4}{additional piecewise constant event rate for more complex crossover}
  \item{rate5}{additional piecewise constant event rate for more complex crossover}
  \item{ratec}{censoring piecewise constant event rate}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates \code{rate1} to \code{ratec} and tchange must have the same length.}
  \item{type}{type of crossover, i.e. markov, semi-markov and hybrid}
  \item{rp2}{re-randomization prob}
  \item{eps}{tolerance}
}
\details{This is to calculate the function (and its derivative)
\deqn{\xi(t)=\int_0^t G_E(t-s)\widetilde{f}(s)S_C(s)ds,}
where \eqn{G_E} is the accrual function defined by \code{taur}, \code{u} and \code{ut}, \eqn{S_C} is the piecewise exponential survival function of the censoring time, defined by \code{tchange} and \code{ratec}, and \eqn{\widetilde{f}} is the density for the event distribution subject to crossover defined by \code{tchange}, \code{rate1} to \code{rate5} and \code{type}.
}
\value{\item{du}{the function}
       \item{duprime}{its derivative}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\seealso{\code{\link{rpwe}}}
\examples{
taur<-2
u<-c(0.6,0.4)
ut<-c(1,2)
r1<-c(0.6,0.3)
r2<-c(0.6,0.6)
r3<-c(0.1,0.2)
r4<-c(0.5,0.4)
r5<-c(0.4,0.5)
rc<-c(0.5,0.6)
exu<-pwecxpwu(t=seq(0,10,by=0.5),taur=taur,u=u,ut=ut,
        rate1=r1,rate2=r2,rate3=r3,rate4=r4,rate5=r5,ratec=rc,
        tchange=c(0,1),type=1,eps=1.0e-2)
c(exu$du,exu$duprime)
}
\keyword{piecewise exponential}
\keyword{treatment crossover}
\keyword{piecewise uniform}
\keyword{random number generator}
