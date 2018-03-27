\name{rmstutil}
\alias{rmstutil}
\title{A utility function to calculate the true restricted mean survival time (RMST) and its variance
     account for delayed treatment, discontinued treatment and non-uniform entry}
\description{A utility function to calculate the true restricted mean survival time (RMST) and its variance
     account for delayed treatment, discontinued treatment and non-uniform entry}
\usage{
rmstutil(tcut=2.0,tstudy=5.0,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),
        rate1=c(1,0.5),rate2=rate1,rate3=c(0.7,0.4),
        rate4=rate2,rate5=rate2,ratec=c(0.5,0.6),
        tchange=c(0,1),type=1,rp2=0.5,
        eps=1.0e-2,veps=1.0e-2)
}
\arguments{
  \item{tcut}{time point at which rmst is calculated}
  \item{tstudy}{the study time point from first patient in, it must be not smaller than \code{tcut}.}
  \item{taur}{Recruitment time}
  \item{u}{Piecewise constant recuitment rate}
  \item{ut}{Recruitment intervals}
  \item{rate1}{piecewise constant event rate before crossover}
  \item{rate2}{piecewise constant event rate after crossover}
  \item{rate3}{piecewise constant event rate for crossover}
  \item{rate4}{additional piecewise constant event rate for more complex crossover}
  \item{rate5}{additional piecewise constant event rate for more complex crossover}
  \item{ratec}{Hazard for time to censoring}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates \code{rate1} to \code{ratec} and tchange must have the same length.}
  \item{type}{type of crossover, 1=markov, 2=semi-markov, 3=hybrid}
  \item{rp2}{re-randomization probability to receive the rescue treatment when semi-markov crossover occurs.
        When it happens, the overall hazard will be rp2*r2(t-s)+(1-rp2)*r4(t), where r2 is the hazard 
             for the semi-markov rescue treatment and r4 is hazard for the markov rescue treatment.} 
 \item{eps}{A small number representing the error tolerance when calculating the utility function
              \deqn{\Phi_l(x)=\frac{\int_0^x s^l e^{-s}ds}{x^{l+1}}} with \eqn{l=0,1,2}.}
\item{veps}{A small number representing the error tolerance when calculating the variance.}
}
\details{More details}
\value{
  \item{tcut}{time point at which rmst is calculated}
  \item{tstudy}{the study time point from first patient in, it must be not smaller than \code{tcut}}
  \item{rmst}{rmst at cut-point \code{tcut}}
  \item{var}{the variance of \code{rmst}}
  \item{vadd}{the additional variance term of \code{rmst}}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This calculates the "true" variance of restricted mean survival times}
\examples{
r1<-c(0.6,0.3)
r2<-c(0.6,0.6)
r3<-c(0.1,0.2)
r4<-c(0.5,0.4)
r5<-c(0.4,0.5)
rc<-c(0.1,0.1)
rmt<-rmstutil(tcut=2.0,tstudy=5.0,taur=5,
        rate1=r1,rate2=r2,rate3=r3,
        rate4=r4,rate5=r5,ratec=rc,
        tchange=c(0,1),type=1,rp2=0.5,
        eps=1.0e-2,veps=1.0e-2)
rmt
}
\keyword{piecewise exponential}
\keyword{treatment crossover}
\keyword{restricted mean survival times}
\keyword{variance}% __ONLY ONE__ keyword per line
