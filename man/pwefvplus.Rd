\name{pwefvplus}
\alias{pwefvplus}
\title{A utility functon}
\description{This will calculate the more complex integration accounting for crossover}
\usage{
pwefvplus(t=seq(0,5,by=0.5),rate1=c(0,5,0.8),rate2=rate1,
                   rate3=c(0.1,0.2),rate4=rate2,rate5=rate2,
                   rate6=c(0.5,0.3),tchange=c(0,3),type=1,
                   rp2=0.5,eps=1.0e-2)
}
\arguments{
  \item{t}{A vector of time points}
  \item{rate1}{piecewise constant event rate}
  \item{rate2}{piecewise constant event rate}
  \item{rate3}{piecewise constant event rate}
  \item{rate4}{additional piecewise constant}
  \item{rate5}{additional piecewise constant}
  \item{rate6}{piecewise constant event rate for censoring}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which event rate changes. The first element of tchange must be zero. The above rates and tchange must have the same length.}
  \item{type}{type of the crossover, markov, semi-markov and hybrid}
  \item{rp2}{re-randomization prob}
  \item{eps}{tolerance}
}
\details{Let \eqn{h_1,\ldots,h_6} correspond to \code{rate1},...,\code{rate6}, and \eqn{H_1,\ldots,H_6} be the corresponding survival functions. Also let \eqn{\pi_2=\code{rp2}}.
when \code{type}=1, we calculate
\deqn{\int_0^t s^k h_2(s)H_2(s)H_6(s)\int_0^s h_3(u)H_1(u)H_3(u)/H_2(u)duds;}
when \code{type}=2, we calculate
\deqn{\int_0^t s^kH_6(s)\int_0^s h_3(u)H_1(u)H_3(u)h_2(s-u)H_2(s-u)duds;}
when \code{type}=3, we calculate the sum of 
\deqn{\pi_2\int_0^t s^kH_4^{1-\pi_2}(s)H_6(s)\int_0^s h_3(u)H_1(u)H_3(u)h_2(s-u)H_2^{\pi_2}(s-u)/H_4^{1-\pi_2}(u)duds}
and 
\deqn{(1-\pi_2)\int_0^t s^kh_4(s)H_4^{1-\pi_2}(s)H_6(s)\int_0^s h_3(u)H_1(u)H_3(u)H_2^{\pi_2}(s-u)/H_4^{1-\pi_2}(u)duds;}
when \code{type}=4, we calculate the sum of 
\deqn{\pi_2\int_0^t s^kH_6(s)\int_0^s h_3(u)H_1(u)H_3(u)h_2(s-u)H_2(s-u)duds}
and
\deqn{(1-\pi_2)\int_0^t s^kh_4(s)H_4(s)H_6(s)\int_0^s h_3(u)H_1(u)H_3(u)/H_4(u)duds;}
when \code{type}=5, we calculate the sum of
\deqn{\pi_2\int_0^t s^kH_6(s)\int_0^s h_3(u)H_1(u)H_3(u)h_2(s-u)H_2(s-u)duds}
and
\deqn{(1-\pi_2)\int_0^t s^kH_6(s)\int_0^s h_3(u)H_1(u)H_3(u)h_4(s-u)H_4(s-u)duds.}
}
\value{\item{f0}{values when \eqn{k=0}}
       \item{f1}{values when \eqn{k=1}}
       \item{f2}{values when \eqn{k=2}}
}
\references{
Luo et al. (2018) Design and monitoring of survival trials in complex scenarios, Statistics in Medicine <doi: https://doi.org/10.1002/sim.7975>.
}
\author{Xiaodong Luo}
\note{This provides the result of the complex integration}
\seealso{\code{\link{rpwe}}}
\examples{
r1<-c(0.6,0.3)
r2<-c(0.6,0.6)
r3<-c(0.1,0.2)
r4<-c(0.5,0.4)
r5<-c(0.4,0.5)
r6<-c(0.4,0.5)
tchange<-c(0,1.75)
pwefun<-pwefvplus(t=seq(0,5,by=0.5),rate1=r1,rate2=r2,rate3=r3,
                 rate4=r4,rate5=r5,rate6=r6,
                 tchange=c(0,3),type=1,eps=1.0e-2)
pwefun
}
\keyword{piecewise exponential}
\keyword{treatment crossover}

