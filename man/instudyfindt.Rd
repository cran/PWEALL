\name{instudyfindt}
\alias{instudyfindt}
\title{calculate the timeline in study when some or all subjects have entered}
\description{This will calculate the timeline from some timepoint in study when some/all subjects have entered accouting for staggered entry, delayed treatment effect,  treatment crossover and loss to follow-up.}
\usage{
instudyfindt(target=400,y=exp(rnorm(300)),z=rbinom(300,1,0.5),
                  d=rep(c(0,1,2),each=100),
                  tcut=2,blinded=1,type0=1,type1=type0,
                  rp20=0.5,rp21=0.5,tchange=c(0,1),
                  rate10=c(1,0.7),rate20=c(0.9,0.7),rate30=c(0.4,0.6),rate40=rate20,
                  rate50=rate20,ratec0=c(0.3,0.3),
                  rate11=rate10,rate21=rate20,rate31=rate30,
                  rate41=rate40,rate51=rate50,ratec1=ratec0,
                  withmorerec=1,
                  ntotal=1000,taur=5,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
                  ntype0=1,ntype1=1,
                  nrp20=0.5,nrp21=0.5,ntchange=c(0,1),
                  nrate10=rate10,nrate20=rate20,nrate30=rate30,nrate40=rate40,
                  nrate50=rate50,nratec0=ratec0,
                  nrate11=rate10,nrate21=rate20,nrate31=rate30,nrate41=rate40,
                  nrate51=rate50,nratec1=ratec0,
                  eps=1.0e-2,init=tcut*1.1,epsilon=0.001,maxiter=100)
}
\arguments{
  \item{target}{target number of events}
  \item{y}{observed times}
  \item{z}{observed treatment indicator when \code{blinded}=0, z=1 denotes the treatment group and 0 the control group}
  \item{d}{event indicator, 1=event, 0=censored, 2=no event or censored up to \code{tcut}, the data cut-point}
  \item{tcut}{the data cut-point}
  \item{blinded}{blinded=1 if the data is blinded,=0 if it is unblinded}
  \item{type0}{type of the crossover for the observed data in the control group}
  \item{type1}{type of the crossover for the observed data in the treatment group}
  \item{rp20}{re-randomization prob for the observed data in the control group}
  \item{rp21}{re-randomization prob for the observed data in the treatment group}
  \item{tchange}{A strictly increasing sequence of time points at which the event rates changes. The first element of tchange                 must be zero. It must have the same length as \code{ratejk}, j=1,2,3,4,5,c; k=0,1}
  \item{rate10}{Hazard before crossover for the old subjects in the control group}
  \item{rate20}{Hazard after crossover for the old subjects in the control group}
  \item{rate30}{Hazard for time to crossover for the old subjects in the control group}
  \item{rate40}{Hazard after crossover for the old subjects in the control group for complex case}
  \item{rate50}{Hazard after crossover for the old subjects in the control group for complex case}
  \item{ratec0}{Hazard for time to censoring for the old subjects in the control group}
  \item{rate11}{Hazard before crossover for the old subjects in the treatment group}
  \item{rate21}{Hazard after crossover for the old subjects in the treatment group}
  \item{rate31}{Hazard for time to crossover for the old subjects in the treatment group}
  \item{rate41}{Hazard after crossover for the old subjects in the treatment group for complex case}
  \item{rate51}{Hazard after crossover for the old subjects in the treatment group for complex case}
  \item{ratec1}{Hazard for time to censoring for the old subjects in the treatment group}
  \item{withmorerec}{withmorerec=1 if more subjects are needed to be recruited; =0 otherwise}
  \item{ntotal}{total number of the potential new subjects}
  \item{taur}{recruitment time for the potential new subjects}
  \item{u}{Piecewise constant recuitment rate for the potential new subjects}
  \item{ut}{Recruitment intervals for the potential new subjects}
  \item{pi1}{Allocation probability to the treatment group for the potential new subjects}
  \item{ntype0}{type of the crossover for the potential new subjects in the control group}
  \item{ntype1}{type of the crossover for the potential new subjects in the treatment group}
  \item{nrp20}{re-randomization prob for the potential new subjects in the control group}
  \item{nrp21}{re-randomization prob for the potential new subjects in the treatment group}
  \item{ntchange}{A strictly increasing sequence of time points at which the event rates changes. The first element of ntchange must be zero. It must have the same length as \code{nratejk}, j=1,2,3,4,5,c; k=0,1}
  \item{nrate10}{Hazard before crossover for the potential new subjects in the control group}
  \item{nrate20}{Hazard after crossover for the potential new subjects in the control group}
  \item{nrate30}{Hazard for time to crossover for the potential new subjects in the control group}
  \item{nrate40}{Hazard after crossover for the potential new subjects in the control group for complex case}
  \item{nrate50}{Hazard after crossover for the potential new subjects in the control group for complex case}
  \item{nratec0}{Hazard for time to censoring for the potential new subjects in the control group}
  \item{nrate11}{Hazard before crossover for the potential new subjects in the treatment group}
  \item{nrate21}{Hazard after crossover for the potential new subjects in the treatment group}
  \item{nrate31}{Hazard for time to crossover for the potential new subjects in the treatment group}
  \item{nrate41}{Hazard after crossover for the potential new subjects in the treatment group for complex case}
  \item{nrate51}{Hazard after crossover for the potential new subjects in the treatment group for complex case}
  \item{nratec1}{Hazard for time to censoring for the potential new subjects in the treatment group}
  \item{eps}{A small number representing the error tolerance when calculating the utility function
              \deqn{\Phi_l(x)=\frac{\int_0^x s^l e^{-s}ds}{x^{l+1}}} with \eqn{l=0,1,2}.}
  \item{init}{initital value of the timeline estimate}
  \item{epsilon}{A small number representing the error tolerance when calculating the timeline.}
  \item{maxiter}{Maximum number of iterations when calculating the timeline}
}
\details{
The hazard functions corresponding to \code{rate11},...,\code{rate51},\code{ratec1}, \code{rate10},...,\code{rate50},\code{ratec0} are all piecewise constant function taking the form  \eqn{\lambda(t)=\sum_{j=1}^m \lambda_j I(t_{j-1}\le t<t_j)}, where \eqn{\lambda_1,\ldots,\lambda_m} are the corresponding elements of the rates and \eqn{t_0,\ldots,t_{m-1}} are the corresponding elements of \cite{tchange}, \eqn{t_m=\infty}. Note that all the rates must have the same \code{tchange}.
The hazard functions corresponding to \code{nrate11},...,\code{nrate51},\code{nratec1}, \code{nrate10},...,\code{nrate50},\code{nratec0} are all piecewise constant functions and all must have the same \code{ntchange}.
}
\value{\item{t1}{the calculated timeline}
       \item{dvalue}{the number of events}
       \item{dvprime}{the derivative of the event cummulative function at time \code{t1}}
       \item{tvar}{the variance of the timeline estimator}
       \item{ny}{total number of subjects that could be in the study}
       \item{eps}{final tolerance}
       \item{iter}{Number of iterations performed}
       \item{t1hist}{the history of the iteration for timeline}
       \item{dvaluehist}{the history of the iteration for the event count}
       \item{dvprimehist}{the history of the iteration for the derivative of event count with respect to time}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{Version 1.0 (7/19/2016)}
\seealso{\code{\link{pwe}},\code{\link{rpwe}},\code{\link{qpwe}},\code{\link{pwecxpwufindt}}}
\examples{
n<-1000
target<-550
ntotal<-1000
pi1<-0.5
taur<-2.8
u<-c(1/taur,1/taur)
ut<-c(taur/2,taur)
r11<-c(1,0.5)
r21<-c(0.5,0.8)
r31<-c(0.7,0.4)
r41<-r51<-r21
rc1<-c(0.5,0.6)
r10<-c(1,0.7)
r20<-c(0.5,1)
r30<-c(0.3,0.4)
r40<-r50<-r20
rc0<-c(0.2,0.4)
tchange<-c(0,1.873)
tcut<-2

####generate the data
E<-T<-C<-Z<-delta<-rep(0,n)
E<-rpwu(nr=n,u=u,ut=ut)$r
Z<-rbinom(n,1,pi1)
n1<-sum(Z)
n0<-sum(1-Z)
C[Z==1]<-rpwe(nr=n1,rate=rc1,tchange=tchange)$r
C[Z==0]<-rpwe(nr=n0,rate=rc0,tchange=tchange)$r
T[Z==1]<-rpwecx(nr=n1,rate1=r11,rate2=r21,rate3=r31,
                rate4=r41,rate5=r51,tchange=tchange,type=1)$r
T[Z==0]<-rpwecx(nr=n0,rate1=r10,rate2=r20,rate3=r30,
                rate4=r40,rate5=r50,tchange=tchange,type=1)$r
y<-pmin(pmin(T,C),tcut-E)
y1<-pmin(C,tcut-E)
delta[T<=y]<-1
delta[C<=y]<-0
delta[tcut-E<=y & tcut-E>0]<-2
delta[tcut-E<=y & tcut-E<=0]<--1

ys<-y[delta>-1]
Zs<-Z[delta>-1]
ds<-delta[delta>-1]

nplus<-sum(delta==-1)
nd0<-sum(ds==0)
nd1<-sum(ds==1)
nd2<-sum(ds==2)


ntaur<-taur-tcut
nu<-c(1/ntaur,1/ntaur)
nut<-c(ntaur/2,ntaur)

###calculate the timeline at baseline
xt<-pwecxpwufindt(target=target,ntotal=n,taur=taur,u=u,ut=ut,pi1=pi1,
              rate11=r11,rate21=r21,rate31=r31,ratec1=rc1,
              rate10=r10,rate20=r20,rate30=r30,ratec0=rc0,
              tchange=tchange,eps=0.001,init=taur,epsilon=0.000001,maxiter=100)
###calculate the timeline in study
yt<-instudyfindt(target=target,y=ys,z=Zs,d=ds,
                       tcut=tcut,blinded=0,type1=1,type0=1,tchange=tchange,
                       rate10=r10,rate20=r20,rate30=r30,ratec0=rc0,
                       rate11=r11,rate21=r21,rate31=r31,ratec1=rc1,
                       withmorerec=1,
                       ntotal=nplus,taur=ntaur,u=nu,ut=nut,pi1=pi1,
                       ntype1=1,ntype0=1,ntchange=tchange,
                       nrate10=r10,nrate20=r20,nrate30=r30,nratec0=rc0,
                       nrate11=r11,nrate21=r21,nrate31=r31,nratec1=rc1,
                       eps=1.0e-2,init=2,epsilon=0.001,maxiter=100)
##timelines                       
c(yt$t1,xt$t1)
##standard errors of the timeline estimators 
c(sqrt(yt$tvar/yt$ny),sqrt(xt$tvar/n))
###95 percent CIs
c(yt$t1-1.96*sqrt(yt$tvar/yt$ny),yt$t1+1.96*sqrt(yt$tvar/yt$ny))
c(xt$t1-1.96*sqrt(xt$tvar/n),xt$t1+1.96*sqrt(xt$tvar/n))

}
\keyword{piecewise exponential}
\keyword{timeline}
\keyword{piecewise uniform}
\keyword{delayed treatment effect}
\keyword{treatment crossover}

