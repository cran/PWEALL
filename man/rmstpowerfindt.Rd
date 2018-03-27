\name{rmstpowerfindt}
\alias{rmstpowerfindt}
\title{Calculating the timepoint where a certain power of mean difference of RMSTs is obtained}
\description{This will calculate the timepoint where a certain power of the mean difference of RMSTs  is obtained accouting for staggered entry, delayed treatment effect,  treatment crossover and loss to follow-up.}
\usage{
rmstpowerfindt(power=0.9,alpha=0.05,twosided=1,tcut=2,tupp=5,tlow=3.0,taur=1.2,
           u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
           rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
           rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
           rate10=rate11,rate20=rate10,rate30=rate31,
           rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
           tchange=c(0,1),type1=1,type0=1,
           rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,
           n=1000,maxiter=20,itereps=0.001)
}
\arguments{
  \item{power}{the desired power}
  \item{alpha}{type-1 error}
  \item{twosided}{twoside test or not}
  \item{tcut}{time point at which rmst is calculated}
  \item{tupp}{an upper study time point where the power should be larger than \code{power}}
  \item{tlow}{a lower study time point where the power should be smaller than \code{power}, \code{tlow} must be not smaller than \code{tcut}}
  \item{taur}{recruitment time}
  \item{u}{Piecewise constant recuitment rate}
  \item{ut}{Recruitment intervals}
  \item{pi1}{Allocation probability for the treatment group}
  \item{rate11}{Hazard before crossover for the treatment group}
  \item{rate21}{Hazard after crossover for the treatment group}
  \item{rate31}{Hazard for time to crossover for the treatment group}
  \item{rate41}{Hazard after crossover for the treatment group for complex case}
  \item{rate51}{Hazard after crossover for the treatment group for complex case}
  \item{ratec1}{Hazard for time to censoring for the treatment group}
  \item{rate10}{Hazard before crossover for the control group}
  \item{rate20}{Hazard after crossover for the control group}
  \item{rate30}{Hazard for time to crossover for the control group}
  \item{rate40}{Hazard after crossover for the control group for complex case}
  \item{rate50}{Hazard after crossover for the control group for complex case}
  \item{ratec0}{Hazard for time to censoring for the control group}
  \item{tchange}{A strictly increasing sequence of time points at which the event rates changes. The first element of tchange                 must be zero. It must have the same length as \code{rate11}, \code{rate21}, \code{rate31}, etc.}
  \item{type1}{Type of crossover in the treatment group}
  \item{type0}{Type of crossover in the control group}
  \item{rp21}{re-randomization prob in the treatment group}
  \item{rp20}{re-randomization prob in the control group}
  \item{eps}{error tolerence}
  \item{veps}{error tolenrence for calculating variance}
  \item{n}{total number of subjects}
  \item{maxiter}{maximum number of bi-section iterations}
  \item{itereps}{error tolerance of \code{power}}
}
\details{
The hazard functions corresponding to \code{rate11},...,\code{rate51},\code{ratec1}, \code{rate10},...,\code{rate50},\code{ratec0} are all piecewise constant function taking the form  \eqn{\lambda(t)=\sum_{j=1}^m \lambda_j I(t_{j-1}\le t<t_j)}, where \eqn{\lambda_1,\ldots,\lambda_m} are the corresponding elements of the rates and \eqn{t_0,\ldots,t_{m-1}} are the corresponding elements of \cite{tchange}, \eqn{t_m=\infty}. Note that all the rates must have the same \code{tchange}.
}

\value{
\item{time}{time calculated when the iterations stop}
\item{power}{the power at \code{time}}
\item{err}{distance from the desired power}
\item{k}{number of bi-section iterations performed}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{Version 1.0 (8/8/2017)}
\seealso{\code{\link{pwe}},\code{\link{rpwe}},\code{\link{qpwe}},\code{\link{ovbeta}},\code{\link{innervar}}}
\examples{
tcut<-3.0
tstudy<-seq(3,6,by=0.2)
taur<-2
u<-c(0.3,0.7)
ut<-c(taur/2,taur)
r11<-c(0.2,0.1)
r21<-r11
r31<-c(0.03,0.02)
r41<-r51<-r21
rc1<-c(0.05,0.04)
r10<-c(0.22,0.22)
r20<-r10
r30<-c(0.02,0.01)
r40<-r50<-r20
rc0<-c(0.04,0.05)
ntotal<-1200
getrmst<-rmstpower(tcut=tcut,tstudy=tstudy,alpha=0.05,twosided=1,
        taur=taur,u=u,ut=ut,pi1=0.5,
        rate11=r11,rate21=r21,rate31=r31,rate41=r41,rate51=r51,ratec1=rc1,
        rate10=r10,rate20=r20,rate30=r30,rate40=r40,rate50=r50,ratec0=rc0,
        tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,n=ntotal)
###powers at each time point
cbind(tstudy,getrmst$power)

###90 percent power should be in (3,4)
gettime<-rmstpowerfindt(power=0.9,alpha=0.05,twosided=1,tcut=tcut,tupp=4,tlow=3.0,taur=taur,
          u=u,ut=ut,pi1=0.5,rate11=r11,rate21=r21,rate31=r31,rate41=r41,rate51=r51,ratec1=rc1,
          rate10=r10,rate20=r20,rate30=r30,rate40=r40,rate50=r50,ratec0=rc0,
          tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,eps=1.0e-2,veps=1.0e-2,
          n=ntotal,maxiter=20,itereps=0.0001)
gettime
}
\keyword{piecewise exponential}
\keyword{piecewise uniform}
\keyword{delayed treatment effect}
\keyword{treatment crossover}
\keyword{timeline for certain power}
\keyword{mean difference of RMSTs}
\keyword{power}


