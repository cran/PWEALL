\name{rmstpower}
\alias{rmstpower}
\title{Calculate powers at different cut-points based on difference of restricted mean survival times (RMST)}
\description{A function to calculate powers at different cut-points based on difference of restricted mean survival times (RMST)
account for delayed treatment, discontinued treatment and non-uniform entry}
\usage{
rmstpower(tcut=2,tstudy=seq(tcut,tcut+2,by=0.5),alpha=0.05,twosided=1,
          taur=1.2,u=c(1/taur,1/taur),ut=c(taur/2,taur),pi1=0.5,
          rate11=c(1,0.5),rate21=rate11,rate31=c(0.7,0.4),
          rate41=rate21,rate51=rate21,ratec1=c(0.5,0.6),
          rate10=rate11,rate20=rate10,rate30=rate31,
          rate40=rate20,rate50=rate20,ratec0=c(0.6,0.5),
          tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,
          eps=1.0e-2,veps=1.0e-2,n=1000)
}
\arguments{
  \item{tcut}{timepoint at which rmst is calculated}
  \item{tstudy}{a vector of study time points, which must be not smaller than \code{tcut}}
  \item{alpha}{type-1 error rate}
  \item{twosided}{twosided test=1 or not}
  \item{taur}{Recruitment time}
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
  \item{rp21}{re-randomization prob for the treatment group}
  \item{rp20}{re-randomization prob for the control group}
  \item{eps}{error tolerence}
  \item{veps}{error tolenrence for calculating variance}
  \item{n}{total number of subjects, both groups combined}
}
\details{
The hazard functions corresponding to \code{rate11},...,\code{rate51},\code{ratec1}, \code{rate10},...,\code{rate50},\code{ratec0} are all piecewise constant function taking the form  \eqn{\lambda(t)=\sum_{j=1}^m \lambda_j I(t_{j-1}\le t<t_j)}, where \eqn{\lambda_1,\ldots,\lambda_m} are the corresponding elements of the rates and \eqn{t_0,\ldots,t_{m-1}} are the corresponding elements of \cite{tchange}, \eqn{t_m=\infty}. Note that all the rates must have the same \code{tchange}.
}
\value{
\item{power}{power}
\item{rmst1}{rmst in the treatment group}
\item{se1}{standard error of the rmst in the treatment group}
\item{rmst0}{rmst in the control group}
\item{se0}{standard error of the rmst in the control group}
\item{drmst}{\code{rmst1}-\code{rmst0}}
\item{sed}{standard error of the mean difference}
}
\references{
Luo, et al. (2017)
}
\author{Xiaodong Luo}
\note{This calculates the restricted mean survival times between the treatment and control groups and their standard errors}
\examples{
tcut<-3.0
tstudy<-seq(3,6,by=1)
taur<-1.2
u<-c(1/taur,1/taur)
ut<-c(taur/2,taur)
r11<-c(0.2,0.1)
r21<-r11
r31<-c(0.03,0.02)
r41<-r51<-r21
rc1<-c(0.01,0.02)
r10<-c(0.2,0.2)
r20<-r10
r30<-c(0.02,0.01)
r40<-r50<-r20
rc0<-c(0.02,0.01)
getrmst<-rmstpower(tcut=tcut,tstudy=tstudy,alpha=0.05,twosided=1,
          taur=taur,u=u,ut=ut,pi1=0.5,
          rate11=r11,rate21=r21,rate31=r31,rate41=r41,rate51=r51,ratec1=rc1,
          rate10=r10,rate20=r20,rate30=r30,rate40=r40,rate50=r50,ratec0=rc0,
          tchange=c(0,1),type1=1,type0=1,rp21=0.5,rp20=0.5,n=1000)
###powers at each time point
cbind(tstudy,getrmst$power)
}
\keyword{piecewise exponential distribution}
\keyword{crossover effect}
\keyword{delayed treatment effect}
\keyword{restricted mean survival time}
\keyword{mean difference}
\keyword{variance}% __ONLY ONE__ keyword per line
