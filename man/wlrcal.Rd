\name{wlrcal}
\alias{wlrcal}
\title{A utility function to calculate the weighted log-rank statistics and their varainces given the weights}
\description{A utility function to calculate the weighted log-rank statistics and their varainces given the weights}
\usage{
wlrcal(n=10,te=c(1,2,3),tfix=2.0,dd1=c(1,0,1),dd0=c(0,1,0),r1=c(1,2,3),r0=c(1,2,3),
       weights=matrix(1,nrow=length(te),ncol=1),eps=1.0e-08)
}
\arguments{
  \item{n}{total number of subjects in the study}
  \item{te}{(ascendingly) ordered unique event times from both groups}
  \item{tfix}{time point where weighted log-rank is calcualted}
  \item{dd1}{number of events from treatment group at each \code{te}}
  \item{dd0}{number of events from control group at each \code{te}}
  \item{r1}{number of at-risk subjects from treatment group at each \code{te}}
  \item{r0}{number of at-risk subjects from control group at each \code{te}}
  \item{weights}{user specified weights, each column is a set of weights at each \code{te}}
  \item{eps}{tolerence when comparing event times}
}
\details{More details}
\value{
  \item{test}{unscaled test statistics}
  \item{var}{variances of the unsclaed test statistics}
  \item{wlr}{weighted log-rank statistics, i.e. scaled test statsitics}
  \item{wlcor}{the correlation matrix of the weighted log-rank statistics}
}
\author{Xiaodong Luo}
\examples{
lr<-wlrcal(n=10,te=c(1,2,3),tfix=2.0,dd1=c(1,0,1),dd0=c(0,1,0),r1=c(1,2,3),r0=c(1,2,3))
lr
}
\keyword{weighted log-rank}
