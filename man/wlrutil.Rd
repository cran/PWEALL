\name{wlrutil}
\alias{wlrutil}
\title{A utility function to calculate some common functions in contructing weights}
\description{A utility function to calculate some common functions in contructing weights}
\usage{
wlrutil(y=c(1,2,3),d=c(1,0,1),z=c(1,0,0),te=c(1,3),eps=1.0e-08)
}
\arguments{
  \item{y}{observed times}
  \item{d}{non-censoring indicators}
  \item{z}{group indicators with \code{z}=1 treatment and \code{z}=0 control}
  \item{te}{(ascendingly) ordered unique event times from both groups}
  \item{eps}{tolerence when comparing event times}
}
\details{More details}
\value{
  \item{mfunc}{various functions in column}
}
\author{Xiaodong Luo}
\examples{
ww<-wlrutil(y=c(1,2,3),d=c(1,0,1),z=c(1,0,0),te=c(1,3),eps=1.0e-08)
ww
}
\keyword{weighted log-rank}
