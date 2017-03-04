\name{cp}
\alias{cp}
\title{Conditional power given observed log hazard ratio}
\description{This will calculate the conditional power given the observed log hazard ratio based on Cox model}
\usage{
cp(Dplan=300,alpha=0.05,two.sided=1,pi1=0.5,Obsbeta=log(seq(1,0.6,by=-0.01)),
   BetaD=log(0.8),Beta0=log(1),prop=seq(0.1,0.9,by=0.1))
}
\arguments{
  \item{Dplan}{Planned number of events at study end}
  \item{alpha}{Type 1 error rate}
  \item{two.sided}{=1 two-sided test and =0 one-sided test}
  \item{pi1}{Allocation probability for the treatment group}
  \item{Obsbeta}{observed log hazard ratio}
  \item{BetaD}{designed log hazard ratio, i.e. under alternative hypothesis}
  \item{Beta0}{null log hazard ratio, i.e. under null hypothesis}
  \item{prop}{proportion of \code{Dplan} observed}
}
\details{This is to calculated conditional power at time point when certain percent of target number of event has been observed and an observed log hazard ratio is provided.
}
\value{
   \item{CPT}{Conditional power under current trend}
   \item{CPN}{Conditional power under null hypothesis}
   \item{CPD}{Conditional power according to design, i.e. under alternative hypothesis}
}
\references{
Halperin, Lan, Ware, Johnson and DeMets (1982). Controlled Clinical Trials.
}
\author{Xiaodong Luo}
\note{This will calculate the conditional power given the observed log hazard ratio based on Cox model}
\seealso{\code{\link{cpboundary}},\code{\link{cpstop}}}
\examples{
###Calculate the CP at 10-90 percent of the target 300 events when the observed HR 
###are seq(1,0.6,by=-0.01) with 2:1 allocation 
###ratio between the treatment group and the control group
cp(pi1=2/3)
}
\keyword{conditional power}

