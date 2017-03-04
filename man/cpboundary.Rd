\name{cpboundary}
\alias{cpboundary}
\title{The stopping boundary based on the conditional power criteria}
\description{This will calculate the stopping boundary based on the conditional power criteria, i.e. if observed HR is above the boundary, the conditional power will be lower than the designated level. All the calculation is based on the proportional hazards assumption and the Cox model.}
\usage{
cpboundary(Dplan=300,alpha=0.05,two.sided=1,pi1=0.5,cpcut=c(0.2,0.3,0.4),
            BetaD=log(0.8),Beta0=log(1),prop=seq(0.1,0.9,by=0.1))
}
\arguments{
  \item{Dplan}{Planned number of events at study end}
  \item{alpha}{Type 1 error rate}
  \item{two.sided}{=1 two-sided test and =0 one-sided test}
  \item{pi1}{Allocation probability for the treatment group}
  \item{cpcut}{the designated conditional power level}
  \item{BetaD}{designed log hazard ratio, i.e. under alternative hypothesis}
  \item{Beta0}{null log hazard ratio, i.e. under null hypothesis}
  \item{prop}{proportion of \code{Dplan} observed}
}
\details{This will calculate the stopping boundary based on the conditional power criteria, i.e. if observed HR is above the boundary, the conditional power will be lower than the designated level. All the calculation is based on the proportional hazards assumption and the Cox model.
}
\value{
   \item{CPTbound}{Boundary based on the conditional power under current trend}
   \item{CPNbound}{Boundary based on the conditional power under null hypothesis}
   \item{CPDbound}{Boundary based on the conditional power according to design, i.e. under alternative hypothesis}
}
\references{
Halperin, Lan, Ware, Johnson and DeMets (1982). Controlled Clinical Trials.
}
\author{Xiaodong Luo}
\note{This will calculate the stopping boundary based on the conditional power criteria}
\seealso{\code{\link{cp}},\code{\link{cpstop}}}
\examples{
###Calculate the stopping boundary at 10-90 percent of the target 300 events 
###when the condition power are c(0.2,0.3,0.4) with 
###2:1 allocation ratio between the treatment group and the control group
cpboundary(pi1=2/3)
}
\keyword{conditional power}
\keyword{stopping boundary}


