\name{PWEALL-package}
\alias{PWEALL-package}
\alias{PWEALL}
\docType{package}
\title{
\packageTitle{PWEALL}
}
\description{
\packageDescription{PWEALL}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{PWEALL}
\packageIndices{PWEALL}

There are 5 types of crossover considered in the package: (1) Markov crossover, (2) Semi-Markov crosover, (3) Hybrid crossover-1, (4) Hybrid crossover-2 and (5) Hybrid crossover-3. The first 3 types are described in Luo et al. (2018). The fourth and fifth types are added for Version 1.3.0. The crossover type is determined by the hazard function after crossover \eqn{\lambda_2^{\bf x}(t\mid u)}. For Type (1), the Markov crossover,
\deqn{\lambda_2^{\bf x}(t\mid u)=\lambda_2(t).} 
For Type (2), the Semi-Markov crossover,
\deqn{\lambda_2^{\bf x}(t\mid u)=\lambda_2(t-u).}
For Type (3), the hybrid crossover-1,
\deqn{\lambda_2^{\bf x}(t\mid u)=\pi_2\lambda_2(t-u)+(1-\pi_2)\lambda_4(t).}
For Type (4), the hazard after crossover is 
\deqn{\lambda_2^{\bf x}(t\mid u)=\frac{\pi_2\lambda_2(t-u)S_2(t-u)+(1-\pi_2)\lambda_4(t)S_4(t)/S_4(u)}{\pi_2 S_2(t-u)+(1-\pi_2)S_4(t)/S_4(u)}.}
For Type (5), the hazard after crossover is 
\deqn{\lambda_2^{\bf x}(t\mid u)=\frac{\pi_2\lambda_2(t-u)S_2(t-u)+(1-\pi_2)\lambda_4(t-u)S_4(t-u)}{\pi_2 S_2(t-u)+(1-\pi_2)S_4(t-u)}.}
The types (4) and (5) are more closely related to "re-randomization", i.e. when a patient crosses, (s)he will have probability \eqn{\pi_2} to have hazard \eqn{\lambda_2} and probability \eqn{1-\pi_2} to have hazard \eqn{\lambda_4}. The types (4) and (5) differ in having \eqn{\lambda_4} as Markov or Semi-markov.  
}
\author{
\packageAuthor{PWEALL}

Maintainer: \packageMaintainer{PWEALL}
}
\references{
Luo et al. (2018) Design and monitoring of survival trials in complex scenarios, Statistics in Medicine <doi: https://doi.org/10.1002/sim.7975>.
}
\keyword{piecewise exponential}
\keyword{various functions}

