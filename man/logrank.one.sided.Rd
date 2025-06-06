% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logrank.one.sided.R
\name{logrank.one.sided}
\alias{logrank.one.sided}
\title{Perform one-sided logrank test}
\usage{
logrank.one.sided(time, event, group, STRATA = NULL)
}
\arguments{
\item{time}{survival time}

\item{event}{event status (0 = censor, 1 = event)}

\item{group}{group indicator (0 = control, 1 = experimental arm)}

\item{STRATA}{strata variable for stratified log-rank test}
}
\value{
\describe{
\item{z}{Test statistics z value}
\item{p}{one sided p value}
}
}
\description{
This functions performs the one-sided logrank test. The standard logrank test in survival package
only produces two-sided test. This function can facilitate one-sided logrank test. Larger z statistic indicates
better treatment effect.
}
\examples{

n <- 100
time <- c(rexp(n, rate=log(2)/12), rexp(n, rate=log(2)/12*1.2))
event <- rep(1, n*2)
group <- c(rep(0, n), rep(1, n))
STRATA <- rep(c(0,1), n)

logrank.one.sided(time, event, group, STRATA)

}
