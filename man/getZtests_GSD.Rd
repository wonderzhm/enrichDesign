% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getZtests_GSD.R
\name{getZtests_GSD}
\alias{getZtests_GSD}
\title{Get group sequential design test statistics for a given trial.}
\usage{
getZtests_GSD(object, alpha = 0.025, HR.Sc.threshold = 1)
}
\arguments{
\item{object}{An object return by \link{getZstats}.}

\item{alpha}{Type I error, always one-sided.}

\item{HR.Sc.threshold}{Hazard ratio threshold for futility.}
}
\value{
It returns a list of test statistics used for later adjustments.
}
\description{
Get group sequential design test statistics for a given trial.
}
\examples{
d <- simu_enrich_trial(n = 200, prop_S = 0.5, duration = 10)
object <- getZstats(d, targetEvents.S = c(28, 70))
getZtests_GSD(object)
}
