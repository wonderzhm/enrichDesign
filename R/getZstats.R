#' Get test statistics for a given trial.
#'
#' @param dat A time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param targetEvents.S The target number of events in S; currently only support length of 2, where
#' futility or population selection will be given at IA without efficacy testing.
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 200, prop_S = 0.5, duration = 10)
#' getZstats(d, targetEvents.S = c(28, 70))
getZstats <- function(dat, targetEvents.S){
  targetEvents.S <- ceiling(targetEvents.S)
  ## three datasets
  dS <- dat %>% filter(.data$subgroup==1)
  dSc <- dat %>% filter(.data$subgroup==0)
  dF <- dat
  ## to save
  z.S <- z.F <- rep(NA, 2) # z-statistics
  p.S <- p.F <- rep(NA, 2) # p-values
  obsEvents.S <- obsEvents.F <- rep(NA, 2) # observed number of events
  ## cut data at IA
  # population S
  d <- cut_by_event(dS, targetEvents = targetEvents.S[1])
  IA_time <- d$calendarCutoff[1]
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = NULL)
  z.S[1] <- res$z # non-adjusted Z statistic
  p.S[1] <- res$p
  obsEvents.S[1] <- sum(res$obs)
  hr.S.IA <- exp(-(res$z)/sqrt(1/sum(1/res$obs)))
  # population F
  d <- cut_by_date(dF, cut_time = IA_time)
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = d$subgroup)
  z.F[1] <- res$z # non-adjusted Z statistic
  p.F[1] <- res$p
  obsEvents.F[1] <- sum(res$obs)
  hr.F.IA <- exp(-(res$z)/sqrt(1/sum(1/res$obs)))
  # population Sc
  d <- cut_by_date(dSc, cut_time = IA_time)
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = NULL)
  hr.Sc.IA <- exp(-(res$z)/sqrt(1/sum(1/res$obs)))
  nSc.IA <- nrow(d)

  ## cut data at FA
  # population S
  d <- cut_by_event(dS, targetEvents = targetEvents.S[2])
  FA_time <- d$calendarCutoff[1]
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = NULL)
  z.S[2] <- res$z # non-adjusted Z statistic
  p.S[2] <- res$p
  obsEvents.S[2] <- sum(res$obs)
  nS <- nrow(d)
  # population F
  d <- cut_by_date(dF, cut_time = FA_time)
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = d$subgroup)
  z.F[2] <- res$z # non-adjusted Z statistic
  p.F[2] <- res$p
  obsEvents.F[2] <- sum(res$obs)
  nF <- nrow(d)

  # resulted test statistic
  return(list(z.S = z.S, z.F = z.F, p.S = p.S, p.F = p.F,
              hr.Sc.IA = hr.Sc.IA, hr.S.IA = hr.S.IA, hr.F.IA = hr.F.IA,
              obsEvents.S = obsEvents.S, obsEvents.F = obsEvents.F,
              nSc.IA = nSc.IA, nS.FA = nS, nF.FA = nF))
}
