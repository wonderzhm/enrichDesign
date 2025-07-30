#' Simulate an oncology trial with sub-population and full population
#'
#' This function will simulate a trial where patients will be equally randomized to treatment and control within sub-population S and its complement sub-population Sc.
#' @param n total number of subjects for full population.
#' @param prop_S proportion for sub-population group S.
#' @param ratio randomization ratio r:1, where r refers to treatment arm; for equal randomization, r=1.
#' @param duration enrollment duration in months.
#' @param hazard_S hazard rates (h_control, h_treatment) for S
#' @param hazard_Sc hazard rates (h_control, h_treatment) for Sc
#' @param dropout_S dropout hazard rates (h_control, h_treatment) for S
#' @param dropout_Sc dropout hazard rates (h_control, h_treatment) for Sc
#' @param w weight parameter in cumulative enrollment pattern.
#'
#' @return A time-to-event dataset consists of
#' \itemize{
#' \item \code{subgroup}: 1=S; 0=Sc.
#' \item \code{trt}: 1=treatment; 0=control
#' \item \code{enterTime}: subject entry time relative to the study starting date
#' \item \code{survTime}: observed event time
#' \item \code{event}: 0=censored; 1=event
#' \item \code{calendarTime}: total time on study, i.e. \code{enterTime} + \code{survTime}.
#' }
#' @importFrom rlang .data
#' @importFrom dplyr mutate bind_rows
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 100, prop_S = 0.5, ratio = 1)
#' head(d)
simu_enrich_trial <- function(n = 100, prop_S = 0.5, ratio = 1, duration = 5,
                              hazard_S = c(0.1, 0.07), hazard_Sc = c(0.11, 0.088),
                              dropout_S = c(0, 0), dropout_Sc = c(0, 0), w = 1){
  ## preparation
  nS <- round(n*prop_S)
  nSc <- n - nS
  r <- ratio/(ratio+1)
  nS_trt <- round(nS*r)
  nS_con <- nS - nS_trt
  nSc_trt <- round(nSc*r)
  nSc_con <- nSc - nSc_trt
  ## simulate for S and control
  dS_con <- simu_single_arm(n = nS_con, duration = duration, hazard = hazard_S[1],
                            dropout = dropout_S[1], w = w) %>%
    mutate(subgroup = 1, trt = 0, .before = .data$enterTime)
  ## simulate for S and trt
  dS_trt <- simu_single_arm(n = nS_trt, duration = duration, hazard = hazard_S[2],
                            dropout = dropout_S[2], w = w) %>%
    mutate(subgroup = 1, trt = 1, .before = .data$enterTime)
  ## simulate for Sc and control
  dSc_con <- simu_single_arm(n = nSc_con, duration = duration, hazard = hazard_Sc[1],
                            dropout = dropout_Sc[1], w = w) %>%
    mutate(subgroup = 0, trt = 0, .before = .data$enterTime)
  ## simulate for Sc and trt
  dSc_trt <- simu_single_arm(n = nSc_trt, duration = duration, hazard = hazard_Sc[2],
                            dropout = dropout_Sc[2], w = w) %>%
    mutate(subgroup = 0, trt = 1, .before = .data$enterTime)

  ## Combine the data
  d <- bind_rows(dS_con, dS_trt, dSc_con, dSc_trt)
  return(d)
}
