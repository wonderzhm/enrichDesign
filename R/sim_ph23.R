#' Simulate a seamless phase 2/3 trial
#'
#' This function will simulate a phase 2/3 trial where patients will be randomized to all arms in stage 1
#' and only the selected arm and control arm will go into stage 2. Since we don't know the arm selection in
#' advance, we will simulated patients for all arms in stage 2 but with adjusted accrual rate assuming stage 2
#' only contains the selected arm and the control.
#' @param hazardC hazard rate for control arm
#' @param hazardT vector of hazard rates for treatment arms
#' @param orrC ORR for control arm
#' @param orrT vector of ORRs for control arm
#' @param rho the correlation between ORR and survival based on Gaussian copula.
#' @param accrual_rate_stage1 number of subjects enrolled per month for stage 1.
#' @param accrual_rate_stage2 number of subjects enrolled per month for stage 2.
#' @param n1_per_arm total number of subjects per arm in stage 1.
#' @param n2_per_arm total number of subjects per arm in stage 2.
#' @param dropoutC dropout hazard rate for control.
#' @param dropoutT vector of dropout hazard rates for treatment arms.
#'
#' @return A time-to-event dataset consists of
#' \itemize{
#' \item \code{stage}: 1=phase II; 2=phase III.
#' \item \code{trt}: 0=control; 1=dose 1; 2=dose 2; ...
#' \item \code{enterTime}: subject entry time relative to the study starting date
#' \item \code{response}: 0=no response; 1=response
#' \item \code{survTime}: observed event time
#' \item \code{event}: 0=censored; 1=event
#' \item \code{calendarTime}: total time on study, i.e. \code{enterTime} + \code{survTime}.
#' }
#' @importFrom rlang .data
#' @importFrom dplyr mutate
#' @importFrom stats pnorm qnorm rexp
#' @importFrom simtrial rpwexp_enroll
#' @export
#'
#' @examples
#' d <- sim_ph23(n1_per_arm = 20, n2_per_arm = 30)
#' head(d)
sim_ph23 <- function(hazardC = 0.1, hazardT = c(0.09, 0.08, 0.07),
                     orrC = 0.2, orrT = c(0.25, 0.30, 0.40), rho = 0.7,
                     accrual_rate_stage1 = 25, accrual_rate_stage2 = 25,
                     n1_per_arm = 50, n2_per_arm = 150,
                     dropoutC = 0, dropoutT = c(0, 0, 0)){
  ## Simulate survival times and tumor responses
  num_trt <- length(hazardT)
  n_per_arm <- n1_per_arm + n2_per_arm
  trt <- rep(0:num_trt, each = n_per_arm)
  ntotal <- length(trt)
  z <- MASS::mvrnorm(ntotal, mu = c(0,0), Sigma = matrix(c(1, rho, rho, 1), 2, 2))
  z_surv <- z[,1]
  z_orr <- z[,2]
  exp_rates <- rep(c(hazardC, hazardT), each = n_per_arm)
  orrs <- rep(c(orrC, orrT), each = n_per_arm)
  surv_time <- -pnorm(q = z_surv, log.p = TRUE)/exp_rates
  response <- (z_orr <= qnorm(orrs))+0
  drop_rates <- rep(c(dropoutC, dropoutT), each = n_per_arm)
  drop_rates[drop_rates<=0] <- 1e-10
  cen_time <- rexp(ntotal, rate = drop_rates)

  ## Simulate enrollment
  ## Stage 1: 50 patients per arm at rate 25 patients per month
  enroll_time <- rep(NA, ntotal)
  stage <- rep( c(rep(1, n1_per_arm), rep(2, n2_per_arm)), num_trt+1)
  stage1.indx <- which(stage==1)
  stage1.size <- length(stage1.indx)
  stage1.enroll <- rpwexp_enroll(
    n = stage1.size, enroll_rate = data.frame(
      duration = c(n1_per_arm*(num_trt+1)/accrual_rate_stage1),
      rate = accrual_rate_stage1))
  enroll_time[stage1.indx] <- sample(stage1.enroll)
  stage2.indx <- which(stage==2)
  stage2.size <- length(stage2.indx)
  stage2.enroll <- rpwexp_enroll(
    n = stage2.size, enroll_rate = data.frame(
      duration = c(n2_per_arm*2/accrual_rate_stage2),
      rate = accrual_rate_stage2/2*(num_trt+1))) + max(stage1.enroll)
  enroll_time[stage2.indx] <- sample(stage2.enroll)

  ## Combine the data
  dat <- data.frame(stage = stage, trt = trt, response = response, enterTime = enroll_time,
                    surv_time = surv_time, cen_time = cen_time) %>%
    mutate(survTime = ifelse(surv_time<=cen_time, surv_time, cen_time)) %>%
    mutate(calendarTime = .data$enterTime + .data$survTime) %>%
    mutate(event = ifelse(surv_time<=cen_time, 1, 0)) %>%
    dplyr::select(-surv_time, -cen_time)
  return(dat)
}
