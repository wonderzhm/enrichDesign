#' Simulate a single arm survival data with uniform enrollment distribution
#'
#' @param n total number of subjects for full population.
#' @param duration enrollment duration in months; must be an integer.
#' @param hazard hazard rate
#' @param dropout dropout hazard rate
#' @param w weight parameter in cumulative enrollment pattern.
#'
#' @return A time-to-event dataset consists of
#' \itemize{
#' \item \code{enterTime}: subject entry time relative to the study starting date
#' \item \code{survTime}: observed event time
#' \item \code{event}: 0=censored; 1=event
#' \item \code{calendarTime}: total time on study, i.e. \code{enterTime} + \code{survTime}.
#' }
#' @importFrom rlang .data
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom stats rexp runif
#' @export
#'
#' @examples
#' d <- simu_single_arm(n = 20, duration = 5)
#' head(d)
simu_single_arm <- function(n = 100, duration = 5, hazard = 0.1, dropout = 0, w = 1){
  ## Simulate survival time
  surv_time <- rexp(n = n, rate = hazard)
  cen_time <- rexp(n = n, rate = max(dropout, 1e-10))

  ## Simulate enrollment
  duration <- round(duration)
  if(duration<=1){
    n_per_month <- n
  }else{
    n_per_month <- rep(NA, duration) # number of pts per month
    tmp <- 0
    for (i in 1:duration) {
      #ith month: cumulative #pts
      cN0i <- max(round((i/duration)^w * n), 1)
      n_per_month[i] <- max(cN0i - tmp, 1)
      tmp = cN0i
    }
    n_per_month[duration] = n - sum(n_per_month[1:(duration-1)])
  }
  enroll_time <- rep(NA, n)
  enroll_time[1:n_per_month[1]] = runif(n_per_month[1], min=0, max=1)
  if(duration>1){
    for (j in 2:duration){
      LL = sum(n_per_month[1:(j-1)])+1
      UU = sum(n_per_month[1:j])
      enroll_time[LL:UU] = runif(n_per_month[j], min=0, max=1) + (j - 1)
    }
  }

  ## Combine the data
  d <- data.frame(enterTime = enroll_time, surv_time = surv_time, cen_time = cen_time) %>%
    mutate(survTime = ifelse(surv_time<=cen_time, surv_time, cen_time)) %>%
    mutate(event = ifelse(surv_time<=cen_time, 1, 0)) %>%
    mutate(calendarTime = .data$enterTime + .data$survTime) %>%
    dplyr::select(-surv_time, -cen_time)
  return(d)
}
