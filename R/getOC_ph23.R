#' Get operating characteristics via simulations for the phase 2/3 design
#'
#' @param ith ith cluster used for parallel computing
#' @param seed random seed for reproducibility
#' @param nsim number of replicates
#' @param test.method method for getting adjusted p-value: \code{"dunnett"} or \code{"simes"}.
#' @param hazardC hazard rate for control arm
#' @param hazardT vector of hazard rates for treatment arms
#' @param orrC ORR for control arm
#' @param orrT vector of ORRs for control arm
#' @param rho the correlation between ORR and survival based on Gaussian copula.
#' @param dropoutC dropout hazard rate for control.
#' @param dropoutT vector of dropout hazard rates for treatment arms.
#' @param accrual_rate_stage1 number of subjects enrolled per month for stage 1.
#' @param accrual_rate_stage2 number of subjects enrolled per month for stage 2.
#' @param n1_per_arm total number of subjects per arm in stage 1.
#' @param n2_per_arm total number of subjects per arm in stage 2.
#' @param targetEventsIA_all target number of events at IA for all subjects from both stage 1 and 2.
#' @param targetEventsFA_all target number of events at FA for all subjects from both stage 1 and 2.
#' @param wei_IA weights (w1, w2) used in combination test at IA.
#' @param wei_FA weights (w1, w2) used in combination test at IA.
#' @param bound_z_IA z-scale rejection boundary at IA.
#' @param bound_z_FA z-scale rejection boundary at FA.
#' @param alpha Type I error, always one-sided.
#' @param update_bound_FA whether to re-calculate the FA rejection boundary using the updated
#' correlation between IA and FA test statistics.
#'
#' @return It returns a matrix with each row corresponding to the analysis results for each trial,
#' where the analysis results include: "Selected dose", "IA time", "IA rejct", "FA reject", "FA time",
#' "Study duration", "Total sample size", "observed correlation", and "actual FA boundary".
#' @importFrom dplyr filter group_by summarise
#' @importFrom rlang .data
#' @importFrom gsDesign gsDesign sfLDOF
#' @export
#'
#' @examples
#' res <- getOC_ph23(seed = 24232, nsim=10)
#' apply(res, 2, mean, rm=TRUE)
getOC_ph23 <- function(ith = 1, seed = 2024, nsim = 1000, test.method = "dunnett",
                       hazardC = 0.1, hazardT = c(0.1, 0.1, 0.1),
                       orrC = 0.2, orrT = c(0.2, 0.2, 0.2), rho = 0.7,
                       dropoutC = 0, dropoutT = c(0, 0, 0),
                       accrual_rate_stage1 = 25, accrual_rate_stage2 = 25,
                       n1_per_arm = 50, n2_per_arm = 150,
                       targetEventsIA_all = 240, targetEventsFA_all = 320,
                       wei_IA = c(0.8, 0.6), wei_FA = c(0.8, 0.6),
                       bound_z_IA = 2.44, bound_z_FA = 2,
                       alpha = 0.025, update_bound_FA = TRUE){
  ## check weights
  if((abs(sum(wei_IA^2)-1)+abs(sum(wei_FA^2)-1))>0.000001) stop("w1^2 + w2^2 should be equal to one!")

  ## Start simulation
  set.seed(seed*ith)
  results.names <- c("Selected dose", "IA time", "IA rejct", "FA reject", "FA time",
                     "Study duration", "Total sample size", "observed correlation",
                     "actual FA boundary")
  results <- matrix(NA, nrow=nsim, ncol=length(results.names))
  colnames(results) <- results.names

  for(sim in 1:nsim){
    ## Simulate survival times and tumor responses
    num_trt <- length(hazardT)
    dat <- sim_ph23(hazardC = hazardC, hazardT = hazardT,
                    orrC = orrC, orrT = orrT, rho = rho,
                    accrual_rate_stage1 = accrual_rate_stage1,
                    accrual_rate_stage2 = accrual_rate_stage2,
                    n1_per_arm = n1_per_arm, n2_per_arm = n2_per_arm,
                    dropoutC = dropoutC, dropoutT = dropoutT)

    ## dose selection
    d1 <- dat %>% filter(.data$stage==1)
    orr_est <- d1 %>% filter(.data$trt!=0) %>%
      group_by(.data$trt) %>% summarise(orr = mean(.data$response))
    selected <- orr_est$trt[which.max(orr_est$orr)]

    ## IA
    res_IA <- getZstat(dat = dat, w = wei_IA, selected = selected, targetEvents = targetEventsIA_all,
                       test.method = test.method)
    IA_time <- res_IA$cut_time
    obsEventsIA_1 <- res_IA$obsEvents_stage1
    IA_sample_size <- res_IA$sample_size
    z_IA <- res_IA$z
    IA_reject <- (z_IA >= bound_z_IA)

    ## FA
    if(IA_reject){
      FA_reject <- IA_reject
      FA_time <- NA
      Study_duration <- IA_time
      total_sample_size <-  IA_sample_size
      obs_cor <- NA
      if(update_bound_FA) bound_z_FA_obs <- NA
    }else{
      res_FA <- getZstat(dat = dat, w = wei_FA, selected = selected, targetEvents = targetEventsFA_all,
                         test.method = test.method)
      FA_time <- res_FA$cut_time
      obsEventsFA_1 <- res_FA$obsEvents_stage1
      FA_sample_size <- res_FA$sample_size
      z_FA <- res_FA$z
      obs_cor <- wei_IA[1]*wei_FA[1]*sqrt(obsEventsIA_1/obsEventsFA_1) +
        wei_IA[2]*wei_FA[2]*sqrt((targetEventsIA_all-obsEventsIA_1)/(targetEventsFA_all-obsEventsFA_1))
      if(update_bound_FA){
        obs_info_fraction <- obs_cor^2
        du <- gsDesign(k = 2, test.type = 1, alpha = alpha, sfu = sfLDOF,
                       n.I = c(targetEventsIA_all, targetEventsIA_all/obs_info_fraction),
                       maxn.IPlan = targetEventsFA_all)
        bound_z_FA_obs <- du$upper$bound[2]
        FA_reject <- (z_FA >= bound_z_FA_obs)
      }else{
        FA_reject <- (z_FA >= bound_z_FA)
      }
      Study_duration <- FA_time
      total_sample_size <-  FA_sample_size
    }

    ## Save results
    results[sim, 1] <- selected
    results[sim, 2] <- IA_time
    results[sim, 3] <- IA_reject
    results[sim, 4] <- FA_reject
    results[sim, 5] <- FA_time
    results[sim, 6] <- Study_duration
    results[sim, 7] <- total_sample_size
    results[sim, 8] <- obs_cor
    if(update_bound_FA){
      results[sim, 9] <- bound_z_FA_obs
    }else{
      results[sim, 9] <- bound_z_FA
    }
  }
  return(results)
}
