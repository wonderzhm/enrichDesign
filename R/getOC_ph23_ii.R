#' Get operating characteristics via simulations for the phase 2/3 design
#'
#' This function uses the independent incremental approach
#' instead of the disjoint-subject approach used in \link{getOC_ph23}.
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
#' @param w weight parameter for stage 1 data.
#' @param targetEventsIA_all target number of events at IA for all subjects from both stage 1 and 2.
#' @param targetEventsFA_all target number of events at FA for all subjects from both stage 1 and 2.
#' @param bound_z_IA z-scale rejection boundary at IA.
#' @param bound_z_FA z-scale rejection boundary at FA.
#' @param alpha Type I error, always one-sided.
#' @param update_bound_FA whether to re-calculate the FA rejection boundary using the updated
#' correlation between IA and FA test statistics.
#' @param nonselected_max_followup The final data cutoff time for non-selected arms.
#' The default \code{NULL} means the non-selected arms will not be followed up after dose selection.
#' @param dose_selection_endpoint either "ORR" or "Survival"
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
#' res <- getOC_ph23_ii(seed = 24232, nsim=10)
#' apply(res, 2, mean, rm=TRUE)
getOC_ph23_ii <- function(ith = 1, seed = 2024, nsim = 1000, test.method = "dunnett",
                          hazardC = 0.1, hazardT = c(0.1, 0.1, 0.1),
                          orrC = 0.2, orrT = c(0.2, 0.2, 0.2), rho = 0.7,
                          dropoutC = 0, dropoutT = c(0, 0, 0),
                          accrual_rate_stage1 = 25, accrual_rate_stage2 = 25,
                          n1_per_arm = 50, n2_per_arm = 150,
                          w = NULL, targetEventsIA_all = 240, targetEventsFA_all = 320,
                          bound_z_IA = 2.44, bound_z_FA = 2,
                          alpha = 0.025, update_bound_FA = TRUE,
                          nonselected_max_followup = NULL,
                          dose_selection_endpoint = "ORR"){

  ## Start simulation
  set.seed(seed*ith)
  results.names <- c("Selected dose", "IA time", "IA reject", "FA reject", "FA time",
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
    if(dose_selection_endpoint == "ORR"){
      orr_est <- d1 %>% filter(.data$trt!=0) %>%
        group_by(.data$trt) %>% summarise(orr = mean(.data$response))
      selected <- orr_est$trt[which.max(orr_est$orr)]
    }else{
      IAd_time <- max(d1$enterTime) + 1e-10
      d1IAd <- cut_by_date(d1, cut_time = IAd_time)
      pvalues <- rep(NA, num_trt)
      for(i in 1:num_trt){
        d1IAi <- d1IAd %>% filter(.data$trt%in%c(0, i))
        res <- nph::logrank.test(time = d1IAi$survTimeCut, event = d1IAi$eventCut,
                                 group = as.factor(d1IAi$trt), alternative = c("greater"),
                                 rho = 0, gamma = 0, event_time_weights = NULL)
        pvalues[i] <- res$test$p
      }
      selected <- which.min(pvalues)
    }

    ## non-selected arms will be cut at dose selection time
    dat2 <- cut_nonselected_arms(dat = dat, selected = selected,
                                 nonselected_max_followup = nonselected_max_followup)

    ## IA
    res_IA <- getZ1(dat = dat2, w = w, selected = selected, targetEvents = targetEventsIA_all,
                    test.method = test.method)
    IA_time <- res_IA$cut_time
    obsEventsIA <- res_IA$obsEventsIA1
    IA_sample_size <- res_IA$sample_size
    z_IA_tilde <- res_IA$Z1_tilde
    z_IA <- res_IA$Z1
    IA_reject <- (z_IA_tilde >= bound_z_IA)

    ## FA
    if(IA_reject){
      FA_reject <- IA_reject
      FA_time <- NA
      Study_duration <- IA_time
      total_sample_size <-  IA_sample_size
      obs_cor <- NA
      if(update_bound_FA) bound_z_FA_obs <- NA
    }else{
      d <- dat2 %>% filter(.data$trt%in%c(0, selected))
      dFA <- cut_by_event(d, targetEvents = targetEventsFA_all)
      FA_time <- dFA$calendarCutoff[1]
      obsEventsFA <- sum(dFA$eventCut)
      res <- nph::logrank.test(time =dFA$survTimeCut, event = dFA$eventCut,
                               group = as.factor(dFA$trt), alternative = c("greater"),
                               rho = 0, gamma = 0, event_time_weights = NULL)
      z_FA <- res$test$z
      z_FA_tilde <- z_FA + sqrt(obsEventsIA/obsEventsFA)*(z_IA_tilde-z_IA)
      FA_sample_size <- nrow(dFA) + (num_trt-1)*n1_per_arm
      obs_cor <- sqrt(obsEventsIA/obsEventsFA)
      if(update_bound_FA){
        obs_info_fraction <- obs_cor^2
        du <- gsDesign(k = 2, test.type = 1, alpha = alpha, sfu = sfLDOF,
                       n.I = c(targetEventsIA_all, targetEventsIA_all/obs_info_fraction),
                       maxn.IPlan = targetEventsFA_all)
        bound_z_FA_obs <- du$upper$bound[2]
        FA_reject <- (z_FA_tilde >= bound_z_FA_obs)
      }else{
        FA_reject <- (z_FA_tilde >= bound_z_FA)
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
