#' Get operating characteristics via simulations for an enrichment design
#'
#' @param seed random seed for reproducibility
#' @param nsim number of replicates
#' @param n total number of subjects for full population.
#' @param prop_S proportion for sub-population group S.
#' @param duration enrollment duration in months.
#' @param targetEvents.S The target number of events in S; currently only support length of 2, where
#' futility or population selection will be given at IA without efficacy testing.
#' @param HR.Sc.threshold Hazard ratio threshold for futility.
#' @param HR.S.threshold Hazard ratio threshold for S.
#' @param HR.F.threshold Hazard ratio threshold for F.
#' @param hazard_S hazard rates (h_control, h_treatment) for S
#' @param hazard_Sc hazard rates (h_control, h_treatment) for Sc
#' @param dropout_S dropout hazard rates (h_control, h_treatment) for S
#' @param dropout_Sc dropout hazard rates (h_control, h_treatment) for Sc
#' @param w weight parameter in cumulative enrollment pattern.
#' @param ratio randomization ratio r:1, where r refers to treatment arm; for equal randomization, r=1.
#' @param alpha Type I error, always one-sided.
#'
#' @return It returns a list.
#' @export
#'
#' @examples
#' res <- getOC(seed = 24232, nsim=10)
#' lapply(res, function(x) mean(apply(x[,1:2], 1, any, na.rm=TRUE)))
getOC <- function(seed = 2024, nsim = 1000, n = 800, prop_S = 0.5, duration = 25,
                  targetEvents.S = c(112, 280), HR.Sc.threshold = 1, HR.S.threshold = 0.8,
                  HR.F.threshold = 0.8, hazard_S = c(0.1, 0.1), hazard_Sc = c(0.12, 0.12),
                  dropout_S = c(0, 0), dropout_Sc = c(0, 0), w = 1, ratio = 1, alpha = 0.025){
  ## Start simulation
  set.seed(seed)
  AED.simes <- matrix(NA, nrow = nsim, ncol = 6)
  colnames(AED.simes) <- c("S.reject", "F.reject", "S.selected", "F.selected", "samplesize", "Sc.futile")
  AED.dunnett <- AED2.dunnett <- AED2.simes <- AED.simes
  GSD.ht <- GSD.simes <- GSD.dunnett <- AED.simes

  for(sim in 1:nsim){
    ## Simulate survival times and tumor responses
    d <- simu_enrich_trial(n = n, prop_S = prop_S, ratio = ratio, duration = duration,
                           hazard_S = hazard_S, hazard_Sc = hazard_Sc,
                           dropout_S = dropout_S, dropout_Sc = dropout_Sc, w = w)
    ## get Z statistics
    Zstats <- getZstats(dat = d, targetEvents.S = targetEvents.S)

    ## get test results: AED
    res <- getZtests_AED(object = Zstats, alpha = alpha, HR.Sc.threshold = HR.Sc.threshold,
                         HR.S.threshold = HR.S.threshold, HR.F.threshold = HR.F.threshold)
    AED.dunnett[sim, ] <- res$reject.dunnett
    AED.simes[sim, ] <- res$reject.simes

    ## get test results: AED and pick the winner
    res <- getZtests_AED2(object = Zstats, alpha = alpha, HR.Sc.threshold = HR.Sc.threshold)
    AED2.dunnett[sim, ] <- res$reject.dunnett
    AED2.simes[sim, ] <- res$reject.simes

    ## get test results: GSD
    res <- getZtests_GSD(object = Zstats, alpha = alpha, HR.Sc.threshold = HR.Sc.threshold)
    GSD.ht[sim,] <- res$reject.ht
    GSD.dunnett[sim, ] <- res$reject.dunnett
    GSD.simes[sim, ] <- res$reject.simes
  }
  return(list(AED.dunnett = AED.dunnett, AED.simes = AED.simes,
              GSD.ht = GSD.ht, GSD.dunnett = GSD.dunnett, GSD.simes = GSD.simes,
              AED2.dunnett = AED2.dunnett, AED2.simes = AED2.simes))
}
