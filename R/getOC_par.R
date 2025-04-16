#' Get operating characteristics via simulations for an enrichment design
#'
#' @param ncore number of cores for parallel computing.
#' @param seed random seed for reproducibility; each core will receive a seed in \code{(1:ncore)*seed}
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
#' @importFrom parallel makeCluster detectCores clusterCall clusterExport parLapply stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
#'
#' @examples
#' res <- getOC_par(ncore = 2, seed = 24232, nsim=10)
#' lapply(res, function(x) mean(apply(x[,1:2], 1, any, na.rm=TRUE)))
getOC_par <- function(ncore = 2, seed = 2024, nsim = 1000, n = 800, prop_S = 0.5, duration = 25,
                      targetEvents.S = c(112, 280), HR.Sc.threshold = 1, HR.S.threshold = 0.8,
                      HR.F.threshold = 0.8, hazard_S = c(0.1, 0.1), hazard_Sc = c(0.12, 0.12),
                      dropout_S = c(0, 0), dropout_Sc = c(0, 0), w = 1, ratio = 1, alpha = 0.025){
  ## Start simulation
  if(ncore > detectCores()) stop("ncore is too big!")
  if(ncore <=1) stop("ncore needs to be at least 2!")
  cl <- makeCluster(ncore, type = "PSOCK")
  registerDoParallel(cl)
  ## Export Functions to the Cluster
  tmp1 <- clusterCall(cl, function() {library(enrichDesign)})
  ## Function input
  seeds <- (1:ncore)*seed
  results <- parLapply(
    cl, seeds, fun = getOC, nsim = nsim, n = n, prop_S = prop_S, duration = duration,
    targetEvents.S = targetEvents.S, HR.Sc.threshold = HR.Sc.threshold, HR.S.threshold = HR.S.threshold,
    HR.F.threshold = HR.F.threshold, hazard_S = hazard_S, hazard_Sc = hazard_Sc, dropout_S = dropout_S,
    dropout_Sc = dropout_Sc, w = w, ratio = ratio, alpha = alpha)
  res <- results[[1]]
  for(i in 2:ncore){
    res$AED.dunnett <- rbind(res$AED.dunnett, results[[i]]$AED.dunnett)
    res$AED.simes <- rbind(res$AED.simes, results[[i]]$AED.simes)
    res$GSD.ht <- rbind(res$GSD.ht, results[[i]]$GSD.ht)
    res$GSD.dunnett <- rbind(res$GSD.dunnett, results[[i]]$GSD.dunnett)
    res$GSD.simes <- rbind(res$GSD.simes, results[[i]]$GSD.simes)
    res$AED2.dunnett <- rbind(res$AED2.dunnett, results[[i]]$AED2.dunnett)
    res$AED2.simes <- rbind(res$AED2.simes, results[[i]]$AED2.simes)
  }
  ## Close Cluster
  stopCluster(cl)
  return(res)
}
