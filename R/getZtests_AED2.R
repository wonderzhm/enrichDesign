#' Get adaptive enrichment design test statistics for a given trial.
#'
#' This function always pick one population moving forward.
#' @param alpha Type I error, always one-sided.
#' @param HR.Sc.threshold Hazard ratio threshold for futility of Sc.
#' @param object An object return by \link{getZstats}.
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom stats pnorm qnorm
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 200, prop_S = 0.5, duration = 10)
#' object <- getZstats(d, targetEvents.S = c(28, 70))
#' getZtests_AED2(object)
getZtests_AED2 <- function(object, alpha = 0.025, HR.Sc.threshold = 1){
  # rejection boundary
  rb <- qnorm(1-alpha)
  # object values
  z.S <- object$z.S
  z.F <- object$z.F
  p.S <- object$p.S
  p.F <- object$p.F
  hr.S.IA <- object$hr.S.IA
  hr.F.IA <- object$hr.F.IA
  hr.Sc.IA <- object$hr.Sc.IA
  obsEvents.S <- object$obsEvents.S
  obsEvents.F <- object$obsEvents.F
  nSc.IA <- object$nSc.IA
  nS.FA <- object$nS.FA
  nF.FA <- object$nF.FA

  # independent incremental statistics and p-values
  zii.S <- (sqrt(obsEvents.S[2])*z.S[2] - sqrt(obsEvents.S[1])*z.S[1])/
    (sqrt(obsEvents.S[2]-obsEvents.S[1]))
  pii.S <- 1- pnorm(zii.S)
  zii.F <- (sqrt(obsEvents.F[2])*z.F[2] - sqrt(obsEvents.F[1])*z.F[1])/
    (sqrt(obsEvents.F[2]-obsEvents.F[1]))
  pii.F <- 1- pnorm(zii.F)

  # temp function to do adjustment
  adjustment2 <- function(ap1.SF){
    # Population selection at IA
    if(hr.Sc.IA>=HR.Sc.threshold){
      Sc.futile <- TRUE
      S.selected <- TRUE
      F.selected <- FALSE
      w <- sqrt(obsEvents.S[1]/obsEvents.S[2])
      z.tilde <- w*qnorm(1-ap1.SF) + sqrt(1-w^2)*zii.S
      S.reject <- (z.tilde>rb) & (z.S[2]>rb)
      F.reject <- FALSE
      samplesize <- nSc.IA + nS.FA
    }else if(hr.F.IA<=hr.S.IA){
      ## Using (obsEvents.F[1]*log(hr.F.IA))<=(obsEvents.S[1]*log(hr.S.IA)) does not work well
      Sc.futile <- FALSE
      S.selected <- FALSE
      F.selected <- TRUE
      w <- sqrt(obsEvents.F[1]/obsEvents.F[2])
      z.tilde <- w*qnorm(1-ap1.SF) + sqrt(1-w^2)*zii.F
      F.reject <- (z.tilde>rb) & (z.F[2]>rb)
      S.reject <- F.reject & (z.S[2]>rb)
      samplesize <- nF.FA
    }else{
      Sc.futile <- FALSE
      S.selected <- TRUE
      F.selected <- FALSE
      w <- sqrt(obsEvents.S[1]/obsEvents.S[2])
      z.tilde <- w*qnorm(1-ap1.SF) + sqrt(1-w^2)*zii.S
      S.reject <- (z.tilde>rb) & (z.S[2]>rb)
      F.reject <- S.reject & (z.F[2]>rb)
      samplesize <- nF.FA
    }
    res <- c(S.reject = S.reject, F.reject = F.reject, S.selected = S.selected,
             F.selected = F.selected, samplesize = samplesize, Sc.futile = Sc.futile)
    return(res)
  }

  ####### Simes method ##########
  # adjusted p-values
  ap1.SF <- simes(c(p.S[1], p.F[1]))
  reject.simes <- adjustment2(ap1.SF)

  ####### Dunnett method ##########
  # adjusted p-values
  rho <- sqrt(obsEvents.S[1]/obsEvents.F[1])
  cr <- matrix(c(1, rho, rho, 1), 2)
  ap1.SF <- dunnett(c(p.S[1], p.F[1]), cr = cr)
  reject.dunnett <- adjustment2(ap1.SF)

  # resulted test statistic
  return(list(reject.dunnett = reject.dunnett, reject.simes = reject.simes))
}
