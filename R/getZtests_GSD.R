#' Get group sequential design test statistics for a given trial.
#'
#' @param alpha Type I error, always one-sided.
#' @param HR.Sc.threshold Hazard ratio threshold for futility.
#' @param object An object return by \link{getZstats}.
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom stats pnorm qnorm
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 200, prop_S = 0.5, duration = 10)
#' object <- getZstats(d, targetEvents.S = c(28, 70))
#' getZtests_GSD(object)
getZtests_GSD <- function(object, alpha = 0.025, HR.Sc.threshold = 1){
  # rejection boundary
  rb <- qnorm(1-alpha)
  # object values
  z.S <- object$z.S
  z.F <- object$z.F
  p.S <- object$p.S
  p.F <- object$p.F
  hr.F.IA <- object$hr.F.IA
  hr.Sc.IA <- object$hr.Sc.IA
  obsEvents.S <- object$obsEvents.S
  obsEvents.F <- object$obsEvents.F
  nSc.IA <- object$nSc.IA
  nS.FA <- object$nS.FA
  nF.FA <- object$nF.FA

  ####### Hierarchical testing ##########
  if(hr.Sc.IA>=HR.Sc.threshold){
    Sc.futile <- TRUE
    S.selected <- TRUE
    F.selected <- FALSE
    S.reject <- (z.S[2]>rb)
    F.reject <- FALSE
    samplesize <- nSc.IA + nS.FA
  }else{
    Sc.futile <- FALSE
    S.selected <- TRUE
    F.selected <- FALSE
    S.reject <- (z.S[2]>rb)
    F.reject <- S.reject & (z.F[2]>rb)
    samplesize <- nF.FA
  }
  reject.ht <- c(S.reject = S.reject, F.reject = F.reject, S.selected = S.selected,
                 F.selected = F.selected, samplesize = samplesize, Sc.futile = Sc.futile)

  ####### Simes method ##########
  if(hr.Sc.IA>=HR.Sc.threshold){
    Sc.futile <- TRUE
    S.selected <- TRUE
    F.selected <- FALSE
    ap.SF <- simes(c(p.S[2], p.F[1])) # adjusted p-value
    z.tilde <- qnorm(1-ap.SF)
    S.reject <- (z.tilde>rb) & (z.S[2]>rb)
    F.reject <- FALSE
  }else{
    Sc.futile <- FALSE
    S.selected <- TRUE
    F.selected <- TRUE
    ap.SF <- simes(c(p.S[2], p.F[2])) # adjusted p-value
    z.tilde <- qnorm(1-ap.SF)
    S.reject <- (z.tilde>rb) & (z.S[2]>rb)
    F.reject <- (z.tilde>rb) & (z.F[2]>rb)
  }
  reject.simes <- c(S.reject = S.reject, F.reject = F.reject, S.selected = S.selected,
                    F.selected = F.selected, samplesize = samplesize, Sc.futile = Sc.futile)

  ####### Dunnett method ##########
  if(hr.Sc.IA>=HR.Sc.threshold){
    Sc.futile <- TRUE
    S.selected <- TRUE
    F.selected <- FALSE
    rho <- obsEvents.S[1]/sqrt(obsEvents.S[2]*obsEvents.F[1])
    cr <- matrix(c(1, rho, rho, 1), 2)
    ap.SF <- dunnett(c(p.S[2], p.F[1]), cr = cr) # adjusted p-value
    z.tilde <- qnorm(1-ap.SF)
    S.reject <- (z.tilde>rb) & (z.S[2]>rb)
    F.reject <- FALSE
  }else{
    Sc.futile <- FALSE
    S.selected <- TRUE
    F.selected <- TRUE
    rho <- sqrt(obsEvents.S[2]/obsEvents.F[2])
    cr <- matrix(c(1, rho, rho, 1), 2)
    ap.SF <- dunnett(c(p.S[2], p.F[2]), cr = cr) # adjusted p-value
    z.tilde <- qnorm(1-ap.SF)
    S.reject <- (z.tilde>rb) & (z.S[2]>rb)
    F.reject <- (z.tilde>rb) & (z.F[2]>rb)
  }
  reject.dunnett <- c(S.reject = S.reject, F.reject = F.reject, S.selected = S.selected,
                      F.selected = F.selected, samplesize = samplesize, Sc.futile = Sc.futile)

  # resulted test statistic
  return(list(reject.ht = reject.ht, reject.dunnett = reject.dunnett,
              reject.simes = reject.simes))
}
