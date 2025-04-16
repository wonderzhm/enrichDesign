#' Get hazard ratio thresholds for adaptive enrichment designs
#'
#' @param n total number of subjects for full population.
#' @param prop_S proportion for sub-population group S.
#' @param duration enrollment duration in months.
#' @param targetEvents.S The target number of events in S; currently only support length of 2, where
#' futility or population selection will be given at IA without efficacy testing.
#' @param CP.Sc.threshold conditional power threshold for futility.
#' @param CP.S.threshold conditional power threshold for S.
#' @param CP.F.threshold conditional power threshold for F.
#' @param hazard_S hazard rates (h_control, h_treatment) for S
#' @param hazard_Sc hazard rates (h_control, h_treatment) for Sc
#' @param dropout_S dropout hazard rates (h_control, h_treatment) for S
#' @param dropout_Sc dropout hazard rates (h_control, h_treatment) for Sc
#' @param ratio randomization ratio r:1, where r refers to treatment arm; for equal randomization, r=1.
#' @param alpha Type I error, always one-sided.
#'
#' @return It returns a list.
#' @importFrom rpact getDesignGroupSequential getPowerSurvival getEventProbabilities
#' @export
#'
#' @examples
#' getDesignParameters()
getDesignParameters <- function(
    n = 800, prop_S = 0.5, duration = 20, targetEvents.S = c(112, 280),
    CP.Sc.threshold = 0.3, CP.S.threshold = 0.8, CP.F.threshold = 0.8,
    hazard_S = c(0.1, 0.07), hazard_Sc = c(0.11, 0.088),
    dropout_S = c(0, 0), dropout_Sc = c(0, 0), ratio = 1, alpha = 0.025){
  # prepare variables
  nF <- n
  nS <- n*prop_S
  nSc <- nF - nS
  hrS <- hazard_S[2]/hazard_S[1]
  hrSc <- hazard_Sc[2]/hazard_Sc[1]
  hrF <- exp(prop_S*log(hrS)+(1-prop_S)*log(hrSc))
  medS_con <- log(2)/hazard_S[1]
  medSc_con <- log(2)/hazard_Sc[1]
  medF_con <- medS_con*prop_S+medSc_con*(1-prop_S)
  dropoutTime <- 12
  dropoutRate1S = 1-exp(-dropout_S[2]*dropoutTime)
  dropoutRate2S = 1-exp(-dropout_S[1]*dropoutTime)
  dropoutRate1Sc = 1-exp(-dropout_Sc[2]*dropoutTime)
  dropoutRate2Sc = 1-exp(-dropout_Sc[1]*dropoutTime)
  dropoutRate1F = dropoutRate1S*prop_S + dropoutRate1Sc*(1-prop_S)
  dropoutRate2F = dropoutRate2S*prop_S + dropoutRate2Sc*(1-prop_S)
  targetEvents_S <- targetEvents.S
  ## GSD design as a reference
  d <- getDesignGroupSequential(kMax = 2,
                                alpha = alpha,
                                sided = 1,
                                informationRates = targetEvents_S/targetEvents_S[2],
                                typeOfDesign = "noEarlyEfficacy")
  # power for S
  dS <- getPowerSurvival(design = d, maxNumberOfEvents = targetEvents_S[2],
                         lambda1 = hazard_S[2],
                         lambda2 = hazard_S[1],
                         allocationRatioPlanned = ratio,
                         accrualTime = c(0, duration),
                         accrualIntensity = nS/duration,
                         directionUpper = FALSE,
                         dropoutRate1 = dropoutRate1S,
                         dropoutRate2 = dropoutRate2S,
                         dropoutTime = dropoutTime)
  dSpower <- dS$overallReject
  times <- dS$analysisTime
  # expected number of events for S^c
  probEvent <- getEventProbabilities(time = times, maxNumberOfSubjects = nSc,
                                     lambda1 = hazard_Sc[2],
                                     lambda2 = hazard_Sc[1],
                                     allocationRatioPlanned = ratio,
                                     accrualTime = c(0, duration),
                                     accrualIntensity = nSc/duration,
                                     dropoutRate1 = dropoutRate1Sc,
                                     dropoutRate2 = dropoutRate2Sc,
                                     dropoutTime = dropoutTime)
  targetEvents_Sc <- ceiling(probEvent$overallEventProbabilities*nSc)
  # expected number of events for F
  targetEvents_F <- (targetEvents_Sc+targetEvents_S)
  # power for F
  dF <- getPowerSurvival(design = d, hazardRatio = hrF,
                         maxNumberOfEvents = targetEvents_F[2],
                         lambda2 = log(2) / medF_con,
                         allocationRatioPlanned = ratio,
                         accrualTime = c(0, duration),
                         accrualIntensity = nF/duration,
                         dropoutRate1 = dropoutRate1F,
                         dropoutRate2 = dropoutRate2F,
                         dropoutTime = dropoutTime,
                         directionUpper = FALSE)
  dFpower <- dF$overallReject
  # find HR threshold for S
  HRs <- seq(0.5, 1.5, 0.01)
  CPs <- rep(NA, length(HRs))
  for(i in 1:length(HRs)){
    D <- targetEvents_S # number of events at IA and FA
    delta <- -log(hrS) # assumed delta: negative log hazard ratio
    hr.est <- HRs[i]
    z <- -log(hr.est)/sqrt(4/D[1])
    b <- -log(dS$criticalValuesEffectScale[2])/sqrt(4/D[2]) # boundary at FA for S
    CPs[i] <- 1-pnorm( b*sqrt(D[2]/(D[2]-D[1])) - z*sqrt(D[1]/(D[2]-D[1])) - delta/sqrt(4/(D[2]-D[1])) )
  }
  HR.S.threshold <- max(HRs[CPs>=CP.S.threshold])
  # find HR threshold for Sc
  HRs <- seq(0.5, 1.5, 0.01)
  CPs <- rep(NA, length(HRs))
  for(i in 1:length(HRs)){
    D <- targetEvents_F # number of events at IA and FA
    delta <- -log(hrF) # assumed delta: negative log hazard ratio
    hr.est <- HRs[i]
    z <- -log(hr.est)/sqrt(4/D[1])
    b <- -log(dF$criticalValuesEffectScale[2])/sqrt(4/D[2]) # boundary at FA for S
    CPs[i] <- 1-pnorm( b*sqrt(D[2]/(D[2]-D[1])) - z*sqrt(D[1]/(D[2]-D[1])) - delta/sqrt(4/(D[2]-D[1])) )
  }
  HR.F.threshold <- max(HRs[CPs>=CP.F.threshold])
  # find HR threshold for Sc
  HRs <- seq(0.5, 1.5, 0.01)
  CPs <- rep(NA, length(HRs))
  for(i in 1:length(HRs)){
    D <- targetEvents_Sc # number of events at IA and FA
    delta <- -log(hrSc) # assumed delta: negative log hazard ratio
    hr.est <- HRs[i]
    z <- -log(hr.est)/sqrt(4/D[1])
    b <- -log(dF$criticalValuesEffectScale[2])/sqrt(4/D[2]) # boundary at FA for S
    CPs[i] <- 1-pnorm( b*sqrt(D[2]/(D[2]-D[1])) - z*sqrt(D[1]/(D[2]-D[1])) - delta/sqrt(4/(D[2]-D[1])) )
  }
  HR.Sc.threshold <- min(HRs[CPs<=CP.Sc.threshold])

  return(list(HR.S.threshold = HR.S.threshold,
              HR.F.threshold = HR.F.threshold,
              HR.Sc.threshold = HR.Sc.threshold,
              GSD.reject.S = dSpower,
              GSD.reject.F = dFpower,
              dS = dS,
              dF = dF))
}
