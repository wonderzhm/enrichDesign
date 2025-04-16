#' Cut a dataset for analysis at a specified date
#'
#' @param data A time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param cut_time Date relative to start of randomization (\code{calendarTime} from
#' input dataset) at which dataset is to be cut off for analysis.

#'
#' @return A data frame ready for survival analysis including columns \code{calendarCutoff},
#' \code{survTimeCut}, \code{eventCut}.
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 100, prop_S = 0.5, ratio = 1)
#' dcut <- cut_by_date(d, 10)
cut_by_date <- function (data, cut_time) {
  data0 = data
  data0$calendarCutoff = cut_time
  data0$survTimeCut = ifelse(data0$calendarTime <= data0$calendarCutoff,
                             data0$survTime, data0$calendarCutoff - data0$enterTime)
  data0$eventCut = ifelse(data0$calendarTime <= data0$calendarCutoff,
                          data0$event, 0)
  return(data0[data0$enterTime<data0$calendarCutoff, ])
}
