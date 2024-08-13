#' Cut a dataset for analysis at a specified date
#'
#' @param data A time-to-event dataset returned from \link{sim_ph23}.
#' @param cut_time Date relative to start of randomization (\code{calendarTime} from
#' input dataset) at which dataset is to be cut off for analysis.

#'
#' @return A data frame ready for survival analysis including columns \code{calendarCutoff},
#' \code{survTimeCut}, \code{eventCut}.
#' @export
#'
#' @examples
#' d <- sim_ph23(n1_per_arm = 20, n2_per_arm = 30)
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
