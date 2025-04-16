#' Cut a dataset for analysis at a specified event count
#'
#' @param data A time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param targetEvents Event count at which data cutoff is to be made.
#'
#' @return A data frame ready for survival analysis including columns \code{calendarCutoff},
#' \code{survTimeCut}, \code{eventCut}.
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 100, prop_S = 0.5, ratio = 1)
#' dcut <- cut_by_event(d, 10)
cut_by_event <- function (data, targetEvents) {
  data0 = data
  data0.order <- data0[order(data0$calendarTime), ]
  data.event <- data0.order[data0.order$event == 1, ]
  data.event$event.seq <- seq.int(nrow(data.event))
  data0$calendarCutoff = data.event$calendarTime[data.event$event.seq ==
                                                   targetEvents]
  data0$survTimeCut = ifelse(data0$calendarTime <= data0$calendarCutoff,
                             data0$survTime, data0$calendarCutoff - data0$enterTime)
  data0$eventCut = ifelse(data0$calendarTime <= data0$calendarCutoff,
                          data0$event, 0)
  return(data0[data0$enterTime<data0$calendarCutoff, ])
}
