#' Adjusted p value for testing H_J using Simes method
#'
#' This functions provides the adjusted p value for testing a familywise hypothesis H_J based on the p values of individual raw p values in J.
#'
#' @param  p A vectory of individual raw p values in J.
#'
#' @return Adjusted p value using simes procedure
#'
#' @examples
#' #Example (1):
#'
#' p = c(0.01, 0.02, 0.03, 0.013)
#' simes(p)
#'
#' @export
simes = function(p) {
    r = rank(p)
    ans = min(min(length(p) * p/r), 1)
    return(ans)
}

