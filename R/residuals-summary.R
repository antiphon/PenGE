#' Summary for CV Residuals
#'
#' @export
summary.cvresiduals <- function(x, ...) {
  # the minimums:
  apply(x$ave, 2, which.min)
}

