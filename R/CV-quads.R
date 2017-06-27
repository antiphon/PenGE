#' Solve What Kind of CV Split Can Be Done
#'
#' How many quadrats can we have given border correction range and tolerated loss of data?
#'
#' @param R_frac  border range / window side length - fraction
#' @param loss  max tolerated lost data, between (0, 1)
#'
#' @export
cv_border_solve_split <- function(R_frac = 0.02, acceptable_loss = 0.5){
  n <- (1 - sqrt(1-acceptable_loss))/(2 * R_frac)
  n
}

#' Loss of Data
#'
#' How much data is wasted due to border correction in a specific nx x nx split of a square window
#'
#' @param R_frac border range / window side length - fraction
#' @param nx nx x nx split
#' @export
cv_border_expected_loss <- function(R_frac = 0.05, nx = 5) {
  1-(1-2*R_frac*nx)^2
}

#' Split Obs Window to Small Windows
#' @param bbox bbox in range-column form
#' @param W owin-object
#' @param ... passed to spatstat::quadrats
#'
#' @details
#' If W not given,  try to convert bbox to owin.
#' @export
#' @import spatstat
#'
split_window <- function(bbox, W, ...){
  if(missing(W)) W <- as.owin(c(bbox))
  QW <- quadrats(W, ...)
  lapply(1:QW$n, function(i) spatstat::as.owin(QW[i]))
}


