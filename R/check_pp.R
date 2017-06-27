#' Internal checks
#'
check_pp <- function(x){
  if(!verifyclass(x, "ppp")) stop("Sorry, x not a ppp class")
  #if(!is.rectangle(x$window)) stop("Sorry, works only for rectangle windows.")
  x
}

get_bbox <- function(x) {
  w <- x$window
  cbind(w$xrange, w$yrange)
}

get_coords <- function(x) {
  as.matrix(coords(x))
}

#' Parse marks of a point pattern
#'
#'
parse_marks <- function(x){
  m <- x$marks
  if(length(m) != x$n |!is.factor(m)) stop("x should be marked with factors.")
  m
}
