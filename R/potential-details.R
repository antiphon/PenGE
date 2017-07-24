#' Estimate Potential
#'
#' Details of the estimated interaction potential between two types
#'
#' @param fit Model fit returned by fitGlbin_CV
#' @param i Type 1 index, one of 1, ..., p
#' @param j Type 2 index, one of 1, ..., p
#'
#' @details The lasso-path coefficients of the interaction parameters between types i and j.
#' @export
potential <- function(fit, i, j, k = NULL) {
  Qpars <- fit$Qpars
  p <- fit$datainfo$p
  beta <- fit$fullfit$beta
  types <- fit$datainfo$types
  out <- list(model = "Multi-scale saturation potential", i=types[i], j=types[j], ij = c(i,j) )
  #
  if(is.null(k)) k <- 1:ncol(beta)
  if(i != j) {
    # solve k
    ii <- min(i,j) - 1
    jj <- max(i,j) - 1
    ijk <- (ii*(2*p - ii - 3) + 2*jj - 2)/2 + 1
    out$range <- Qpars$ranges2[[ijk]]
    out$sat <- Qpars$sat2[[ijk]]
    # which row:
    ri <- grep(paste0(types[ii+1], "v", types[[jj+1]]), rownames(beta))
    out$coef <- beta[ri, k]
  }
  else{
    out$range <- Qpars$ranges1[[i]]
    out$sat <- Qpars$sat1[[i]]
    ri <- grep(paste0("intra_", types[ii+1]), rownames(beta))
    out$coef <- beta[ri, k]
  }
  out
}
