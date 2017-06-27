#' Interaction Matrix for a Fitted Multivariate Multiscale Gibbs Model
#'
#' @param fit A model fit
#' @param k which penalty level, index in lambda
#' @param group_criteria per interaction step-function estimate vector, how to determine overall !=0
#' @param signed If TRUE, m[i,j]= 0 all 0 ; +1 all ranges non-negative; -1 all non-positive; +2 mix of pos. and negative (and possibly 0).
#'
#' @export
imatrix <- function(fit, k = NULL, group_criteria = any, signed = FALSE){
  # If CV wrapper given
  if(!is.null(fit$fullfit)){
    fit0 <- fit
    fit <- fit0$fullfit
    types <- fit0$datainfo$types
  }
  else{
    types <- NULL # TODO: parse
  }
  # determine the penalised
  idx <- fit$index
  # if k=NULL, use 0.5AIC
  if(is.null(k)) k <- round(0.5 * which.min(fit$aic[!is.na(fit$aic)]))
  #
  coef <- fit$beta[,k]
  nnz <- 1 * (coef != 0)   # non zero's
  cn <- rownames(fit$beta)
  idx1 <-  grep("intra_", cn)  # which are intra
  idx2 <- grep("inter_", cn)
  idx1_grouped <- split(idx1, idx[idx1])
  idx2_grouped <- split(idx2, idx[idx2])

  #
  if(signed){
    p1 <- sapply(idx1_grouped, function(i) group_criteria(nnz[i] == 1) * sign_it(coef[i]) )
    p2 <- sapply(idx2_grouped, function(i) group_criteria(nnz[i] == 1) * sign_it(coef[i]) )
  }
  else{
    p1 <- sapply(idx1_grouped, function(i) group_criteria(nnz[i] == 1)   )
    p2 <- sapply(idx2_grouped, function(i) group_criteria(nnz[i] == 1)   )
  }
  # make the matrix
  M <- diag(p1)
  # inter
  M[lower.tri(M)] <- p2
  M <- t(M)
  M[lower.tri(M)] <- p2
  # Add var names
  colnames(M) <- rownames(M) <- types
  # for plotting, set a class
  #class(M) <- c("imatrix", is(M))
  M
}


#' Image for imatrix
#'
#' @param x imatrix
#' @param colors to use
#'
#' @details
#'
#' Colors for binary: (0,1)
#'
#' Colors for signed: (-1,0,1,2)
#'
#' @export
image.imatrix <- function(x, cols = c("black", "white"), ...) {
  p <- ncol(x)
  image(1:p, 1:p, x, ann=F, axes=F, asp=1, col = cols, ...)
  axis(1, 1:p, rownames(x), ...)
  axis(2, 1:p, rownames(x), ...)
}



# sign of a coefficient vector
sign_it <- function(v) {
  haspos <- any(v>0)
  hasneg <- any(v<0)
  if(haspos) if(hasneg) 2 else 1 else if(hasneg) -1 else 0
}

#' Covariate estimate Matrix
#'
#'
covariate_matrix <- function(fit, pref="covariate") {
  beta <- fit$beta
  idxc <- grep(pref, rownames(beta))
  if(length(idxc)==0) return(NULL)
  nan <- rownames(beta)[idxc]
  # split by cov type
  covn <- unique(gsub("_[0-99,A-Z]*$", "",nan))
  tab <- NULL
  for(cn in covn) tab <- rbind(tab, beta[idxc[grep(cn, nan)]]   )
  rownames(tab) <- covn
  colnames(tab) <- gsub("[^0-99,A-Z]", "", colnames(tab))
  tab
}

#
# covmatrix_from_out <- function(out) {
#   beta <- out$theta0
#   idxc <- grep("covariate", names(beta))
#   if(length(idxc)==0) return(NULL)
#   nan <- names(beta)[idxc]
#   # split by cov type
#   covn <- unique(gsub("_[0-99]*$", "",nan))
#   tab <- NULL
#   for(cn in covn) tab <- rbind(tab, beta[idxc[grep(cn, nan)]]   )
#   rownames(tab) <- covn
#   colnames(tab) <- gsub("[^0-99]", "", colnames(tab))
#   tab
# }
#
#
#
# # For the "out" object (older format)
# imatrix_from_out <- function(out, group_criteria = any, signed = FALSE){
#     #
#     coef <- out$theta
#     nnz <- 1 * (coef != 0)   # non zero's
#     cn <- names(coef)
#     idx <- out$index # the grouping vector
#     idx1 <-  grep("intra_", cn)  # which are intra
#     idx2 <- grep("inter_", cn)
#     idx1_grouped <- split(idx1, idx[idx1])
#     idx2_grouped <- split(idx2, idx[idx2])
#     #
#     if(signed){
#       p1 <- sapply(idx1_grouped, function(i) group_criteria(nnz[i] == 1) * sign_it(coef[i]) )
#       p2 <- sapply(idx2_grouped, function(i) group_criteria(nnz[i] == 1) * sign_it(coef[i]) )
#     }
#     else{
#       p1 <- sapply(idx1_grouped, function(i) group_criteria(nnz[i] == 1)   )
#       p2 <- sapply(idx2_grouped, function(i) group_criteria(nnz[i] == 1)   )
#     }
#     # make the matrix
#     M <- diag(p1)
#     # inter
#     M[lower.tri(M)] <- p2
#     M <- t(M)
#     M[lower.tri(M)] <- p2
#     M
#   }
