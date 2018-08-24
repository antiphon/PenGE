#' Orthonormalisation for the design matrix
#'
#' @param X design matrix
#' @param index grouping
#' @export

orthonormalise <- function(X, index) {
  # standardize
  S <- standardize(X)
  # orthogonalise groupwise
  OG <- orthogonalise(S$X, index)
  list(X = OG$X,
       centers = S$centers, scales  = S$scale,
       QL = OG$QL,
       index = index )
}

#' Un-orthonormalisation of Coefficients
#'
#' @param X design matrix
#' @param ON orthonormalisation setup object from 'orthonormalise'
#' @export
unorthonormalise <- function(beta, ON) {
  bs <- unorthogonalise(beta, ON)
  unstandardize(bs, cs = ON)
}


#' Orthogonalise Grouped Columns of X
#'
#' @param X matrix of coefficients (rows: coefs, cols: penalisations)
#' @param index grouping
#'
#' @details
#'
#' For columns with index == 0, do nothing. For each group of columns defined by index == j,
#' transform the matrix such that t(Z)%*%Z/n = I, Z=X[,index==j], n = nrow(X), I = diag(1, n).
#'
#' @export
orthogonalise <- function(X, index)  {
  n <- nrow(X)
  p <- ncol(X)
  XX <- matrix(0,  n, p)
  colnames(XX) <- colnames(X)
  XX[, index == 0] <- X[, index == 0]
  J <- max(index)
  QL <- vector("list", J)
  for( g in 1:J ) {
    i <- which( index == g )
    deco <- svd(X[, i, drop=FALSE], nu = 0)
    ql <- t( t(deco$v) * sqrt(n) / deco$d )
    #microbenchmark::microbenchmark(ql <- sweep(deco$v, 2, sqrt(n)/deco$d, "*"))
    QL[[g]] <- ql
    XX[, i] <- X[, i] %*% ql
  }

  OG <- list(X=XX, QL = QL, index = index)
  OG
}

#' Unorthonogalise Coefficients
#'
#' @param beta matrix of coefficients (rows: coefs, cols: penalisations)
#' @param OG orthogonalisation object from 'orthogonalise'
#' @import Matrix
#' @export
unorthogonalise <- function(beta, OG) {
  # check if intercept was added
  p0 <- length(OG$index)
  p1 <- nrow(beta)
  if(abs(p0-p1) > 1) stop("dimension mismatch")
  intercept <- 1 * (p1 - p0 == 1)
  idx <- c( if(intercept) 0 else NULL, OG$index)
  # use a single block diagonal matrix
  b <- matrix(0, nrow(beta), ncol(beta) )
  rownames(b) <- rownames(beta)
  B <- Matrix::bdiag( OG$QL )
  not <- which(idx == 0)
  if(length(not)) {
    b[not, ] <- beta[not,]
    b[-not, ] <-  as.matrix(  B %*% beta[-not,]  )
  }
  else
    b <-  as.matrix(  B %*% beta  )
  b
}



#' Standardize X
#'
#' Same as running 'scale' apart from result object being a list.
#' @param X design matrix. Will freak out if contains constant columns.
#'
#' @export
standardize <- function(X) {
  n <- nrow(X)
  a <- apply(X, 2, function(z) c(mean(z), sd(z)))
  cent <- a[1,]
  scal <- a[2,]
  if(any(scal == 0)) stop(paste0("Constants columns in X: ", paste0(which(scal == 0), collapse=","), "" ))
  Y <- t((t(X) - cent) / scal)
  list(X = Y, centers = cent, scales = scal)
}

#' Unstandardise Coefficients
#'
#' @param beta coef matrix (1 row per coef)
#' @param centers col means
#' @param scales col sds
#' @param cs alternative list with elements 'centers' and 'scales'
#'
#' @export
unstandardize <- function(beta, centers, scales, cs) {
  if(!missing(cs)) {
    centers <- cs$centers
    scales <- cs$scales
  }
  # check if intercept was added
  p0 <- length(centers)
  p1 <- nrow(beta)
  if(abs(p0-p1)>1) stop("dimension mismatch")
  intercept <- 1 * (p1-p0 == 1)
  std_i <- (1 + 1*intercept):nrow(beta)
  b <- beta
  b[std_i,] <- b[std_i,]/scales
  shift <- c( crossprod(centers, b[std_i,]) )
  b[-std_i,] <- b[-std_i,] - shift # adjust the global intercept
  b
}


