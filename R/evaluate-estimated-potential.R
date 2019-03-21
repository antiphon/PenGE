#' Evaluate Estimated Potential
#'
#' @param fit a group penalised fit object
#' @param k which of the estimates to use
#' @param i species 1
#' @param j species 2 (can be i)
#' @param r range over which to evaluate the potential
#' @export

eval_estimated_potential <- function(at, fit, k, i = 1, j = 1) {
  #
  if(!is.null(fit$fullfit)){
    types <- fit$datainfo$types
    fit0 <- fit
    ffit <- fit$fullfit
  }
  else{
    stop("can't parse types")
  }
  # determine the penalised
  idx <- ffit$index_out
  if(is.null(idx)) idx <- ffit$index # old version 1
  # if k=NULL, use 0.5AIC
  aic <- ffit$aic
  if(is.null(k)) k <- round(0.5 * which.min(!aic[is.na(aic)]))

  ii <- i
  i <- min(ii,j)
  j <- max(ii,j)
  ntypes <- length(types)
  if(is.numeric(i)) i <- types[i]
  if(is.numeric(j)) j <- types[j]
  ni <- which(i == types)
  nj <- which(j == types)
  # assuming stepper function. parse basis
  Qp <- fit$Qpars
  # fix this
  NN <- matrix(0, ntypes, ntypes)
  NN[upper.tri(NN)] <- 1:(ntypes*(ntypes-1)/2)
  nij <- NN[ni,nj]
  #

  # find coefficients
  B <- ffit$beta[, k]
  if(i == j) b <- B[  grep(paste0("intra_", i), names(B))  ]
  else b <- B[  grep(paste0("inter_", i, "v", j), names(B))  ]
  b <- unname(b)
  # Here depends on model
  type <- Qp$type
  if(is.null(type)) type <- "stepper" # old version support
  #browser()
  if(type == "stepper") {
    #
    if(i == j) rv <- Qp$ranges1[[ni]]
    else rv <- Qp$ranges2[[nij]]
    # evaluate
    rve <- c(0,rv)
    s <- findInterval(at, rve)
    # compute
    out <- c(b,0)[s]
  }
  else if(type == "genstepper") {
    if(i == j) rmat <- Qp$ranges1[[ni]]
    else rmat <- Qp$ranges2[[nij]]
    # evaluate
    V <- apply(rmat, 2, function(ab) c(0,1,0)[1+findInterval(at, ab)]   )
    # compute weighted sum
    out <- c(V %*% b )
  }
  else if(type == "splines") {
    if(i == j) kn <- Qp$knots1[[ni]]
    else kn <- Qp$knots2[[nij]]
    nk <- length(kn)
    so <- Qp$spline_order
    np <- length(b)  #length(kn) + so - 1 # in case dropped the last wonky basis function
    V <- sapply(1:np, rgeyer::spline_eval, x=at/kn[nk], p = so, knots = kn[-c(1, nk)]/kn[nk], ver=3)
    out <- c( V%*%b  )
  }
  # done
  out
}



