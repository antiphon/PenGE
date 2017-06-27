#' Group Lasso Penalised Logistic Regression Estimation
#'
#' Uses the simplified LCD
#'
#' @param Q Model object
#' @param add_intercept Add overall intercept?
#' @param dfmax stop if dfmax significant terms reached
#' @param lambda penalty vector. if missing, will be computed
#' @param lambda.min If lambda is to be computed, min. lambda will be max.lambda * lambda.min
#' @param nlambda steps to generate for lambda vector
#' @param log_linear_lambda log-linear or linear
#' @param use_glm_in_lmax use R-s own glm function (non-sparse) to determine max. lambda
#' @param quads The cross validation quadrats
#' @param mc.cores mclapply mc.cores
#' @param sparse use Sparse Matrix
#' @param verb Verbosity
#' @param ... passed on to
#'
#' @details
#'
#' Basically just a wrapper for package glbinc for estimating group penalised binomial regression. If quads are given, fits the CV-models as well. Here the mc.cores -parameter is used.
#'
#'
#' @import looptimer glbinc parallel
#' @export

fitGlbin_CV <- function(Q,
                        dfmax = ncol(Q$X) * 1.1,
                        lambda = NULL,
                        lambda.min = 0.001,
                        nlambda = 100,
                        log_linear_lambda = FALSE,
                        verb = 0,
                        add_intercept=TRUE,
                        quads = NULL,
                        mc.cores = 3,
                        sparse = FALSE,
                        use_glm_in_lmax = TRUE,
                        ...) {

  # main subset operation here
  ok <- Q$subset
  X <- if(sparse) Q$X[ok,] else as.matrix( Q$X[ok,] )
  y <- Q$y[ok]
  o <- Q$offset[ok]
  # the penalisation index vector
  index <- Q$index
  # if we add the intercept, drop the first to reduce singularities
  if(add_intercept){
    X <- X[,-1]
    index <- index[-1]
  }
  #
  # Determine penalty vector
  if(is.null(lambda)){
    lmax <- max_penalty(X, y, index, o, use_glm = use_glm_in_lmax)
    lvec <- seq(lmax, lmax * lambda.min, l = nlambda) # linear
    if(log_linear_lambda) lvec <- exp( seq(log(lmax), log(lmax * lambda.min), l = nlambda ) )
  }
  else lvec <- sort(lambda, decreasing=TRUE)
  #
  fitter <- if(sparse) glbin_lcd_c_sparse else glbin_lcd_c
  #
  # add the full window as quad #1
  quads <-  append(list(0), quads)
  #
  # The full and CV fits in parallel blocks, using mc.cores:
  nq <- length(quads)
  nblocks <- ceiling(nq/mc.cores)
  bl <- rep(1:nblocks, each=mc.cores)[1:nq]
  blocks <- split(1:nq, bl)
  #
  #
  # function to fit the model in a given window
  fitone <- function(quad){
    z <- if(!is.list(quad)) rep(TRUE, sum(ok)) else !inside.owin(Q$locations[ok,1], Q$locations[ok,2], quad)
    fit <- fitter(X=X[z,], y=y[z], index = index, offset=o[z],
                  lambda = lvec,
                  dfmax = dfmax, verb = verb,
                  add_intercept = add_intercept,
                  ...)
    fit # done
  }
  #
  # lets go
  t0 <- looptimer(n = nblocks, prefix = "* CV ", when_ready = FALSE)
  fitl <- list()

  for(bl in blocks) {
    fitbl <- mclapply(quads[bl], fitone, mc.cores=mc.cores)
    fitl[bl] <- fitbl
    if(verb>0) print(t0 <- looptimer(t0))
  }
  # done
  # return the results:
  out <- list(fullfit = fitl[[1]], cvfits = fitl[-1], lambda = lvec, quads = quads[-1],
          sparse = sparse, Qpars = Q$parameters, cidx = Q$cidx, datainfo = Q$datainfo)
  class(out) <- c("cvfit", is(out))
  out
}



