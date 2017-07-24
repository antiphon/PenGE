#' Compute Spatial Residuals from CV Fits
#'
#' @param fit cv-fit
#' @param x original data pattern used to generate the design matrix
#'
#'@export
residuals_cv <- function(fit, x,
                         dum_min = 1000, dum_max = 5*dum_min, rho_factor = 10,
                         dum_grid_dimyx = NULL, # in case given use grid dummies
                         verb = FALSE, mc.cores = 1, weight_by_n = FALSE, ...) {
  if(length(fit$quads) < 1|length(fit$cvfits) < 1) stop("CV fit not detected (fit$quads or fit$cvfits of zero length).")
  # the lambda
  l <- fit$lambda
  nl <- length(fit$lambda)
  # cv fits
  cvfits <- fit$cvfits
  if(verb){
    require(looptimer)
    t0 <- looptimer(n = length(cvfits))
  }
  # cv quads
  quads <- fit$quads
  # setup
  Qpars <- fit$Qpars
  # we need data to compute the integrals
  xpp <- x
  # Parse data
  x <- check_pp(xpp)
  bbox <- get_bbox(xpp)
  marks <- parse_marks(xpp)
  pat <- get_coords(xpp)
  types <- unique(marks)
  pat <- pat[,1:2]

  # compute risk over CV splits:
  nq <- length(cvfits)
  nblocks <- ceiling(nq/mc.cores)
  bl <- rep(1:nblocks, each=mc.cores)[1:nq]
  blocks <- split(1:nq, bl)

  # for one quadrat
  for_one <- function(k) {
    fitk <- cvfits[[k]]
    if(!is.list(fitk)){
      warning(paste0("CV fit ", k, "is not of correct type."))
      NULL
    }
    else{
      W <- quads[[k]]
      # test data pattern
      xin <- inside.owin(pat[,1], pat[,2], W)
      if(sum(xin) == 0){
        warning(paste("quadrat", k, "is empty"))
        res <- NULL
        count <- 0
      }
      xm <- list(x=pat[xin,], bbox = cbind(W$xrange, W$yrange))
      m <- marks[xin]
      count <- length(m)
      # check if all types present
      z <- unique(m)
      mis <- setdiff(types,z)
      # add a point of the missing type to keep data structure intact for downstream computations.
      # add it outisde the window so it will not influence anything.
      for(mi in mis) {
        #xm$x <- rbind(xm$x, as.matrix(coords(runifpoint(1, W))))
        xnew <- c(0, W$xrange[2] + Qpars$border_r + 1)
        xm$x <- rbind(xm$x, xnew)
        m <- c(m, mi)
      }
      # check if grid dummies for integration
      if(!is.null(dum_grid_dimyx)){
        dummies <- dummy_grid(types, nx = dum_grid_dimyx[2], ny = dum_grid_dimyx[1],
                              bbox = xm$bbox)

      } else{
        dummies <- NULL
      }
      #
      # compute new Q in the quad window, idea is to estimate the integral by dummy mean.
      xpp <- ppp(xm$x[,1], xm$x[,2], marks = m, window = W)
      Qk <- make_Q_stepper_multi(xpp, ranges1 = Qpars$ranges1, ranges2 = Qpars$ranges2,
                                 rho_factor = rho_factor,
                                 auto_sat = Qpars$auto_sat,
                                 sat1l = Qpars$sat1, sat2l = Qpars$sat2,
                                 dum_min = dum_min, dum_max = dum_max,
                                 dummies = dummies,
                                 save_locations = TRUE,
                                 ...)
      #
      # compute residuals
      res <- resid_by_Q(fitk, Qk, W = erosion.owin(W, Qpars$border_r), ...)
      res_ave <- average_residuals_over_types(res)
      # add NA in case AIC stop rule was used:
      ni <- nrow(res_ave)
      if(ni < nl)  res_ave[(ni+1):nl,] <- NA
      list(res=res_ave, n=count, res_by_type = res)
    }
  }# eo for_one(...)

  Rk <- list()
  # run through the blocks
  t0 <- looptimer(n = nblocks, prefix = "* resid CV")
  counts <- NULL
  res_by_type <- list()
  for(b in blocks) {
    resl <- mclapply(b, for_one, mc.cores = mc.cores)
    for(r in resl) for(n in names(r$res)) Rk[[n]] <- rbind(Rk[[n]], r$res[[n]])
    counts <- c(counts, sapply(resl, getElement, "n"))
    res_by_type[b] <- lapply(resl, getElement, "res_by_type")
    if(verb) print(t0 <- looptimer(t0))
  }

  # Estimate of risk per residual type as the mean of squared residuals
  #
  ww <- if(weight_by_n) counts/sum(counts) else 1/length(counts)
  #browser()
  Rest <- data.frame( sapply(Rk, function(Rkl) colSums(ww*(Rkl^2), na.rm=TRUE)) )
  #
  # Average Risk over residual types, weighted by standard deviation:
  Rsd <- sapply(Rk, function(Rkl) apply(Rkl, 2, sd, na.rm=TRUE))
  Rsd <- Rsd/rowSums(Rsd)
  Rest[["mean"]] <- rowSums(Rest*Rsd)
  #
  # Average of inv and pearson residuals only:
  Rsd <- Rsd[,-1]/rowSums(Rsd[,-1])
  Rest[["mean_no_raw"]] <- rowSums(Rest[,-1]*Rsd)
  #
  # done
  Rest
  out <- list(ave=Rest, by_resid = Rk, counts = counts, res_by_type = res_by_type)
  class(out) <- c("cvresiduals")
  out
}



#' Compute various h-residuals from group lasso fits from the Q object.
#'
#' The integral is estimated using MC expectation over a set of dummies.
#' @param fit glbinc fit
#' @param Q corresponding design object
#' @param W window for border correction
#' @export
resid_by_Q <- function(fit, Q, W=NULL, ...) {
  if(is.null(Q$locations)) stop("Q needs to have locations saved (save_locations = TRUE)")
  if(is.null(W))
    W <- as.owin(c( Q$bbox + c(1,-1) * Q$parameters$border_r) ) # border correction!
  # the prediction points
  these <- inside.owin(Q$locations[,1], Q$locations[,2], W)
  # the prediction types
  types <- unique(Q$locations[these, 3])
  # the estimates
  beta <- fit$beta[, !is.na(fit$aic), drop = FALSE]
  # re-level intercepts to match indicators:
  for(i in seq_along(types)[-1]) beta[i,] <- beta[i,] + beta[1,]

  # compute papangelou
  lam <- as.matrix(  Q$X[these,] %*% beta  )
  #
  V <- area(W)
  raw <- inv <- pearson <- NULL
  # Residual per type
  for(i in types) { # per type.
    lda <- exp(lam[Q$y[these] == 1 & Q$locations[these, 3] == i , , drop = FALSE ])
    ldu <- exp(lam[Q$y[these] == 0 & Q$locations[these, 3] == i , , drop = FALSE ])
    # raw residuals
    raw <- rbind(raw, nrow(lda) - apply(ldu, 2, mean) * V)
    # inverse
    inv <- rbind(inv, colSums(1/lda) - apply(ldu>0, 2, mean) * V)
    # pearson
    pearson <- rbind(pearson, colSums(1/sqrt(lda)) - apply(sqrt(ldu), 2, mean) * V)
  }
  #
  n <- c( table(Q$locations[Q$y == 1, 3]) )
  list(raw=raw, inv=inv, pearson=pearson, n = n)
}







#' Sum or average residuals over types
#' @export
average_residuals_over_types <- function(res, weight_by_count = FALSE,
                                         sqrt_weights = FALSE, trim = 0){
  if(trim == 0){
    n <- res$n
    if(sqrt_weights) n <- sqrt(n)
    if(weight_by_count) w <- n/sum(n) else w <- 1
    # total sums over types, weight by point count.
    data.frame(raw = colSums(w * res$raw),
               inverse = colSums(w * res$inv),
               pearson = colSums(w * res$pearson))
  }
  else{
    mf <- function(x) apply(x, 2, mean, trim = trim, na.rm=TRUE)
    data.frame(raw = mf(res$raw),
               inverse = mf(res$inv),
               pearson = mf(res$pearson))
  }
}
