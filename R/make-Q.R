#' Multivariate Stepper Gibbs Model Pseudolikelihood Design Matrix Generation
#'
#' Construct the model matrix of multi-level multi-variate saturation model for penalised logistic regression pseudo-likelihood.
#'
#' @param x Point pattern data, preferrably a ppp-object. Say of p-types.
#' @param ranges1 Intra-type ranges, p vectors
#' @param ranges2 Inter-type ranges, p*(p-1)/2 vectors
#' @param sat1 Intra-type saturations, vectors (for Saturation model)
#' @param sat2 Inter-type saturations, vectors (for Saturation model)
#' @param sat1l More particular sat1, matrices (for Saturation model)
#' @param sat2l More particular sat2, matrices (for Saturation model)
#' @param auto_sat Determine saturations automatically?
#' @param covariates list of im-objects to use as covariates
#' @param penalise_covariates If FALSE, covariates are estimated without penalisation.
#' @param cov_noise_sd If > 0, add normal noise to covariates.
#' @param rho_factor Dummy pattern intensity w.r.t. data
#' @param rho_N Use this many dummies (vector pf length p)
#' @param dum_max Upper limit for dummy count, per type (after rho_factor and rho_N)
#' @param dum_min Lower limit for dummy count, per type (after rho_factor and rho_N)
#' @param dummies Optional: (x,y,m) matrix of dummy locations to use.
#' @param save_locations Store the locations in the output object? Used for CV
#' @param fatal_col0 Should we error on zero interaction columns (singularities).
#' @param border_r Border correction range.
#' @param verb Verbosity
#' @param ... Omitted.
#'
#' @details
#' [todo]
#'
#'
#' @return
#' A list suitable for inputting to penLogi-function
#'
#'
#' @import spatstat Matrix rgeyer
#' @export
#'
#'
make_Q_stepper_multi <- function(x, # data
                                 ranges1, ranges2, # spatial range vectors in lists
                                 sat1 = 1, sat2 = 1, sat1l = NULL, sat2l = NULL,# saturation model
                                 auto_sat = TRUE, # set saturation levels based on intensities and ranges
                                 covariates, penalise_covariates = TRUE,
                                 rho_factor = 4, rho_N=NULL, # count of dummies
                                 dum_max = 10000, dum_min = 100, # limit the dummies
                                 dummies = NULL, # or just use these dummies?
                                 save_locations = TRUE, # keep coordinates per row
                                 ...,
                                 verb=FALSE, #
                                 cov_noise_sd=0, # add a bit of noise to covariate images?
                                 fatal_col0=FALSE, # will error if singular columns
                                 border_r=0 # border correction, creates a subset-vector.
) {
  xpp <- x
  # Parse data
  x <- check_pp(xpp)
  bbox <- get_bbox(xpp)
  marks <- parse_marks(xpp)
  x <- get_coords(xpp)
  #
  V <- prod( apply(bbox, 2, diff) ) # window area
  m <- as.integer(mf <- factor(marks)) - 1
  mark_levels <- levels(mf)
  mark_types <- unique(m)
  #
  cat2 <- if(verb) cat else function(...) NULL
  #
  # check covariates
  ncovs <- 0
  if(!missing(covariates)) {
    if(!is.list(covariates)) stop("covariantes should be a list of im-objects.")
    ncovs <- length(covariates)
    cov_names <- names(covariates)
  }
  #
  # dimensions
  N <- table(m)
  p <- length(N)
  npairs <- p*(p-1)/2
  rn_intra <- sapply(ranges1, length)
  p_intra <- sum(rn_intra)
  cs_intra <- c(0, cumsum(rn_intra))
  # some checks
  if(p != length(rn_intra)) stop("ranges1 length does not match found types.")
  #
  # dummy intensities
  # limit
  #
  if(is.null(dummies)){
    Nd <- if(is.null(rho_N)) rho_factor * N else rho_N
    Nd <- round(Nd)
    Nd <- sapply(Nd, function(n) max(min(n, dum_max), dum_min)  )
    # Dummy sampling, use spatstat's straitified random
    nd <- ( ceiling(sqrt( sum(Nd) )) )
    B <- boundingbox(as.owin(c(bbox)))
    dums <- stratrand(B, nd, nd)
    dums <- cbind(dums$x, dums$y)[sample(1:length(dums$x), sum(Nd)),]
    d <- cbind(dums, sample(rep(1:length(Nd)-1, Nd)))
    dummiesm <- d
  }
  else{
    dummiesm <- dummies
  }
  o <- order(dummiesm[,3])
  dm <- dummiesm[o,3]
  dm <- as.integer(dmf <- factor(dm)) - 1
  dmark_types <- unique(dm)
  # check dummy labels match data
  if(!all(dmark_types %in% mark_types)) stop("dummies and data have different mark levels.")

  dummies <- dummiesm[o,-3]
  Nd <- c(table(dm))
  rho <- Nd/V
  #
  # Cumulative sum for indexing types, incl. dummies
  cs_data <- c(0, cumsum(Nd+N))
  #
  # Setup the intra-type interactions, covariates and intercept.
  Xdf_intercept <- Matrix(0, nrow = sum(N) + sum(Nd), ncol = p) # for speed in the first stage
  Xdf <- Matrix(0, nrow = sum(N) + sum(Nd), ncol = p * ncovs + p_intra, sparse=TRUE) # Intra
  # auxiliary data flag
  y <- NULL
  #
  if(save_locations) locations <- NULL
  #
  # border correction vector
  subset <- rep(TRUE, nrow(Xdf))
  # check saturation lists in case of non-symmetric sat's used.
  sat1l_given <- !is.null(sat1l)
  sat2l_given <- !is.null(sat2l)
  if(sat1l_given){
    if(length(sat1l) != p) stop("length(sat1l) != p")
  } else sat1l <- list()
  if(sat2l_given){
    if(length(sat2l) != npairs) stop("length(sat2l) != p * (p-1)/2")
  } else sat2l <- list()
  #
  # final bits and bobs
  first_names <- first_cov_names <- intra_names <- NULL

  #####################################################################################
  cat2("Intra components:\n")
  for(i in 1:p) {
    x_i <- x[m == i-1,, drop=FALSE]
    d_i <- dummies[dm == i-1,, drop=FALSE]
    # row location in the the data matrix
    idx <- 1:(N[i]+Nd[i]) + cs_data[i]
    #
    # determine saturation levels:
    if(auto_sat){
      if(sat1l_given) sats <- sat1l[[i]]
      else {
        dr2 <- diff(c(0, ranges1[[i]])^2)
        sats <- pmax(1, qpois(.99, dr2 * pi * nrow(x_i)/V))
      }
    }
    else{
      if(sat1l_given) sats <- sat1l[[i]]
      else sats <- rep(sat1, rn_intra[i]) # just a constant
    }
    # components: data to data
    xk <- stepper_components(x_i, to = NULL, r=ranges1[[i]], sat = sats, bbox = bbox, ...)
    # data to dummy
    dk <- stepper_components(x_i, d_i, r=ranges1[[i]], sat = sats, bbox = bbox, ...)
    #
    # Covariates
    if(ncovs) {
      require(spatstat)
      kc <- 1
      for(cov in covariates) {
        vx <- cov[data.frame(x_i), drop=FALSE]
        vd <- cov[data.frame(d_i), drop=FALSE]
        vc <- c(vx, vd)
        if(cov_noise_sd>0) vc <- vc + rnorm(length(vc), cov_noise_sd)
        Xdf[idx, (i-1) * ncovs + kc] <- vc
        kc <- kc + 1
      }
      first_cov_names <- c(first_cov_names, paste0(cov_names, "_", mark_levels[i]))
    }
    #
    # Store the values for this type
    #
    #
    # Border correction?
    if(border_r > 0) {
      bd <- bdistsxy(rbind(x_i, d_i), bbox)
      subset[idx] <- bd > border_r
    }
    #
    # intercept
    #browser()
    Xdf_intercept[idx, i] <- 1.0
    # interactions
    v <- rbind(xk, dk)
    # column index
    cidx <- ncovs * p + 1 + cs_intra[i] + (1:rn_intra[i]-1)
    Xdf[idx, cidx] <- v
    # binary observation indicators
    yi <- rep(1:0, c(N[i], Nd[i]))
    y <- c(y, yi)
    # used saturations
    sat1l[[i]]  <- sats
    # column names
    first_names <- c(first_names, paste0("intercept", mark_levels[i]))
    intra_names <- c(intra_names, paste0(paste0("intra_", mark_levels[i], "_", 1:rn_intra[i])))
    # if we want to keep the locations
    if(save_locations) locations <- rbind(locations,  cbind( rbind(x_i, d_i),  mark=i-1, is_data = yi ) )
    # done for type i
    cat2("      \r", i,"/", p)
  }
  cat2("\n")
  # combine with intercept matrix
  Xdf <- cBind(Xdf_intercept, Xdf)
  #####################################################################################
  # Inter-type interactions
  pair_names <- NULL
  #
  if(npairs > 0){
    k <- 1
    cat2("Inter components:\n")
    # some useful numbers
    rn_inter <- sapply(ranges2, length)
    p_inter <- sum(rn_inter)
    cs_inter <- c(0, cumsum(rn_inter))
    # check
    if(npairs != length(rn_inter)) stop("n. of pairs does not length of ranges2")
    # extend the data frame
    Xdf <- cbind(Xdf, Matrix(0, nrow=nrow(Xdf), ncol=p_inter, sparse=TRUE))
    for(i in 1:(p-1)) {
      x_i <- x[m == i-1,, drop=FALSE]
      d_i <- dummies[dm == i-1,, drop=FALSE]
      Ni <- N[i]+Nd[i]
      # location of type i points in the data frame
      idx_i <- 1:Ni + cs_data[i]
      for(j in (i+1):p){
        x_j <- x[m == j-1,, drop=FALSE]
        d_j <- dummies[dm == j-1,, drop=FALSE ]
        Nj <- N[j] + Nd[j]
        xij <- rbind(cbind(x_i,i-1), cbind(x_j,j-1))
        dij <- rbind(cbind(d_i,i-1), cbind(d_j,j-1))
        idx_j <- 1:Nj + cs_data[j]
        #
        # determine saturation levels and compute the interactions
        if(!auto_sat){
          # just a constant
          sats <- rep(sat2, rn_inter[k])
          # symmetric c's from i->j and j->i
          xk <- stepper_biv_components(xij, to = NULL, r=ranges2[[k]], sat = sats, bbox = bbox, ...)
          dk <- stepper_biv_components(xij, dij, r=ranges2[[k]], sat = sats, bbox = bbox, ...)
          sats <- cbind(ij=sats, ji=sats)
        }
        else{
          if(sat2l_given){
            sats <- sat2l[[k]]
            sats_ij <- sats[,1]
            sats_ji <- sats[,2]
          }
          else{
            dr2 <- diff(c(0, ranges2[[k]])^2)
            sats_ij <- pmax(1, qpois(.99, dr2 * pi * nrow(x_j)/V))
            sats_ji <- pmax(1, qpois(.99, dr2 * pi * nrow(x_i)/V))
          }
          sats_same <- all(sats_ij == sats_ji)
          # need to compute in two steps if different:
          xk_ij <- stepper_biv_components(xij, to = NULL, r=ranges2[[k]], sat = sats_ij, bbox = bbox, ...)
          xk_ji <- if(sats_same) xk_ij else stepper_biv_components(xij, to = NULL, r=ranges2[[k]], sat = sats_ji, bbox = bbox, ...)
          # dummies, same story
          dk_ij <- stepper_biv_components(xij, dij, r=ranges2[[k]], sat = sats_ij, bbox = bbox, ...)
          dk_ji <- if(sats_same) dk_ij else stepper_biv_components(xij, dij, r=ranges2[[k]], sat = sats_ji, bbox = bbox, ...)
          # collect right bits
          xk <- rbind(xk_ij[1:N[i],, drop=FALSE],  xk_ji[ N[i]+1:N[j],, drop=FALSE])
          dk <- rbind(dk_ij[1:Nd[i],, drop=FALSE], dk_ji[Nd[i]+1:Nd[j],, drop=FALSE])
          # for storing
          sats <- cbind(ij=sats_ij, ji=sats_ji)
        }
        # store saturation levels
        sat2l[[k]]  <- sats
        #
        # store values
        #browser()
        v <- rbind(xk, dk)
        # column index
        cidx <- p + ncovs * p + p_intra + 1 + cs_inter[k] + (1:rn_inter[k]-1)
        Xdf[idx_i, cidx] <- rbind( xk[1:N[i],, drop=FALSE], dk[1:Nd[i],, drop=FALSE] )
        Xdf[idx_j, cidx] <- rbind( xk[1:N[j]+N[i],, drop=FALSE], dk[1:Nd[j]+Nd[i],, drop=FALSE] )
        # column names
        pair_names <- c(pair_names,
                        paste0(paste0("inter_", mark_levels[i], "v", mark_levels[j] ),  "_", 1:rn_inter[k]))
        # done for pair (i,j)
        cat2("      \r", k, "/", npairs)
        k <- k + 1
      }
    }
  } # done pairs of types
  #######################################
  # done computing interactions
  #
  cat2("\n")
  ########################################
  # Name the columns for sanity
  colnames(Xdf) <- c(first_names, first_cov_names, intra_names, pair_names)
  #
  # offset
  H <- rep(-log(rho), N+Nd)
  #
  ############
  # The penalisation grouping vector. Denote by 0 not penalised:
  # dont penalise intercepts
  index <- rep(0, p)
  # covariates?
  cidx_covariate <- NULL
  if(ncovs){
    cidx_covariate <- p + 1:(p*ncovs)
    penc <- if(penalise_covariates) 1:(p*ncovs) else rep(0, p * ncovs) # no grouping
    index <- c(index, penc)
  }
  # Interactions:
  m <- max(index)
  #browser()
  # first order
  pen1  <- rep(m + 1:p, rn_intra)
  index <- c(index, pen1)
  cidx_intra <- p * (1+ncovs) + 1:length(pen1)
  # second order
  m2 <- max(index)
  pen2 <- if(npairs > 0) rep(m2 + 1:npairs, rn_inter) else NULL
  index <- c(index, pen2)
  cidx_inter <- p * (1+ncovs) + length(cidx_intra) + 1:length(pen2)
  #######################################
  # Collect to a list object
  out <- list(y=y, X=Xdf, offset=H, index = index, subset = subset,
              cidx = list(intercept = 1:p,
                          covariates = cidx_covariate,
                          intra = cidx_intra,
                          inter = cidx_inter
              ),
              datainfo = list(p = p, N=N, types = mark_levels)
  )
  #
  # store computation variables
  pars <- list()
  pars$border_r <- border_r
  pars$bbox <- bbox
  pars$ranges1 <- ranges1
  pars$ranges2 <- ranges2
  pars$sat1 <- sat1l
  pars$sat2 <- sat2l
  pars$auto_sat <- auto_sat
  pars$penalise_covariates <- penalise_covariates
  pars$dum_min <- dum_min
  pars$dum_max <- dum_max
  pars$rho_factor  <- rho_factor
  pars$ncovariates <- ncovs
  out$parameters <- pars
  #
  # all locations?
  if(save_locations) {
    out$locations <- locations
  }
  ############################
  # Check 0's
  out$nonsingular <- TRUE
  zeros <- sum(check_Q_sums(out)[,3])
  if(zeros>0){
    ooo <- if(fatal_col0) stop else warning
    ooo("some interactions are singular (0 columns)!")
    out$nonsingular <- FALSE
  }
  #done with it.
  out
}

###########
#####
# # Dummy sampling, use spatstat's straitified random
# sample_dummies <- function(Nd, bbox){
#   N <- sum(Nd)
#   nd <- ensure2vector( ceiling(sqrt(N)) )
#   B <- boundingbox(as.owin(c(bbox)))
#   dums <- stratrand(B, nd[1], nd[2])
#   dums <- cbind(dums$x, dums$y)[sample(1:length(dums$x), N),]
#   d <- cbind(dums, sample(rep(1:length(Nd)-1, Nd)))
#   d
# }
#
# ######
# Dummies on a grid, same for each type
dummy_grid <- function(types = 0:1, nx = 50, ny = nx, bbox){
  stepsx <- seq(bbox[1,1], bbox[2,1], l=nx)
  stepsy <- seq(bbox[1,2], bbox[2,2], l=ny)
  dumxy <- as.matrix( expand.grid(x = stepsx , y = stepsy) )
  nd <- nrow(dumxy)
  as.matrix( data.frame(dumxy, cbind(rep(types, each = nd))) )
}

# # Border distances
#' Distance xy to bbox
#'
bdistsxy <- function(xy, bbox){
  xd <- pmin( xy[,1]-bbox[1,1], bbox[2,1]-xy[,1])
  yd <- pmin( xy[,2]-bbox[1,2], bbox[2,2]-xy[,2])
  pmin(xd,yd)
}

#' Check Q-matrix singularities
#'
#'
check_Q_sums <- function(Q) {
  a <- Q$y==1
  z <- cbind(data = colSums(Q$X[a & Q$subset, Q$index>0, drop=FALSE]), dummy = colSums(Q$X[!a & Q$subset, Q$index>0, drop=FALSE]))
  z <- cbind(z, zero_column = 1 * (z[,1]==0 & z[,2] == 0))
  z
}
