## ---- results=FALSE, include=FALSE---------------------------------------
library(PenGE)
library(spatstat)

## ---- eval=FALSE---------------------------------------------------------
#  library(PenGE)
#  library(spatstat)

## ---- fig.width=7, fig.height=6------------------------------------------
data(lansing)
plot(split(lansing), cex=.5)

## ------------------------------------------------------------------------
p <- 6 # we have six species
rv <- c(0.025, 0.05) # use these for everything
r1 <- rep( list(rv) , p) # intra-type interaction ranges
r2 <- rep( list(rv), p * (p-1)/2) # inter-type

## ---- fig.width = 6, fig.height=6----------------------------------------
M <- imatrix(fit, signed = TRUE)
image.imatrix(M, col = c("red","black","green","blue"), cex.axis = 0.7)

## ---- fig.width = 6, fig.height=6----------------------------------------
# remove the intercepts which are not penalised, and divide by number of ranges per interaction
nnzero <- (colSums(fit$fullfit$beta != 0) - p )/ 2
plot(fit$lambda, nnzero, xlim = rev(range(fit$lambda)), main ="number of non-zero interactions")

## ------------------------------------------------------------------------
r_border <- 0.05 # define a border zone radius, should be the max. of the model's interaction radii.
win_side_len <- 1
nmax <- cv_border_solve_split(R_frac = r_border/win_side_len, acceptable_loss = 0.5)
nx <- round( sqrt(nmax) ) # nx times nx regular quadrature
quads <- split_window(lansing$window, nx = nx, ny = nx)

## ---- fig.width=6, fig.height=6------------------------------------------
plot(unmark(lansing), cex=.3, main = "CV quadrats")
for(q in quads) plot(erosion(w = q, r = r_border), add=T)

## ---- fig.width=7, fig.height=4------------------------------------------
inv <- Rest$by_resid$inverse
pen <- cvfit$lambda 

par(mfrow=c(1,2))

# The residual per penalty for each CV quadrat
plot(pen, inv[1,], "l", xlim = rev(range(pen)), ylim=range(inv), main="Inverse residual SS by quadrat")
for(i in 2:nx^2) lines(pen, inv[i,], col=i)
# Estimated Risk:
risk_inv <- Rest$ave$inverse
plot(pen, log(risk_inv), "l", xlim = rev(range(pen)), main = "Mean inverse residual SS")

## ---- fig.width = 6, fig.height=6----------------------------------------
kopt <- which.min(risk_inv)
Mcv <- imatrix(cvfit, k = kopt, signed = TRUE)
image.imatrix(Mcv, cols = c(2,1,3,4),  zlim = c(-1,2), cex.axis=.7)

## ------------------------------------------------------------------------
print(Mcv)

## ------------------------------------------------------------------------
print(Q$parameters$sat2[[1]])
print( Q1$parameters$sat2[[1]] )

## ------------------------------------------------------------------------
Q$parameters$ranges2[[1]]

## ------------------------------------------------------------------------
potential(cvfit, i = 1, j = 3, k = kopt)

## ------------------------------------------------------------------------
cov1 <- as.im(function(x,y) (x-.5)^2, W = lansing$window)
cov2 <- as.im(function(x,y) (y-.5)^2, W = lansing$window)
covs <- list(covx = cov1, covy = cov2)

## ------------------------------------------------------------------------
Qc <- make_Q_stepper_multi(x=lansing, ranges1 = r1, ranges2 = r2, border_r = 0.05, covariates = covs, penalise_covariates = FALSE)

## ------------------------------------------------------------------------
cvfit_c <- fitGlbin_CV(Qc, quads = quads, mc.cores = 5)
Rc <- residuals_cv(cvfit_c, x = lansing, mc.cores = 5, covariates = covs)

## ------------------------------------------------------------------------
kopt_c <- which.min(Rc$ave$inverse)
cov_est <- covariate_matrix(cvfit_c, k = kopt_c, pref = "cov")
print(cov_est, digits=2)

