# Dev spline fitting

devtools::load_all()

library(spatstat)
library(parallel)




# univariate
if(0) {
  Rmax <- 0.3
  knots <- seq(0, Rmax, l = 6)
  #x <- spatstat::setmarks(spatstat:::rpoint(1000), factor(1))
  x <- spatstat::setmarks(spatstat:::rThomas(50, 0.02, 5), factor(1))

  Q <- make_Q_splines_multi(x, knots1 = list(knots), knots2=NULL, sat1=20, auto_sat=F)

  z <- fitGlbin3_CV(Q)

  # Check fit

  r <- seq(0, Rmax, l = 100)

  e <- eval_estimated_potential(z, k =  100, at = r)

  plot(r, e, ylim = c(-1,1) * .5)
}


# Bivariate independent stuff
#x <- spatstat::setmarks(spatstat:::rpoint(1000), factor(1))
if(0){
  #x <- spatstat::rmpoint(c(10,10)*100)
  x <- spatstat::superimpose("1"=spatstat::rpoint(100), "2"=spatstat::rThomas(15,0.04,10))

  #
  Rmax <- 0.3
  knots <- seq(0, Rmax, l = 4)
  knots1 <- rep(list(knots), 2)
  knots2 <- list(knots)

  Q <- make_Q_splines_multi(x, knots1 = knots1, knots2=knots2, sat1=200, sat2=2000,
                            auto_sat=T, border_r = 0.2,
                            spline_order = 3, verb=1)

  #
  z <- fitGlbin3_CV(Q, verb=1, maxiter=5000, penalty = "grSCAD", lambda.min = 0.00001)

  # Check fit
  k <- 0.5 * which.min(z$fullfit$aic)
  r <- seq(0, Rmax, l = 100)

  par(mfrow=c(2,2))
  plot(z$fullfit$aic)
  abline(v=k)
  for(i in 1:2) for(j in i:2){
    e <- eval_estimated_potential(z, k =  k, at = r, i = i, j = j)
    plot(r, e, ylim = c(-1,1) * .5, main=paste(i,j, sep="-"))
  }

}


# Bivariate exactly correct model
L <- 100
bb <- cbind(0:1, 0:1)*L
Rmax <- 0.05 * L
knots <- seq(0, Rmax, l = 2 + 2)
knots_i <- knots[-c(1, length(knots))]
knots1 <- rep(list(knots), 2)
knots2 <- list(knots)

# Pontentials
set.seed(12345)
sp <- 3
K <- length(knots_i)  +1 + sp
xs <- seq(0, Rmax*1.2, l = 100)
B <- rgeyer:::spline_eval(xs/Rmax, p = sp, knots = knots_i/Rmax)
theta <- 5:0/2
ys <- B%*%theta
SAT <- 2
#plot(xs, ys, ylim = c(-2,2))

if(!exists("x")) {
  n <- 200
  th0 <- list(knots = c(0,1), theta = c(0,0,0,0), c = c(0,0,0,0) )
  th <- list(knots = knots, theta = theta, c = rep(SAT, length(theta)))
  th1 <- th0
  th2 <- th0
  th12 <- th
  xl <- rgeyer::rstepper_spline_biv(n = c(n,n), bbox = bb,
                                    theta1 = list(th1,th2), iter = 20000,
                                    theta2 = list(th12), dbg = 1)
  x <- as.ppp(xl[,1:2], W = as.owin(c(bb)))
  x <- setmarks(x, factor(xl[,3]))

  r <- seq(0, Rmax * 1.2, l = 100)
  G <- Kcross::pcf_cross_all_box(x, r)
  cat("sim done\n")
}

if(!exists("Q")) {
  quads <- split_window(x$window, nx = 2)
  set.seed(1)
  Q <- make_Q_splines_multi(x, knots1 = knots1, knots2=knots2, sat1=SAT, sat2=SAT,
                            rho_factor = 6,
                            auto_sat=F, border_r = Rmax*0.5,
                            spline_order = sp, verb=1)
}
# Steeper
# r1 <- seq(0, Rmax, l = 4)[-1]
# r11 <- rep(list(r1), 2)
# r12 <- list(r1)
# Q <- make_Q_stepper_multi(x, ranges1 = r11, ranges2=r12, sat1=2, sat2=2, rho_factor = 5,
#                           auto_sat=F, border_r = Rmax*0,
#                           spline_order = sp, verb=1)


#
fit <-  fitGlbin3_CV(Q, verb=1, quads = quads, mc.cores = 5)
fitM <- fitGlbin3_CV(Q, verb=1, quads = quads, mc.cores = 5, penalty = "grMCP")
fits <- list(lasso=fit, mcp = fitM)
res <- lapply(fits, residuals_cv, x = x, weight_by_n = TRUE)
# Check fit
kv <- lapply(res, function(rez) sapply(rez$by_resid, function(v) which.min(colMeans(v))))
M <- lapply(names(fits), function(m)  lapply(kv[[m]], function(k) imatrix(fits[[m]], k = k)))


############################
# Plots
par(mfrow=c(5,3), mar = c(3,4,2,2))

# CV and AIC per fit
for(re in names(fits)) {
  fiz <- fits[[re]]
  rez <- res[[re]]
  lam <- fiz$lambda
  xlx <- rev(range(lam))
  plot(lam, colMeans(rez$by_resid$inverse), ylab="inverse CV", xlim = xlx, "l", main = re,
       ylim = range(rez$by_resid$inverse))
  apply(rez$by_resid$inverse, 1, lines, x = lam, col ="gray90")
  abline(v=lam[kv[[re]][2]])
  plot(lam, colMeans(rez$by_resid$pearson), ylab="pearson CV", xlim = xlx, "l", main = "",
       ylim = range(rez$by_resid$pearson))
  abline(v=lam[kv[[re]][3]])
  apply(rez$by_resid$pearson, 1, lines, x = lam, col ="gray90")
  plot(lam, fiz$fullfit$aic, ylab="AIC", xlim = xlx, "l")

}
# potentials
for(i in 1:2) for(j in i:2){
  thh <- if(i == j) list(th1,th2)[[i]] else th12
  kn <- thh$knots; kni <- kn[-c(1,length(kn))]
  yy <- rgeyer::spline_eval(r/max(kn), p = sp, knots = kni/max(kn)) %*% thh$theta
  plot(xs, yy, ylim = c(-1,1)*2, main=paste("estimated ", i,j, sep="-"), "l", lwd=3)

  # estimate
  e <- eval_estimated_potential(fit, k =  kv$lasso[2], at = r, i = i, j = j)
  lines(r, e, col = 2)
  # max
  e <- eval_estimated_potential(fit, k =  100, at = r, i = i, j = j)
  lines(r, e, col = 2, lty = 3)
  # MCP
  e <- eval_estimated_potential(fitM, k =  kv$mcp[2], at = r, i = i, j = j)
  lines(r, e, col = 3)

  legend("bottom", col=c(1, 2,3), lty=1, c("true", "Lasso", "MCP"), ncol = 3, bty="n")
}
# Misc
plot(lam, fit$fullfit$df, xlim = xlx, "l", ylab="df", col=2)
lines(lam, fitM$fullfit$df, col=3)
# lines(lam, colSums(fit$fullfit$beta!=0), col=2)


plot(fit$fullfit, only_pen = T, main ="lasso", lambda_reverse = T, lambda_log = F)
plot(fitM$fullfit, only_pen = T, main ="mcp", lambda_reverse = T, lambda_log = F)
#
plot(r, G[1,1,], ylim = range(G), ylab = "pcf", "l")
lines(r, G[2,2,], col=2)
lines(r, G[1,2,], col=3)
plot(x, cols=1:2, legend = F)

