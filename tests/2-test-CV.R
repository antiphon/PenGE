# Test CV
devtools::load_all(".")

library(spatstat)
data(lansing)

x <- lansing
if(!exists("Q")){
  p <- 6 # we have six species
  rv <- c(0.05, 0.1) # use this for everything
  r1 <- rep( list(rv) , p) # intra-type interaction ranges
  r2 <- rep( list(rv), p * (p-1)/2) # inter-type
  Q <- make_Q_stepper_multi(x=x, ranges1 <-r1 , ranges2 <- r2, verb = 1)
  # quads
  cv_border_expected_loss(0.1, nx = 2)
  K <- 2
  quads <- split_window(W = x$window, nx=K, ny=K)
}
# this will take a while
fit1 <-fitGlbin_CV(Q, verb=1, quads = quads, mc.cores=5)
fit2 <- fitGlbin3_CV(Q, verb=1, orthonormalise = FALSE, quads = quads, mc.cores = 5)
fit3 <- fitGlbin3_CV(Q, verb=1, orthonormalise = TRUE, quads = quads, mc.cores = 5)


fits <- list(fit1,fit2,fit3)

# Residuals
res1 <- residuals_cv(fit1, x = x, verb=1)
res2 <- residuals_cv(fit2, x = x, verb=1)
res3 <- residuals_cv(fit3, x = x, verb=1)

ress <- list(res1, res2, res3)


par(mfrow=c(3,4))
b <- NULL
for(i in 1:3){
  f <- fits[[i]]
  r <- ress[[i]]

  plot(f$fullfit$aic)
  plot(f$fullfit, only_pen = T)
  plot(f$lambda, ar)
  M <- imatrix(f, signed = TRUE, k = which.min(ar))
  image.imatrix(M, col = c("red","black","green","blue"), zlim=c(-1,2))
}




