# Check the potential evaluation tool

devtools::load_all()

library(spatstat)
data(lansing)

x <- lansing
if(!exists("Q")){
  p <- 6 # we have six species
  rv <- c(0.05, 0.1, 0.2) # use this for everything
  r1 <- rep( list(rv) , p) # intra-type interaction ranges
  r2 <- rep( list(rv), p * (p-1)/2) # inter-type
  Q <- make_Q_stepper_multi(x=x, ranges1 <-r1 , ranges2 <- r2, verb = 1, sat1 = 1000, sat2= 1000)
  # quads
  cv_border_expected_loss(0.1, nx = 2)
  K <- 2
  quads <- split_window(W = x$window, nx=K, ny=K)
}


#fit <- fitGlbin3_CV(Q, verb=1, orthonormalise = TRUE)
fit <- fitGlbin3_CV(Q, verb=1, penalty = "grLasso", orthonormalise = TRUE, maxiter = 5000)


#

r <- seq(0, .22, l = 100)
G <- Kcross::pcf_cross_all_box(x, r)


k <- 50
maxk <- ncol(fit$fullfit$beta) # use the full estimates.
M1 <- imatrix(fit, maxk, signed = F)
M <- imatrix(fit, k, signed = T) * M1
not0 <- sum(diag(M!=0))+sum(M[upper.tri(M)]!=0)


par(mfrow=c(not0+1,2), mar=c(0,2,1,0)+.5)
for(i in 1:6)
  for(j in i:6)if(M[i,j]){
    z <- eval_estimated_potential(fit, maxk, i, j)
    plot(r, z, type="l", ann=F, xaxt="n"); box(); title(main = paste(i,j))
    abline(h = 1, col=2)
    plot(r,G[i,j,], type="l", ann=F, axes=F, ylim = c(0,2)); box(); title("pcf")
    lines(r, exp(z), col=3)
    abline(h = 1, col=2)
}

image.imatrix(M, cols=c(2,1,3,4), zlim=c(-1,2))

