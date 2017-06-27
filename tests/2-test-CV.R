# Test CV
devtools::load_all(".")

library(spatstat)
data(lansing)

x <- lansing
if(1){
  p <- 6 # we have six species
  rv <- c(0.05, 0.1) # use this for everything
  r1 <- rep( list(rv) , p) # intra-type interaction ranges
  r2 <- rep( list(rv), p * (p-1)/2) # inter-type
  Q <- make_Q_stepper_multi(x=x, ranges1 <-r1 , ranges2 <- r2, verb = 1)

  # quads


  # fit
  fit <-fitGlbin_CV(Q)
}

M <- imatrix(fit, signed = TRUE)
image.imatrix(M, col = c("red","black","green"))

