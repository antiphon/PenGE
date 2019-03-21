# Dev genstepper fitting

devtools::load_all()

library(spatstat)
library(parallel)




# univariate
if(0) {
  Rmax <- 0.3
  knots <- seq(0, Rmax, l = 6)
  rmat <- rbind(knots[-6],knots[-1])
  #x <- spatstat::setmarks(spatstat:::rpoint(1000), factor(1))
  x <- spatstat::setmarks(spatstat:::rThomas(50, 0.02, 5), factor(1))

  Q <- make_Q_stepper_multi_rmat(x, ranges1 = list(rmat), ranges2=NULL, sat1=20, auto_sat=F)

  z <- fitGlbin3_CV(Q)

  # Check fit

  r <- seq(0, Rmax, l = 100)

  e <- eval_estimated_potential(z, k =  100, at = r)

  plot(r, e, ylim = c(-1,1) * .5, "l")
}


# Bivariate independent stuff
#x <- spatstat::setmarks(spatstat:::rpoint(1000), factor(1))
if(0){
  #x <- spatstat::rmpoint(c(10,10)*100)
  x <- spatstat::superimpose("1"=spatstat::rpoint(100), "2"=spatstat::rThomas(15,0.04,10))

  #
  Rmax <- 0.3
  knots <- seq(0, Rmax, l = 4)
  rmat <- rbind(knots[-4],knots[-1])

  r1 <- c(0, Rmax/4, Rmax/2)
  rmat <- rbind(r1, r1+Rmax/2)

  knots1 <- rep(list(rmat), 2)
  knots2 <- list(rmat)

  Q <- make_Q_stepper_multi_rmat(x, ranges1 = knots1, ranges2=knots2, sat1=200, sat2=2000,
                            auto_sat=T, border_r = 0.2, verb=1)

  #
  z <- fitGlbin3_CV(Q, verb=0, maxiter=5000, penalty = "grSCAD")#, lambda.min = 0.001)

  # Check fit
  k <- 0.5 * which.min(z$fullfit$aic)
  r <- seq(0, Rmax, l = 100)

  par(mfrow=c(2,2))
  plot(z$fullfit$aic)
  abline(v=k)
  for(i in 1:2) for(j in i:2){
    e <- eval_estimated_potential(z, k =  k, at = r, i = i, j = j)
    plot(r, e, ylim = c(-1,1) * .5, main=paste(i,j, sep="-"), type="l")
  }

}

