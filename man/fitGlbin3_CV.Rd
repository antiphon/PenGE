% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl-bin-fit3.R
\name{fitGlbin3_CV}
\alias{fitGlbin3_CV}
\title{Group Lasso Penalised Logistic Regression Estimation (v3)}
\usage{
fitGlbin3_CV(Q, dfmax = ncol(Q$X) * 1.1, lambda = NULL,
  lambda.min = 0.001, nlambda = 100, log_linear_lambda = FALSE,
  verb = 0, add_intercept = TRUE, quads = NULL, mc.cores = 3,
  sparse = FALSE, orthonormalise = TRUE,
  standardize = !orthonormalise, use_glm_in_lmax = TRUE,
  maxiter = 10000, ...)
}
\arguments{
\item{Q}{Model object}

\item{dfmax}{stop if dfmax significant terms reached}

\item{lambda}{penalty vector. if missing, will be computed}

\item{lambda.min}{If lambda is to be computed, min. lambda will be max.lambda * lambda.min}

\item{nlambda}{steps to generate for lambda vector}

\item{log_linear_lambda}{log-linear or linear}

\item{verb}{Verbosity}

\item{add_intercept}{Add overall intercept?}

\item{quads}{The cross validation quadrats}

\item{mc.cores}{mclapply mc.cores}

\item{sparse}{use Sparse Matrix}

\item{orthonormalise}{Orthonormalise X before computations?}

\item{standardize}{Standardize X before computations? (overriden by orthonormalise)}

\item{use_glm_in_lmax}{use R-s own glm function (non-sparse) to determine max. lambda}

\item{...}{passed on to}
}
\description{
Uses the version 3 of glbinc with orthonormalisation
}
\details{
Basically just a wrapper for package glbinc for estimating group penalised binomial regression.
If 'quads' (list of sub-windows, see 'split_window'-function) are given,
fits the CV-models as well. Here the mc.cores -parameter is used.
}
