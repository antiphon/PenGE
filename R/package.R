#' Penalised Gibbs Estimation
#'
#' @details
#' Estimation of interactions in spatial point patterns using intra- and inter-type interacting Gibbs models.
#'
#' Model inference is based on pseudo-likelihood approximations.
#'
#' Significance of interactions, and therefore the presence of them, is determined by penalised maximum pseudo-likelihood.
#'
#' Default penalisation is Group Lasso. Grouped penalisation allows us to determine if the whole multi-scale potentials, which are given by step-functions/basis-functions, are 0 or not.
#'
#' @aliases PenGE
#' @docType package
"_PACKAGE"
