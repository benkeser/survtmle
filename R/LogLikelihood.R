#' Log-Likelihood
#'
#' Computes the log-likelihood for a model. Used by \code{optim} on occasion.
#'
#' @param beta A vector of coefficients in a logistic GLM.
#' @param X The design matrix.
#' @param Y The outcome.
#'
#' @importFrom stats plogis
#' @return Numeric of the summed negative log-likelihood loss over observations.
#'

LogLikelihood <- function(beta, X, Y){
  pi <- stats::plogis(X %*% beta)
  pi[pi == 0] <- .Machine$double.neg.eps
  pi[pi == 1] <- 1 - .Machine$double.neg.eps
  log_like <- sum(Y * log(pi) + (1 - Y) * log(1 - pi))
  return(-log_like)
}

#' Log-Likelihood Offset
#'
#' Computes the log-likelihood for a logistic regression model with an offset.
#' Used by \code{optim} on occasion.
#'
#' @param beta A vector of coefficients in a logistic GLM.
#' @param Y A vector of the outcome.
#' @param H The matrix of covariates.
#' @param offset The vector of offsets.
#'
#' @importFrom stats plogis
#'
#' @return Numeric of the summed negative log-likelihood loss over observations.
#'

LogLikelihood_offset <- function(beta, Y, H, offset){
  X <- as.matrix(cbind(offset, H))
  pi <- stats::plogis(X %*% as.matrix(c(1, beta)))
  pi[pi == 0] <- .Machine$double.neg.eps
  pi[pi == 1] <- 1 - .Machine$double.neg.eps
  log_like <- sum(Y * log(pi) + (1 - Y) * log(1 - pi))
  return(-log_like)
}
