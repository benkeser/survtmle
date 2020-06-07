#' Log-Likelihood
#'
#' @description Computes the log-likelihood for a model. Used by
#'  \code{\link[stats]{optim}} on occasion.
#'
#' @param beta A \code{numeric} vector of coefficients in a logistic GLM.
#' @param X The design matrix.
#' @param Y The outcome.
#' @param wts A \code{numeric} vector of observation-level weights.
#'
#' @importFrom stats plogis
#'
#' @return Numeric of the summed negative log-likelihood loss over
#'  observations.
LogLikelihood <- function(beta, X, Y, wts = rep(1, length(Y))) {
  pi <- stats::plogis(X %*% beta)
  pi[pi == 0] <- .Machine$double.neg.eps
  pi[pi == 1] <- 1 - .Machine$double.neg.eps
  log_like <- sum(wts * (Y * log(pi) + (1 - Y) * log(1 - pi)))
  return(-log_like)
}

#' Log-Likelihood Offset
#'
#' @description Computes the log-likelihood for a logistic regression model
#'  with an offset. Used by \code{\link[stats]{optim}} on occasion.
#'
#' @param beta A vector of coefficients in a logistic GLM.
#' @param Y A \code{numeric} vector of the outcome.
#' @param H The \code{matrix} of covariates.
#' @param offset The \code{numeric} vector of offsets.
#' @param wts A \code{numeric} vector of observation-level weights.
#'
#' @importFrom stats plogis
#'
#' @return Numeric of the summed negative log-likelihood loss over
#'  observations.
LogLikelihood_offset <- function(beta, Y, H, offset,
                                 wts = rep(1, length(Y))) {
  X <- as.matrix(cbind(offset, H))
  pi <- stats::plogis(X %*% as.matrix(c(1, beta)))
  pi[pi == 0] <- .Machine$double.neg.eps
  pi[pi == 1] <- 1 - .Machine$double.neg.eps
  log_like <- sum(wts * (Y * log(pi) + (1 - Y) * log(1 - pi)))
  return(-log_like)
}
