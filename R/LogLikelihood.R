#' LogLikelihood 
#' 
#' A function that computes the log-likelihood for a model. Used by 
#' \code{optim} on occasion.
#' 
#' @param beta A vector of coefficients in a logistic glm
#' @param X The design matrix 
#' @param Y The outcome
#' @importFrom stats plogis
#' @return Numeric of the summed negative log-likelihood loss over observations. 
#' 

LogLikelihood <- function(beta, X, Y){
  pi <- stats::plogis(X%*%beta)
  pi[pi==0] <- .Machine$double.neg.eps
  pi[pi==1] <- 1-.Machine$double.neg.eps
  logLike <- sum(Y*log(pi) + (1-Y)*log(1-pi))
  return(-logLike)
}
