#' LogLikelihood_offset
#' 
#' A function that computes the log-likelihood for a logistic regression
#' model with an offset. Used by 
#' \code{optim} on occasion.
#' 
#' @importFrom stats plogis
#' @param beta A vector of coefficients in a logistic glm
#' @param Y The outcome vector
#' @param H The covariate matrix
#' @param offset The offset vector
#' @return Numeric of the summed negative log-likelihood loss over observations. 
#' 


LogLikelihood_offset <- function(beta, Y, H, offset){
  X <- as.matrix(cbind(offset, H))
  pi <- stats::plogis(X%*%as.matrix(c(1,beta)))
  pi[pi==0] <- .Machine$double.neg.eps
  pi[pi==1] <- 1-.Machine$double.neg.eps
  logLike <- sum(Y*log(pi) + (1-Y)*log(1-pi))
  return(-logLike)    
}
