#' LogLikelihood.offset
#' 
#' A function that computes the log-likelihood for a logistic regression
#' model with an offset. Used by 
#' \code{optim} on occasion.
#' 
#' @param beta A vector of coefficients in a logistic glm
#' @param X The design matrix 
#' @param Y The outcome
#' 
#' @return Numeric of the summed negative log-likelihood loss over observations. 
#' 


LogLikelihood.offset <- function(beta, Y, H, offset){
  X <- as.matrix(cbind(offset, H))
  pi <- plogis(X%*%as.matrix(c(1,beta)))
  pi[pi==0] <- .Machine$double.neg.eps
  pi[pi==1] <- 1-.Machine$double.neg.eps
  logLike <- sum(Y*log(pi) + (1-Y)*log(1-pi))
  return(-logLike)    
}
