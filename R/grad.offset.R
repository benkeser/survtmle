#' grad_offset
#' 
#' A function that computes the gradient of the for a logistic regression
#' model with an offset term. Used by 
#' \code{optim} on occasion.
#' 
#' @param beta A vector of coefficients in a logistic glm
#' @param H The covariate matrix
#' @param offset The offset vector
#' @param Y The outcome
#' @importFrom stats plogis
#' @return Numeric vector of the gradient of the parameter vector


grad_offset <- function(beta, Y, H, offset=NULL){
  pi <- stats::plogis(cbind(offset,H)%*%c(1,beta))
  pi[pi==0] <- .Machine$double.neg.eps
  pi[pi==1] <- 1-.Machine$double.neg.eps
  gr <- crossprod(H, Y-pi)
  return(-gr)
}
