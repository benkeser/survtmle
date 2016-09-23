#' grad 
#' 
#' A function that computes the gradient of the for a logistic regression
#' model. Used by 
#' \code{optim} on occasion.
#' 
#' @param beta A vector of coefficients in a logistic glm
#' @param X The design matrix 
#' @param Y The outcome
#' 
#' @return Numeric vector of the gradient of the parameter vector

grad <- function(beta, Y, X){
  pi <- plogis(X%*%beta)
  pi[pi==0] <- .Machine$double.neg.eps
  pi[pi==1] <- 1-.Machine$double.neg.eps
  gr <- crossprod(X, Y-pi)
  return(-gr)
}
