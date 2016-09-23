#' print.survtmle
#' 
#' The print method for an object of class \code{survtmle}
#' 
#' @param x An object of class \code{survtmle}
#' 
#' @return Prints the estimates of cumulative incidence and the diagonal
#' of the estimated covariance matrix. 
#' 
#' @export

print.survtmle <- function(x,...){
  if(length(x$est) <= 4){
    print(x[c("est","var")])
  }else{
    out <- list(est=x$est, var=diag(x$var))
    print(out)
  }
}