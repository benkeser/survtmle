#' confint.survtmle
#' 
#' Computes confidence intervals for a fitted \code{survtmle} object.
#' 
#' @param object An object of class \code{survtmle}.
#' @param parm A numeric vector indicating which indices of \code{object$est} to
#' return confidence intervals for (default is to return all).
#' @param level The confidence level requested. 
#' 
#' @return A matrix with columns giving the lower and upper confidence limits
#' for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % 
#' (by default 2.5% and 97.5%).
#' 
#' @importFrom stats format.perc
#' @export
#' 
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' ftime <- round(runif(n,1,4))
#' ftype <- round(runif(n,0,2))
#' trt <- rbinom(n,1,0.5) + rbinom(n,1,0.25)
#' adjustVars <- data.frame(W1 = rnorm(n), W2 = rnorm(n))
#' 
#' # Fit 1
#' # fit a survtmle object 
#' fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", t0=6)
#' # get confidence intervals
#' ci <- confint(fit1)
#' # ci

confint.survtmle <- function(object, parm = 1:length(object$est), level = 0.95){
	estVec <- object$est[parm]
	ses <- sqrt(diag(object$var)[parm])
	a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    pct <- stats:::format.perc(a, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,pct))
    ci[] <- estVec + ses %o% fac
    row.names(ci) <- row.names(object$est)[parm]
    ci
}