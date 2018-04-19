utils::globalVariables(c("."))

#' confint.survtmle
#'
#' Computes confidence intervals for a fitted \code{survtmle} object.
#'
#' @param object An object of class \code{survtmle}.
#' @param parm A numeric vector indicating which indexes of \code{object$est} to
#'  return confidence intervals for (default is to return all).
#' @param level The confidence level requested.
#' @param ... Other arguments. Not currently used.
#'
#' @return A matrix with columns giving the lower and upper confidence limits
#'  for each parameter. These will be labeled as (1-level)/2 and 1 - (1-level)/2
#'  in percent. The default is 2.5% and 97.5%.
#'
#' @export
#'
#' @method confint survtmle
#'
#' @importFrom stats qnorm
#'
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' ftime <- round(runif(n, 1, 4))
#' ftype <- round(runif(n, 0, 2))
#' trt <- rbinom(n, 1, 0.5)
#' adjustVars <- data.frame(W1 = rnorm(n), W2 = rnorm(n))
#'
#' # fit a survtmle object
#' fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
#'                 adjustVars = adjustVars, glm.trt = "W1 + W2",
#'                 glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2",
#'                 method = "mean", t0 = 4)
#' # get confidence intervals
#' ci <- confint(fit)
#' ci
#
confint.survtmle <- function(object,
                             parm = seq_along(object$est),
                             level = 0.95,
                             ...) {
  a <- (1 - level) / 2
  fac <- stats::qnorm(c(a, 1 - a))
  pct <- format.perc(c(a, 1 - a), 3)
  estVec <- object$est[parm]
  ses <- sqrt(diag(object$var)[parm])
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ci[] <- estVec + outer(ses, fac)
  row.names(ci) <- row.names(object$est)[parm]
  return(ci)
}

################################################################################

#' confint.tp.survtmle
#'
#' Computes confidence intervals for a fitted \code{tp.survtmle} object.
#'
#' @param object An object of class \code{tp.survtmle}, as produced by invoking
#'  the function \code{timepoints} on an object produced by \code{survtmle}, for
#'  which a confidence interval is to be computed.
#' @param parm A numeric vector indicating which indexes of \code{object$est} to
#'  return confidence intervals for (default is to return all). NOT USED NOW.
#' @param level A \code{numeric} indicating the level of the confidence interval
#'  to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @return A list of matrices, each with columns giving the lower and upper
#'  confidence limits for each parameter. These will be labeled as (1-level)/2
#'  and 1 - (1-level)/2 in percent. The default is 2.5% and 97.5%. The list
#'  contains as many matrices as their are comparison groups in the input data.
#'
#' @importFrom stats qnorm
#' @importFrom dplyr "%>%"
#'
#' @method confint tp.survtmle
#'
#' @export
#'
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' ftime <- round(runif(n, 1, 4))
#' ftype <- round(runif(n, 0, 2))
#' trt <- rbinom(n, 1, 0.5)
#' adjustVars <- data.frame(W1 = rnorm(n), W2 = rnorm(n))
#'
#' # fit a survtmle object
#' fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
#'                 adjustVars = adjustVars, glm.trt = "W1 + W2",
#'                 glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2",
#'                 method = "mean", t0 = 4)
#' # extract cumulative incidence at each timepoint
#' tpfit <- timepoints(fit, times = seq_len(4))
#' # get confidence intervals
#' ci <- confint(tpfit)
#' ci
#
confint.tp.survtmle <- function(object,
                                parm,
                                level = 0.95,
                                ...) {
  # first, let's get Z <- (1 - alpha)
  a <- (1 - level) / 2
  fac <- c(-1, 1) * abs(stats::qnorm(p = a))

  # extract vectors of point estimates and variance across all timepoints
  est_allt <- list()
  ses_allt <- list()
  for (i in seq_along(object)) {
    est_allt[[i]] <- object[[i]]$est
    ses_allt[[i]] <- sqrt(diag(object[[i]]$var))
  }

  # find number of contrast groups
  n_grps <- est_allt %>%
    lapply(., nrow) %>%
    unique() %>%
    unlist()

  # construct output object and fill across groups
  ci_tables <- list()
  for (i in seq_len(n_grps)) {
    est <- lapply(est_allt, `[[`, i) %>%
      unlist()
    ses <- lapply(ses_allt, `[[`, i) %>%
      unlist() %>%
      outer(., fac)
    ci <- est + ses
    colnames(ci) <- format.perc(c(a, 1 - a), 3)
    ci_tables[[i]] <- ci
  }
  names(ci_tables) <- lapply(est_allt, row.names) %>%
    unlist() %>%
    unique()
  return(ci_tables)
}

################################################################################

#' format.perc
#'
#' Copied from package \code{stats}.
#'
#' @param probs Probabilities
#' @param digits Number of digits to round to
#
format.perc <- function(probs, digits) {
  paste(format(100 * probs,
    trim = TRUE,
    scientific = FALSE, digits = digits
  ), "%")
}
