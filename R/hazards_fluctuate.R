#' Fluctuation for the Method of Cause-Specific Hazards
#'
#' @description This function performs a fluctuation of an initial estimate of
#'  the cause-specific hazard functions via a call to \code{\link[stats]{glm}}
#'  (i.e., a logistic submodel) or a call to \code{\link[stats]{optim}} (to
#'  ensure fluctuations stay within model space). The structure of the function
#'  is specific to how it is called within \code{\link{hazard_tmle}}. In
#'  particular, \code{dataList} must have a very specific structure for this
#'  function to run properly. The list should consist of \code{data.frame}
#'  objects. The first will have the number of rows for each observation equal
#'  to the \code{ftime} corresponding to that observation. Subsequent entries
#'  will have \code{t0} rows for each observation and will set \code{trt}
#'  column equal to each value of \code{trtOfInterest} in turn. The function
#'  will fit a logistic regression with (a scaled version of) \code{Nj} as
#'  outcome, the logit of the current (pseudo-) hazard estimate as offset and
#'  the targeted minimum loss-based estimation "clever covariates". The
#'  function then obtains predictions based on this fit on each of the
#'  \code{data.frame} objects in \code{dataList}.
#'
#' @param dataList A list of \code{data.frame} objects.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param ofInterestJ Numeric vector indicating \code{ftypeOfInterest} that was
#'  passed to \code{\link{hazard_tmle}}.
#' @param nJ The number of unique failure types.
#' @param uniqtrt The values of \code{trtOfInterest} passed to
#'  \code{\link{mean_tmle}}.
#' @param ntrt The number of \code{trt} values of interest.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#' @param verbose A \code{logical} indicating whether the function should print
#'  messages to indicate progress.
#' @param att A \code{boolean} indicating whether to compute the ATT estimate,
#'  instead of treatment specific survival curves. This option only works with 
#'  two levels of \code{trt} that are labeled with 0 and 1.
#' @param ... Other arguments. Not currently used.
#'
#' @return The function returns a list that is exactly the same as the input
#'  \code{dataList}, but with updated columns corresponding with estimated
#'  cumulative incidence at each time and estimated "clever covariates" at each
#'  time.
#'
#' @importFrom Matrix Diagonal
#' @importFrom stats optim
fluctuateHazards <- function(dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt,
                             t0, verbose, att, ...) {
  eps <- NULL
  for (z in uniqtrt) {
    for (j in allJ) {
      # clever covariates
      cleverCovariatesNotSelf <- NULL
      if (length(ofInterestJ[ofInterestJ != j]) > 0) {
        cleverCovariatesNotSelf <- c(
          cleverCovariatesNotSelf,
          paste0(
            "H", ofInterestJ[ofInterestJ != j],
            ".jNotSelf.z", z
          )
        )
      }
      if (j %in% ofInterestJ) {
        cleverCovariatesSelf <- paste0("H", j, ".jSelf.z", z)
      } else {
        cleverCovariatesSelf <- NULL
      }

      # calculate offset term and outcome
      dataList <- lapply(dataList, function(x, j, allJ) {
        x$thisScale <- pmin(x[[paste0("u", j)]], 1 - x[[paste0("hazNot", j)]])
        -x[[paste0("l", j)]]
        x$thisOffset <- stats::qlogis(pmin(
          (x[[paste0("Q", j, "Haz")]] - x[[paste0("l", j)]]) / x$thisScale,
          1 - .Machine$double.neg.eps
        ))
        x$thisOutcome <- (x[[paste0("N", j)]] - x[[paste0("l", j)]]) /
          x$thisScale
        x
      }, j = j, allJ = allJ)

      fluc.mod <- stats::optim(
        par = rep(0, length(c(
          cleverCovariatesNotSelf,
          cleverCovariatesSelf
        ))),
        fn = LogLikelihood_offset,
        Y = dataList[[1]]$thisOutcome,
        H = suppressWarnings(
          as.matrix(Matrix::Diagonal(x = dataList[[1]]$thisScale) %*%
            as.matrix(dataList[[1]][, c(
              cleverCovariatesNotSelf,
              cleverCovariatesSelf
            )]))
        ),
        offset = dataList[[1]]$thisOffset,
        method = "BFGS", gr = grad_offset,
        control = list(reltol = 1e-7, maxit = 50000)
      )

      if (fluc.mod$convergence != 0) {
        warning("Fluctuation convergence failure. Using with initial estimates.
              Proceed with caution")
        beta <- rep(0, length(fluc.mod$par))
      } else {
        beta <- fluc.mod$par
      }
      eps <- c(eps, beta)

      dataList <- lapply(dataList, function(x, j) {
        x[[paste0("Q", j, "PseudoHaz")]][x$trt == z] <-
          plogis(x$thisOffset[x$trt == z] +
            suppressWarnings(
              as.matrix(
                Matrix::Diagonal(x = x$thisScale[x$trt == z]) %*%
                  as.matrix(x[x$trt == z, c(
                    cleverCovariatesNotSelf,
                    cleverCovariatesSelf
                  )])
              ) %*% as.matrix(beta)
            ))
        x[[paste0("Q", j, "Haz")]][x$trt == z] <-
          x[[paste0("Q", j, "PseudoHaz")]][x$trt == z] *
          x$thisScale[x$trt == z] + x[[paste0("l", j)]][x$trt == z]
        x
      }, j = j)

      # update variables based on new haz
      dataList <- updateVariables(
        dataList = dataList, allJ = allJ,
        ofInterestJ = ofInterestJ,
        nJ = nJ, uniqtrt = uniqtrt, ntrt = ntrt,
        verbose = verbose, t0 = t0, att = att
      )
    }
  }
  attr(dataList, "fluc") <- eps
  dataList
}
