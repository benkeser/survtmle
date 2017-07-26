#' Fluctuation for the Method of Iterated Means
#'
#' This function performs a fluctuation of an initial estimate of the
#' G-computation regression at a specified time \code{t} using a call to
#' \code{glm} (i.e., a logistic submodel) or a call to \code{optim} (if bounds
#' are specified). The structure of the function is specific to how it is called
#' within \code{mean_tmle}. In particular, \code{wideDataList} must have a very
#' specific structure for this function to run properly. The list should consist
#' of \code{data.frame} objects. The first should have all rows set to their
#' observed value of \code{trt}. The remaining should in turn have all rows set
#' to each value of \code{trtOfInterest} in the \code{survtmle} call. The latter
#' will be used to obtain predictions that are then mapped into the estimates of
#' the cumulative incidence function at \code{t0}. Currently the code requires
#' each \code{data.frame} to have named columns for each name in
#' \code{names(adjustVars)}, as well as a column named \code{trt}. It must also
#' have a columns named \code{Nj.Y} where j corresponds with the numeric values
#' input in \code{allJ}. These are the indicators of failure due to the various
#' causes before time \code{t} and are necessary for determining who to include
#' in the fluctuation regression. Similarly, each \code{data.frame} should have
#' a column call \code{C.Y} where Y is again \code{t-1}, so that right censored
#' observations are not included in the regressions. The function will fit a
#' logistic regression with \code{Qj.star.t + 1} as outcome (also needed as a
#' column in \code{wideDataList}) with offset \code{qlogis(Qj.star.t)} and
#' number of additional covariates given by \code{length(trtOfInterest)}. These
#' additional covariates should be columns in the each \code{data.frame} in
#' \code{wideDataList} called \code{H.z.t} where \code{z} corresponds to a each
#' unique value of \code{trtOfInterest}. The function returns the same
#' \code{wideDataList}, but with a column called \code{Qj.star.t} added to it,
#' which is the fluctuated initial regression estimate evaluated at the observed
#' data points.
#'
#' @param wideDataList A list of \code{data.frame} objects.
#' @param t The timepoint at which to compute the iterated mean.
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean_tmle}.
#' @param whichJ Numeric value indicating the cause of failure for which
#'        regression should be computed.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#'        Needed only because the naming convention for the regression if
#'        \code{t == t0} is different than if \code{t != t0}.
#' @param bounds A list of bounds to be used when performing the outcome
#'        regression (Q) with the Super Learner algorithm. NOT YET IMPLEMENTED.
#' @param Gcomp A boolean indicating whether \code{mean_tmle} was called to
#'        evaluate the G-computation estimator, in which case this function does
#'        nothing but re-label columns.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats as.formula optim glm qlogis
#' @importFrom Matrix Diagonal
#'
#' @return The function then returns a list that is exactly the same as the
#'         input \code{wideDataList}, but with a column named \code{Qj.star.t}
#'         added to it, which is the fluctuated conditional mean of
#'         \code{Qj.star.t+1} evaluated at the each of the rows of each
#'         \code{data.frame} in \code{wideDataList}.
#'

fluctuateIteratedMean <- function(wideDataList, t, uniqtrt, whichJ, allJ, t0,
                                  Gcomp = FALSE, bounds = NULL, ...) {
  outcomeName <- ifelse(t == t0, paste("N", whichJ, ".", t0, sep = ""),
                        paste("Q", whichJ, "star.", t + 1, sep = ""))

  ## determine who to include in estimation
  include <- rep(TRUE, nrow(wideDataList[[1]]))
  if(t != 1) {
    for(j in allJ) {
      # exclude previously failed subjects
      include[wideDataList[[1]][[paste0("N",j,".",t-1)]] == 1] <- FALSE
    }
    # exclude previously censored subjects
    include[wideDataList[[1]][[paste0("C.",t-1)]]==1] <- FALSE
  }
  if(is.null(bounds)) {
    wideDataList <- lapply(wideDataList, function(x, t) {
      # check for 0's and 1's
      x[[paste0("Q", whichJ, ".", t)]][x[[paste0("Q", whichJ, ".", t)]] <
                                       .Machine$double.neg.eps] <- .Machine$double.neg.eps
      x[[paste0("Q", whichJ, ".", t)]][x[[paste0("Q", whichJ, ".", t)]] >
                                       1 - .Machine$double.neg.eps] <- 1 - .Machine$double.neg.eps
      x
    }, t = t)
 
    flucForm <- paste(outcomeName, "~ -1 + offset(stats::qlogis(Q", whichJ,
                      ".", t, ")) +",
                      paste0("H", uniqtrt, ".", t, collapse = "+"), sep = "")

    if(!Gcomp) {
      # fluctuation model
      suppressWarnings(
        flucMod <- stats::glm(stats::as.formula(flucForm), family = "binomial",
                              data = wideDataList[[1]][include,],
                              start = rep(0, length(uniqtrt)))
      )
      # get predictions back
      wideDataList <- lapply(wideDataList, function(x, t) {
        suppressWarnings(
          x[[paste0("Q", whichJ, "star.", t)]] <- predict(flucMod, newdata = x,
                                                          type = "response")
        )
        x
      }, t = t)
    } else {
      # if Gcomp, just skip fluctuation step and assign this
      wideDataList <- lapply(wideDataList, function(x, t) {
        x[[paste0("Q", whichJ, "star.", t)]] <- x[[paste0("Q", whichJ, ".", t)]]
        x
      }, t = t)
    }
  } else {
    if(!Gcomp) {
      cleverCovariates <- paste0("H", uniqtrt, ".", t)
      lj.t <- paste0("l",whichJ,".",t)
      uj.t <- paste0("u",whichJ,".",t)
      Qtildej.t <- paste0("Qtilde",whichJ,".",t)
      Nj.tm1 <- paste0("N",whichJ,".",t-1)
      Qj.t <- paste0("Q",whichJ,".",t)
      NnotJ.tm1 <- paste0("NnotJ.",t-1)
      # calculate offset term and outcome
      wideDataList <- lapply(wideDataList, function(x) {
        x[["thisOutcome"]] <- (x[[outcomeName]] - x[[lj.t]])/(x[[uj.t]]-x[[lj.t]])
        x[["thisScale"]] <- x[[uj.t]] - x[[lj.t]]
        x[[Qtildej.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
          (x[[Qj.t]] - x[[lj.t]])/x[["thisScale"]]
        x[[Qtildej.t]][x[[Qtildej.t]]==1] <- 1-.Machine$double.neg.eps

        x$thisOffset <- 0
        x$thisOffset[(x[[NnotJ.tm1]] + x[[Nj.tm1]])==0] <- 
          stats::qlogis(x[[Qtildej.t]][(x[[NnotJ.tm1]] + x[[Nj.tm1]])==0])
        x
      })

      fluc.mod <- stats::optim(par = rep(0, length(cleverCovariates)),
                               fn = LogLikelihood_offset,
                               Y = wideDataList[[1]]$thisOutcome[include],
                               H = as.matrix(wideDataList[[1]][include,
                                             cleverCovariates]),
                               offset = wideDataList[[1]]$thisOffset[include],
                               method = "BFGS", gr = grad_offset,
                               control = list(reltol = 1e-7, maxit = 50000))
 
      if(fluc.mod$convergence != 0) {
        stop("fluctuation convergence failure")
      } else {
        beta <- fluc.mod$par

        wideDataList <- lapply(wideDataList, function(x) {
          x[[paste0("Q",whichJ,"star.",t)]] <- x[[Nj.tm1]] + 
            (1 - x[[NnotJ.tm1]] - x[[Nj.tm1]])* (plogis(x$thisOffset + 
              as.matrix(x[, cleverCovariates]) %*% as.matrix(beta))*x$thisScale
            + x[[lj.t]])
          x
        })
      }
    } else {
      wideDataList <- lapply(wideDataList, function(x, t) {
          x[[paste0("Q",whichJ,"star.",t)]] <- x[[Qj.t]]
          x
        }, t = t)
    }
  }
  wideDataList
}
