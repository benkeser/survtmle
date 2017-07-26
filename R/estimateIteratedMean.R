#' Estimation for the Method of Iterated Means
#'
#' This function computes an estimate of the G-computation regression at a
#' specified time \code{t} using either \code{glm} or \code{SuperLearner}. The
#' structure of the function is specific to how it is called within
#' \code{mean_tmle}. In particular, \code{wideDataList} must have a very
#' specific structure for this function to run properly. The list should consist
#' of \code{data.frame} objects. The first should have all rows set to their
#' observed value of \code{trt}. The remaining should in turn have all rows set
#' to each value of \code{trtOfInterest} in the \code{survtmle} call. Currently
#' the code requires each \code{data.frame} to have named columns for each name
#' in \code{names(adjustVars)}, as well as a column named \code{trt}. It must
#' also have a columns named \code{Nj.Y} where j corresponds with the numeric
#' values input in \code{allJ}. These are the indicators of failure due to the
#' various causes before time \code{t} and are necessary for determining who to
#' include in the regression. Similarly, each \code{data.frame} should have a
#' column call \code{C.Y} where Y is again \code{t - 1}, so that right censored
#' observations are not included in the regressions. The function will fit a
#' regression with \code{Qj.star.t+1} (also needed as a column in
#' \code{wideDataList}) on functions of \code{trt} and \code{names(adjustVars)}
#' as specified by \code{glm.ftime} or \code{SL.ftime}.
#'
#' @param wideDataList A list of \code{data.frame} objects.
#' @param t The timepoint at which to compute the iterated mean.
#' @param whichJ Numeric value indicating the cause of failure for which
#'        regression should be computed.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#'        Needed only because the naming convention for the regression if
#'        \code{t == t0} is different than if \code{t != t0}.
#' @param adjustVars Object of class \code{data.frame} that contains the
#'        variables to adjust for in the regression.
#' @param SL.ftime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        outcome regression (either cause-specific hazards or conditional mean).
#'        See \code{?SuperLearner} for more information on how to specify valid
#'        \code{SuperLearner} libraries. It is expected that the wrappers used
#'        in the library will play nicely with the input variables, which will
#'        be called \code{"trt"} and \code{names(adjustVars)}.
#' @param glm.ftime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the outcome regression (either cause-specific hazards or
#'        conditional mean). Ignored if \code{SL.ftime != NULL}. Use \code{"trt"}
#'        to specify the treatment in this formula (see examples). The formula
#'        can additionally include any variables found in
#'        \code{names(adjustVars)}.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param returnModels A boolean indicating whether to return the
#'        \code{SuperLearner} or \code{glm} objects used to estimate the
#'        nuisance parameters. Must be set to \code{TRUE} if the user plans to
#'        use calls to \code{timepoints} to obtain estimates at times other than
#'        \code{t0}. See \code{?timepoints} for more information.
#' @param bounds A list of bounds to be used when performing the outcome
#'        regression (Q) with the Super Learner algorithm. NOT YET IMPLEMENTED.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats as.formula predict model.matrix optim glm
#' @importFrom SuperLearner SuperLearner SuperLearner.CV.control
#'
#' @return The function then returns a list that is exactly the same as the
#'         input \code{wideDataList}, but with a column named \code{Qj.t} added
#'         to it, which is the estimated conditional mean of \code{Qj.star.t+1}
#'         evaluated at the each of the rows of each \code{data.frame} in
#'         \code{wideDataList}.
#'

estimateIteratedMean <- function(wideDataList, t, whichJ, allJ, t0, adjustVars,
                                 SL.ftime = NULL, glm.ftime = NULL, verbose,
                                 returnModels = FALSE, bounds = NULL, ...){
  ## determine who to include in estimation
  include <- rep(TRUE, nrow(wideDataList[[1]]))
  if(t != 1) {
    for(j in allJ) {
      # exclude previously failed subjects
      include[wideDataList[[1]][[paste0("N",j,".",t-1)]]==1] <- FALSE
    }
    # exclude previously censored subjects
    include[wideDataList[[1]][[paste0("C.",t-1)]]==1] <- FALSE
  }

  ## determine the outcome for the regression
  outcomeName <- ifelse(t == t0, paste("N", whichJ, ".", t0, sep = ""),
                        paste("Q", whichJ, "star.", t + 1, sep = ""))

  ## create an indicator of any failure prior to t
  wideDataList <- lapply(wideDataList, function(x, t){
    if(length(allJ) > 1) {
      x[[paste0("NnotJ.",t-1)]] <- 
        rowSums(cbind(rep(0, nrow(x)), x[, paste0('N', allJ[allJ != whichJ], '.', t - 1)]))
    } else {
      x[[paste0("NnotJ.",t-1)]] <- 0
    }
    x
  },t = t)

  lj.t <- paste0("l",whichJ,".",t)
  uj.t <- paste0("u",whichJ,".",t)
  Qtildej.t <- paste0("Qtilde",whichJ,".",t)
  Nj.tm1 <- paste0("N",whichJ,".",t-1)
  Qj.t <- paste0("Q",whichJ,".",t)
  NnotJ.tm1 <- paste0("NnotJ.",t-1)
  ## GLM code
  if(is.null(SL.ftime)) {
    if(is.null(bounds)) { # with no bounds
      Qform <- paste(outcomeName, "~", glm.ftime)
      suppressWarnings({
        Qmod <- stats::glm(as.formula(Qform), family = "binomial",
                           data = wideDataList[[1]][include,])

        wideDataList <- lapply(wideDataList, function(x, whichJ, t) {
          suppressWarnings(
            x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
                predict(Qmod,newdata=x,type="response")
          )
          x
        }, t = t, whichJ = whichJ)
      })
    } else { # with bounds
      Qform <- paste(outcomeName, "~", glm.ftime)
      X <- model.matrix(as.formula(Qform), data = wideDataList[[1]][include, ])
      Ytilde <- (wideDataList[[1]][include,outcomeName] - wideDataList[[1]][[lj.t]][include])/
        (wideDataList[[1]][[uj.t]][include] - wideDataList[[1]][[lj.t]][include])
      Qmod <- optim(par = rep(0, ncol(X)), fn = LogLikelihood, Y = Ytilde,
                    X = X, method = "BFGS", gr = grad,
                    control = list(reltol = 1e-7, maxit = 50000))
      beta <- Qmod$par
      wideDataList <- lapply(wideDataList, function(x, j, t) {
        newX <- model.matrix(as.formula(Qform), data = x)
        x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
          (plogis(newX%*%beta)*(x[[uj.t]]-x[[lj.t]]) + x[[lj.t]])

        x
      }, j = whichJ, t = t)
    }
  } else if(is.null(glm.ftime)) { # Super Learner
    if(is.null(bounds)) { # with no bounds
      # some stability checks
      # number of unique outcome values
      nUniq <- length(unique(wideDataList[[1]][include,outcomeName]))
      cvControl <- SuperLearner::SuperLearner.CV.control()
      if(t == t0) {
        # if there are less than 2 events at t0, just fit regression using only Z
        nE <- sum(wideDataList[[1]][include, outcomeName])
        ignoreSL <- nE <= 2
        if(ignoreSL) {
          suppressWarnings({
            Qmod <- stats::glm(stats::as.formula(paste0(outcomeName, " ~ trt")),
                               data = wideDataList[[1]][include, ])
            wideDataList <- lapply(wideDataList, function(x, whichJ, t) {
              suppressWarnings(
              x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]- x[[Nj.tm1]])*
                predict(Qmod,newdata=data.frame(trt=x$trt))
             )
             x
            }, t = t, whichJ = whichJ)
          })
        } else {
          simplify <- nE <= cvControl$V
          if(simplify) cvControl <- list(V = nE - 1, stratifyCV = TRUE)
          suppressWarnings(
            Qmod <- SuperLearner::SuperLearner(Y = wideDataList[[1]][include,outcomeName],
                                               X = wideDataList[[1]][include, c("trt", names(adjustVars))],
                                               SL.library = SL.ftime,
                                               cvControl = cvControl,
                                               family = "binomial",
                                               verbose = verbose)
          )
          wideDataList <- lapply(wideDataList, function(x, whichJ, t) {
            x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
              predict(Qmod, newdata = x[, c('trt', names(adjustVars))], onlySL = TRUE)$pred
            x
          }, t = t, whichJ = whichJ)
        }
      } else {
        suppressWarnings(
          Qmod <- SuperLearner::SuperLearner(Y = wideDataList[[1]][include, outcomeName],
                                             X=wideDataList[[1]][include, c("trt", names(adjustVars))],
                                             SL.library = SL.ftime,
                                             cvControl = cvControl,
                                             family = "binomial",
                                             verbose = verbose)
        )
        wideDataList <- lapply(wideDataList, function(x, whichJ, t) {
          suppressWarnings(
            x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[Nj.tm1]]-x[[NnotJ.tm1]])*
              predict(Qmod, newdata = x[, c('trt', names(adjustVars))], onlySL = TRUE)$pred
          )
          x
        }, t = t, whichJ = whichJ)
      }
    } else {
      stop("Super Learner code with bounds not written yet")
    }
  }
  out <- list(wideDataList = wideDataList,
              ftimeMod = if(returnModels)
                Qmod
              else
                NULL)
  out
}
