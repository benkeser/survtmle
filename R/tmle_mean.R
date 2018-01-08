#' TMLE for G-Computation of Cumulative Incidence
#'
#' This function estimates the marginal cumulative incidence for failures of
#' specified types using targeted minimum loss-based estimation based on the
#' G-computation representation of cumulative incidence. The function is called
#' by \code{survtmle} whenever \code{method = "mean"} is specified. However,
#' power users could, in theory, make calls directly to this function.
#'
#' @param ftime A numeric vector of failure times. Right-censored observations
#'        should have corresponding \code{ftype} set to 0.
#' @param ftype A numeric vector indicating the type of failure. Observations
#'        with \code{ftype=0} are treated as a right-censored observation. Each
#'        unique value besides zero is treated as a separate type of failure.
#' @param trt A numeric vector indicating observed treatment assignment. Each
#'        unique value will be treated as a different type of treatment.
#'        Currently, only two unique values are supported.
#' @param adjustVars A \code{data.frame} of adjustment variables that will be
#'        used in estimating the conditional treatment, censoring, and failure
#'        (hazard or conditional mean) probabilities.
#' @param t0 The time at which to return cumulative incidence estimates. By
#'        default this is set to \code{max(ftime[ftype > 0])}.
#' @param SL.ftime A character vector or list specification to be passed to the
#'        \code{SL.library} option in the call to \code{SuperLearner} for the
#'        outcome regression (either cause-specific hazards or iterated mean).
#'        See \code{?SuperLearner} for more information on how to specify valid
#'        \code{SuperLearner} libraries. It is expected that the wrappers used
#'        in the library will play nicely with the input variables, which will
#'        be called \code{"trt"}, \code{names(adjustVars)}, and \code{"t"} (if
#'        \code{method = "hazard"}).
#' @param SL.ctime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        estimate of the conditional hazard for censoring. It is expected that
#'        the wrappers used in the library will play nicely with the input
#'        variables, which will be called \code{"trt"} and
#'        \code{names(adjustVars)}.
#' @param SL.trt A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        estimate of the conditional probability of treatment. It is expected
#'        that the wrappers used in the library will play nicely with the input
#'        variables, which will be \code{names(adjustVars)}.
#' @param glm.ftime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the outcome regression. Ignored if \code{SL.ftime} is not equal to
#'        \code{NULL}. Use \code{"trt"} to specify the treatment in this formula
#'        (see examples). The formula can additionally include any variables
#'        found in \code{names(adjustVars)}.
#' @param glm.ctime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the estimate of the conditional hazard for censoring. Ignored if
#'        \code{SL.ctime} is not equal to \code{NULL}. Use \code{"trt"} to
#'        specify the treatment in this formula (see examples). The formula can
#'        additionally include any variables found in \code{names(adjustVars)}.
#' @param glm.trt A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the estimate of the conditional probability of treatment. Ignored
#'        if \code{SL.trt} is not equal to \code{NULL}. The formula can include
#'        any variables found in \code{names(adjustVars)}.
#' @param glm.family The type of regression to be performed if fitting GLMs in
#'        the estimation and fluctuation procedures. The default is "binomial"
#'        for logistic regression. Only change this from the default if there
#'        are justifications that are well understood. This is passed directly
#'        to \code{estimateCensoring}.
#' @param returnIC A boolean indicating whether to return vectors of influence
#'        curve estimates. These are needed for some post-hoc comparisons, so it
#'        is recommended to leave as \code{TRUE} (the default) unless the user
#'        is sure these estimates will not be needed later.
#' @param returnModels A boolean indicating whether to return the
#'        \code{SuperLearner} or \code{glm} objects used to estimate the
#'        nuisance parameters. Must be set to \code{TRUE} if the user plans to
#'        use \code{timepoints} to obtain estimates of incidence at times other
#'        than \code{t0}. See \code{?timepoints} for more information.
#' @param ftypeOfInterest An input specifying what failure types to compute
#'        estimates of incidence for. The default value computes estimates for
#'        values \code{unique(ftype)}. Can alternatively be set to a vector of
#'        values found in \code{ftype}.
#' @param trtOfInterest An input specifying which levels of \code{trt} are of
#'        interest. The default value computes estimates for values
#'        \code{unique(trt)}. Can alternatively be set to a vector of values
#'        found in \code{trt}.
#' @param bounds A \code{data.frame} of bounds on the conditional hazard
#'        function (if \code{method = "hazard"}) or on the iterated conditional
#'        means (if \code{method = "mean"}). The \code{data.frame} should have a
#'        column named \code{"t"} that includes values \code{1:t0}. The other
#'        columns should be names \code{paste0("l",j)} and \code{paste0("u",j)}
#'        for each unique failure type label j, denoting lower and upper bounds,
#'        respectively. See examples.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress. If \code{SuperLearner} is called
#'        internally, this option will additionally be passed to
#'        \code{SuperLearner}.
#' @param Gcomp A boolean indicating whether to compute the G-computation
#'        estimator (i.e., a substitution estimator with no targeting step).
#'        Theory does not support inference for the G-computation estimator if
#'        Super Learner is used to estimate failure and censoring distributions.
#'        The G-computation is only implemented if \code{method = "mean"}.
#' @param gtol The truncation level of predicted censoring survival. Setting to
#'        larger values can help performance in data sets with practical
#'        positivity violations.
#' @param ... Other options. Not currently used.
#'
#' @return An object of class \code{survtmle}.
#' \describe{
#' \item{call}{The call to \code{survtmle}.}
#' \item{est}{A numeric vector of point estimates -- one for each combination of
#'            \code{ftypeOfInterest} and \code{trtOfInterest}.}
#' \item{var}{A covariance matrix for the point estimates.}
#' \item{meanIC}{The empirical mean of the efficient influence function at the
#'               estimated, targeted nuisance parameters. Each value should be
#'               small or the user will be warned that excessive finite-sample
#'               bias may exist in the point estimates.}
#' \item{ic}{The efficient influence function at the estimated, fluctuated
#'           nuisance parameters, evaluated on each of the observations. These
#'           are used to construct confidence intervals for post-hoc
#'           comparisons.}
#' \item{ftimeMod}{If \code{returnModels=TRUE} the fit object(s) for the call to
#'                 \code{glm} or \code{SuperLearner} for the outcome regression
#'                 models. If \code{method="mean"} this will be a list of length
#'                 \code{length(ftypeOfInterest)} each of length \code{t0} (one
#'                 regression for each failure type and for each timepoint). If
#'                 \code{method="hazard"} this will be a list of length
#'                 \code{length(ftypeOfInterest)} with one fit corresponding to
#'                 the hazard for each cause of failure. If
#'                 \code{returnModels = FALSE}, this entry will be \code{NULL}.}
#' \item{ctimeMod}{If \code{returnModels = TRUE} the fit object for the call to
#'                 \code{glm} or \code{SuperLearner} for the pooled hazard
#'                 regression model for the censoring distribution. If
#'                 \code{returnModels = FALSE}, this entry will be \code{NULL}.}
#' \item{trtMod}{If \code{returnModels = TRUE} the fit object for the call to
#'               \code{glm} or \code{SuperLearner} for the conditional
#'               probability of \code{trt} regression model. If
#'               \code{returnModels = FALSE}, this entry will be \code{NULL}.}
#' \item{t0}{The timepoint at which the function was evaluated.}
#' \item{ftime}{The numeric vector of failure times used in the fit.}
#' \item{ftype}{The numeric vector of failure types used in the fit.}
#' \item{trt}{The numeric vector of treatment assignments used in the fit.}
#' \item{adjustVars}{The \code{data.frame} of failure times used in the fit.}
#' }
#'
#' @examples
#'
#' ## Single failure type examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' trt <- rbinom(n,1,0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
#'
#' ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n, 0, 1))
#'
#' # Fit 1 - fit mean_tmle object with GLMs for treatment, censoring, failure
#' fit1 <- mean_tmle(ftime = ftime, ftype = ftype,
#'                   trt = trt, adjustVars = adjustVars,
#'                   glm.trt = "W1 + W2",
#'                   glm.ftime = "trt + W1 + W2",
#'                   glm.ctime = "trt + W1 + W2")
#'
#' @export
#'

mean_tmle <- function(ftime,
                      ftype,
                      trt,
                      t0 = max(ftime[ftype > 0]),
                      adjustVars = NULL,
                      SL.ftime = NULL,
                      SL.ctime = NULL,
                      SL.trt = NULL,
                      glm.ftime = NULL,
                      glm.ctime = NULL,
                      glm.trt = "1",
                      glm.family = "binomial",
                      returnIC = TRUE,
                      returnModels = FALSE,
                      ftypeOfInterest = unique(ftype[ftype != 0]),
                      trtOfInterest = unique(trt),
                      bounds = NULL,
                      verbose = FALSE,
                      Gcomp = FALSE,
                      gtol = 1e-3,
                      ...) {

  # assemble data frame of necessary variables
  n <- length(ftime)
  id <- seq_len(n)
  dat <- data.frame(id = id, ftime = ftime, ftype = ftype, trt = trt)

  if (!is.null(adjustVars)) {
    dat <- cbind(dat, adjustVars)
  }

  # calculate number of failure types
  nJ <- length(ftypeOfInterest)
  allJ <- sort(unique(ftype[ftype != 0]))
  ofInterestJ <- sort(ftypeOfInterest)

  # calculate number of groups
  ntrt <- length(trtOfInterest)
  uniqtrt <- sort(trtOfInterest)

  # estimate trt probabilities
  trtOut <- estimateTreatment(
    dat = dat,
    ntrt = ntrt,
    uniqtrt = uniqtrt,
    adjustVars = adjustVars,
    SL.trt = SL.trt,
    glm.trt = glm.trt,
    returnModels = returnModels,
    gtol = gtol
  )
  dat <- trtOut$dat
  trtMod <- trtOut$trtMod

  # make long version of data sets needed for estimation of censoring
  dataList <- makeDataList(
    dat = dat, J = allJ, ntrt = ntrt, uniqtrt = uniqtrt,
    t0 = t0, bounds = bounds
  )

  # estimate censoring
  censOut <- estimateCensoring(
    dataList = dataList,
    ntrt = ntrt,
    uniqtrt = uniqtrt,
    t0 = t0,
    verbose = verbose,
    adjustVars = adjustVars,
    SL.ctime = SL.ctime,
    glm.ctime = glm.ctime,
    glm.family = glm.family,
    returnModels = returnModels,
    gtol = gtol
  )
  dataList <- censOut$dataList
  ctimeMod <- censOut$ctimeMod

  wideDataList <- makeWideDataList(
    dat = dat, dataList = dataList,
    adjustVars = adjustVars, t0 = t0,
    allJ = allJ, ntrt = ntrt, uniqtrt = uniqtrt
  )

  # estimate/fluctuate iterated means
  timeAndType <- expand.grid(rev(seq_len(t0)), ofInterestJ)

  # empty list for Qmod if returnModels
  ftimeMod <- vector(mode = "list", length = length(ofInterestJ))
  names(ftimeMod) <- paste0("J", ofInterestJ)
  for (j in seq_along(ofInterestJ)) {
    ftimeMod[[j]] <- vector(mode = "list", length = t0)
    names(ftimeMod[[j]]) <- paste0("t", seq_len(t0))
  }

  for (i in seq_len(nrow(timeAndType))) {
    estOut <- estimateIteratedMean(
      wideDataList = wideDataList,
      t = timeAndType[i, 1],
      whichJ = timeAndType[i, 2],
      ntrt = ntrt,
      uniqtrt = uniqtrt,
      allJ = allJ,
      t0 = t0,
      SL.ftime = SL.ftime,
      adjustVars = adjustVars,
      glm.ftime = glm.ftime,
      verbose = verbose,
      returnModels = returnModels,
      bounds = bounds
    )

    wideDataList <- estOut$wideDataList
    eval(parse(text = paste0(
      "ftimeMod$J", timeAndType[i, 2], "$t",
      timeAndType[i, 1], "<-estOut$ftimeMod"
    )))
    wideDataList <- fluctuateIteratedMean(
      wideDataList = wideDataList,
      t = timeAndType[i, 1],
      whichJ = timeAndType[i, 2],
      ntrt = ntrt, uniqtrt = uniqtrt,
      allJ = allJ, t0 = t0,
      SL.ftime = SL.ftime,
      glm.ftime = glm.ftime,
      returnModels = returnModels,
      bounds = bounds,
      Gcomp = Gcomp
    )
  }

  # get point estimate
  est <- rowNames <- NULL
  for (j in ofInterestJ) {
    for (z in seq_along(uniqtrt)) {
      thisEst <- eval(parse(text = paste(
        "mean(wideDataList[[", z + 1, "]]$Q",
        j, "star.1)", sep = ""
      )))
      est <- rbind(est, thisEst)
      rowNames <- c(rowNames, paste(c(uniqtrt[z], j), collapse = " "))
      eval(parse(text = paste(
        "wideDataList[[1]]$Q", j, "star.0.Z", uniqtrt[z],
        " <- rep(thisEst,n)", sep = ""
      )))
      eval(parse(text = paste(
        "wideDataList[[1]]$Q", j, "star.1.Z", uniqtrt[z],
        " <- wideDataList[[(z+1)]]$Q", j, "star.1",
        sep = ""
      )))
    }
  }
  row.names(est) <- rowNames

  # calculate influence function
  for (j in ofInterestJ) {
    for (z in seq_along(uniqtrt)) {
      for (t in rev(seq_len(t0))) {
        outcomeName <- ifelse(t == t0, paste("N", j, ".", t0, sep = ""),
          paste("Q", j, "star.", t + 1, sep = "")
        )
        eval(parse(text = paste(
          "wideDataList[[1]]$D.Z", uniqtrt[z], ".", j,
          "star.", t, " <- wideDataList[[1]]$H",
          uniqtrt[z], ".", t,
          "*(wideDataList[[1]][,outcomeName] - wideDataList[[1]]$Q",
          j, "star.", t, ")", sep = ""
        )))
      }
      eval(parse(text = paste(
        "wideDataList[[1]]$D.Z", uniqtrt[z], ".", j,
        "star.0 <- wideDataList[[1]]$Q", j, "star.1.Z",
        uniqtrt[z], " - wideDataList[[1]]$Q", j,
        "star.0.Z", uniqtrt[z], sep = ""
      )))
      ind <- eval(parse(text = paste(
        "grep('D.Z", uniqtrt[z], ".", j,
        "star', names(wideDataList[[1]]))",
        sep = ""
      )))
      eval(parse(text = paste(
        "wideDataList[[1]]$IC", j, "star.Z", uniqtrt[z],
        " <- rowSums(cbind(rep(0, nrow(wideDataList[[1]])),wideDataList[[1]][,ind]))",
        sep = ""
      )))
    }
  }

  # calculate standard error
  infCurves <- wideDataList[[1]][
    , grep("IC", names(wideDataList[[1]])),
    drop = FALSE
  ]
  meanIC <- apply(infCurves, MARGIN = 2, FUN = mean)
  var <- t(as.matrix(infCurves)) %*% as.matrix(infCurves) / (n ^ 2)
  row.names(var) <- colnames(var) <- rowNames

  out <- list(
    est = est, var = var, meanIC = meanIC, ic = infCurves,
    trtMod = trtMod, ftimeMod = ftimeMod, ctimeMod = ctimeMod,
    ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars
  )
  class(out) <- "survtmle"
  return(out)
}
