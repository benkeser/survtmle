#' TMLE for Cause-Specific Hazard Functions
#'
#' @description This function estimates the marginal cumulative incidence for
#'  failures of specified types using targeted minimum loss-based estimation
#'  based on the initial estimates of the cause-specific hazard functions for
#'  failures of each type. The function is called by \code{\link{survtmle}}
#'  whenever \code{method = "hazard"} is specified. However, power users could,
#'  in theory, make calls directly to this function.
#'
#' @param ftime A numeric vector of failure times. Right-censored observations
#'  should have corresponding \code{ftype} set to 0.
#' @param ftype A numeric vector indicating the type of failure. Observations
#'  with \code{ftype=0} are treated as having been right-censored. Each unique
#'  value besides zero is treated as a separate type of failure.
#' @param trt A numeric vector indicating observed treatment assignment. Each
#'  unique value will be treated as a different type of treatment. Currently,
#'  only two unique values are supported.
#' @param adjustVars A \code{data.frame} of adjustment variables that will be used in
#'  estimating the conditional treatment, censoring, and failure (hazard or
#'  conditional mean) probabilities.
#' @param mediator A \code{data.frame} of mediating variables. If non-\code{NULL}
#'  the code will compute estimates of the counterfactual cumulative incidence
#'  under a set of joint interventions that each in turn set \code{trt} to a 
#'  particular level of trtOfInterest and then set the mediator equal to the value
#'  it would have naturally assumed under \code{trt = mediatorTrtVal}. Alternatively,
#'  the estimates could be interpreted as estimates of "interventional mediation effects"
#'  wherein the latter intervention draws \code{mediator} from the observed data 
#'  distribution of of \code{mediator} given \code{trt = mediatorTrtVal} and 
#'  \code{adjustVars}.
#' @param mediatorTrtVal A \code{numeric} of length 1 indicating the value of
#'  the treatment under which the mediator is set.  
#' @param mediatorSampProb A \code{vector} of non-negative sampling weights for the
#'  mediator/s.
#' @param mediatorInCensMod A \code{boolean} indicating whether to include mediators
#'  in the censoring model.
#' @param t0 The time at which to return cumulative incidence estimates. By
#'  default this is set to \code{max(ftime[ftype > 0])}.
#' @param SL.ftime A character vector or list specification to be passed to the
#'  \code{SL.library} option of \code{\link[SuperLearner]{SuperLearner}} for
#'  the cause-specific hazards. See the documentation of
#'  \code{\link[SuperLearner]{SuperLearner}} for more information on how to
#'  specify valid \code{SuperLearner} libraries. It is expected that the
#'  wrappers used in the library will play nicely with the input variables,
#'  which will be called \code{"trt"}, \code{names(adjustVars)}, and \code{"t"}
#'  if \code{method = "hazard"}.
#' @param SL.ctime A character vector or list specification to be passed to the
#'  \code{SL.library} option of \code{\link[SuperLearner]{SuperLearner}} for
#'  the estimate of the conditional hazard for censoring. It is expected that
#'  the wrappers used in the library will play nicely with the input variables,
#'  which will be called \code{"trt"} and \code{names(adjustVars)}.
#' @param SL.trt A character vector or list specification to be passed to the
#'  \code{SL.library} option of \code{\link[SuperLearner]{SuperLearner}} for
#'  the estimate of the conditional probability of treatment. It is expected
#'  that the wrappers used in the library will play nicely with the input
#'  variables, which will be \code{names(adjustVars)}.
#' @param glm.ftime A character specification of the right-hand side of the
#'  equation passed to the \code{\link[stats]{formula}} option of a call to
#'  \code{\link[stats]{glm}} for the outcome regression. Ignored if
#'  \code{SL.ftime} is not equal to \code{NULL}. Use \code{"trt"} to specify
#'  the treatment in this formula (see examples). The formula can additionally
#'  include any variables found in \code{names(adjustVars)}.
#' @param glm.ctime A character specification of the right-hand side of the
#'  equation passed to the \code{\link[stats]{formula}} option of a call to
#'  \code{\link[stats]{glm}} for the estimate of the conditional hazard for
#'  censoring. Ignored if \code{SL.ctime} is not equal to \code{NULL}. Use
#'  \code{"trt"} to specify the treatment in this formula (see examples). The
#'  formula can additionally include any variables found in
#'  \code{names(adjustVars)}.
#' @param glm.trt A character specification of the right-hand side of the
#'  equation passed to the \code{\link[stats]{formula}} option of a call to
#'  \code{\link[stats]{glm}} for the estimate of the conditional probability of
#'  treatment. Ignored if \code{SL.trt} is not equal to \code{NULL}. The
#'  formula can include any variables found in \code{names(adjustVars)}.
#' @param glm.family The type of regression to be performed if fitting GLMs in
#'  the estimation and fluctuation procedures. The default is "binomial" for
#'  logistic regression. Only change this from the default if there are
#'  justifications that are well understood. This is passed directly to
#'  \code{\link{estimateCensoring}} and \code{\link{estimateHazards}}.
#' @param att A \code{boolean} indicating whether to compute the ATT estimate,
#'  instead of treatment specific survival curves. This option only works with 
#'  two levels of \code{trt} that are labeled with 0 and 1. 
#' @param returnIC A \code{logical} indicating whether to return vectors of
#'  influence curve estimates. These are needed for some post-hoc comparisons,
#'  so it is recommended to leave as \code{TRUE} (the default) unless the user
#'  is sure these estimates will not be needed later.
#' @param returnModels A \code{logical} indicating whether to return the
#'  \code{glm} or \code{SuperLearner} objects used to estimate the nuisance
#'  parameters. Must be set to \code{TRUE} if the user plans to use
#'  \code{\link{timepoints}} to obtain estimates of incidence at times other
#'  than \code{t0}. See the documentation of \code{\link{timepoints}} for more
#'  information.
#' @param ftypeOfInterest An input specifying what failure types to compute
#'  estimates of incidence for. The default value computes estimates for values
#'  \code{unique(ftype)}. Can alternatively be set to a vector of values found
#'  in \code{ftype}.
#' @param trtOfInterest An input specifying which levels of \code{trt} are of
#'  interest. The default value computes estimates for all values in
#'  \code{unique(trt)}. Can alternatively be set to a vector of values found in
#'  \code{trt}.
#' @param cvControl A \code{list} providing control options to be fed directly
#'  into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
#'  the contents of \code{SuperLearner.CV.control} exactly. For details,
#'  consult the documentation of the \pkg{SuperLearner} package. This is
#'  usually passed in through the \code{\link{survtmle}} wrapper function.
#' @param bounds A \code{data.frame} of bounds on the conditional hazard
#'  function. The \code{data.frame} should have a column named \code{"t"} that
#'  includes values \code{seq_len(t0)}. The other columns should be names
#'  \code{paste0("l",j)} and \code{paste0("u",j)} for each unique failure type
#'  label j, denoting lower and upper bounds, respectively. See examples.
#' @param verbose A \code{logical} indicating whether the function should print
#'  messages to indicate progress. If \code{\link[SuperLearner]{SuperLearner}}
#'  is called internally, this option will additionally be passed to
#'  \code{\link[SuperLearner]{SuperLearner}}.
#' @param tol The stopping criteria. The TMLE algorithm performs updates to the
#'  initial estimators until the empirical mean of the efficient influence
#'  function is smaller than \code{tol} or until \code{maxIter} iterations have
#'  been completed. The default (\code{1 / length(ftime)}) is a sensible value.
#'  Larger values can be used in situations where convergence of the algorithm
#'  is an issue; however, this may result in large finite-sample bias.
#' @param maxIter The maximum number of iterations for the algorithm. The
#'  algorithm will iterate until either the empirical mean of the efficient
#'  influence function is smaller than \code{tol} or until \code{maxIter}
#'  iterations have been completed.
#' @param gtol The truncation level of predicted censoring survival. Setting to
#'  larger values can help performance in data sets with practical positivity
#'  violations.
#' @param ... Other options. Not currently used.
#'
#' @return An object of class \code{survtmle}.
#' \describe{
#'   \item{call}{The call to \code{survtmle}.}
#'   \item{est}{A numeric vector of point estimates -- one for each combination
#'     of \code{ftypeOfInterest} and \code{trtOfInterest}.}
#'   \item{var}{A covariance matrix for the point estimates.}
#'   \item{meanIC}{The empirical mean of the efficient influence function at
#'     the estimated, targeted nuisance parameters. Each value should be small
#'     or the user will be warned that excessive finite-sample bias may exist
#'     in the point estimates.}
#'   \item{ic}{The efficient influence function at the estimated, fluctuated
#'     nuisance parameters, evaluated on each of the observations. These are
#'     used to construct confidence intervals for post-hoc comparisons.}
#'   \item{ftimeMod}{If \code{returnModels = TRUE} the fit object(s) for calls
#'     to \code{\link[SuperLearner]{SuperLearner}} or \code{\link[stats]{glm}}
#'     for the outcome regression models. If \code{method="mean"} this will be
#'     a list of length \code{length(ftypeOfInterest)} each of length \code{t0}
#'     (one regression for each failure type and for each timepoint). If
#'     \code{method = "hazard"} this will be a list of length
#'     \code{length(ftypeOfInterest)} with one fit corresponding to the hazard
#'     for each cause of failure. If \code{returnModels = FALSE}, this entry
#'     will be \code{NULL}.}
#'   \item{ctimeMod}{If \code{returnModels = TRUE} the fit object for the call
#'     to \code{\link[stats]{glm}} or \code{\link[SuperLearner]{SuperLearner}}
#'     for the pooled hazard regression model for the censoring distribution.
#'     If \code{returnModels = FALSE}, this entry will be \code{NULL}.}
#'   \item{trtMod}{If \code{returnModels = TRUE} the fit object for the call to
#'     \code{\link[stats]{glm}} or \code{\link[SuperLearner]{SuperLearner}} for
#'     the conditional probability of \code{trt} regression model. If
#'     \code{returnModels = FALSE}, this entry will be \code{NULL}.}
#'   \item{t0}{The timepoint at which the function was evaluated.}
#'   \item{ftime}{The \code{numeric} vector of failure times used in the fit.}
#'   \item{ftype}{The \code{numeric} vector of failure types used in the fit.}
#'   \item{trt}{The \code{numeric} vector of treatment assignments used in the
#'     fit.}
#'   \item{adjustVars}{The \code{data.frame} of failure times used in the fit.}
#' }
#'
#' @examples
#' ## Single failure type examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' trt <- rbinom(n, 1, 0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
#'
#' ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n, 0, 1))
#'
#' # Fit 1 - fit hazard_tmle object with GLMs for treatment, censoring, failure
#' fit1 <- hazard_tmle(
#'   ftime = ftime, ftype = ftype,
#'   trt = trt, adjustVars = adjustVars,
#'   glm.trt = "W1 + W2",
#'   glm.ftime = "trt + W1 + W2",
#'   glm.ctime = "trt + W1 + W2",
#'   returnModels = TRUE
#' )
#' @export

# TO DO: Add input for weight and ph2 indicator
# TO DO: Add option for referent trt arm for mediation? e.g., a_1, a_2
hazard_tmle <- function(ftime,
                        ftype,
                        trt,
                        t0 = max(ftime[ftype > 0]),
                        adjustVars = NULL,
                        mediator = NULL,
                        mediatorTrtVal = NULL,
                        mediatorSampProb = NULL,
                        mediatorInCensMod = FALSE, 
                        mediatorStratify.ftime = TRUE,                       
                        SL.ftime = NULL,
                        SL.ctime = NULL,
                        SL.trt = NULL,
                        SL.trtMediator = NULL,
                        SL.ftimeMissing = NULL,
                        SL.mediator = NULL,
                        SL.eif = NULL,
                        glm.eif = NULL,
                        glm.mediator = NULL,
                        glm.trtMediator = NULL,
                        glm.ftime = NULL,
                        glm.ctime = NULL,
                        glm.trt = "1",
                        glm.ftimeMissing = NULL, 
                        glm.family = "binomial",
                        att = FALSE,
                        returnIC = TRUE,
                        returnModels = FALSE,
                        ftypeOfInterest = unique(ftype[ftype != 0]),
                        trtOfInterest = unique(trt),
                        cvControl = list(),
                        bounds = NULL,
                        verbose = FALSE,
                        tol = 1 / (length(ftime)),
                        maxIter = 5,
                        gtol = 1e-3,
                        gtolCens = 0.05,
                        truncateH = 0.975,
                        ...) {

  # assemble data frame of necessary variables
  n <- length(ftime)
  id <- seq_len(n)
  dat <- data.frame(id = id, ftime = ftime, 
                    ftype = ftype, trt = trt)
  if (!is.null(mediatorSampProb)){
    R <- as.numeric(complete.cases(mediator))
    dat$sampWt <- R / mediatorSampProb
    mediatorSampWt <- R / mediatorSampProb
  }else{
    dat$sampWt <- 1
  }
  if (!is.null(adjustVars)) dat <- cbind(dat, adjustVars)
  if (!is.null(mediator)) dat <- cbind(dat, mediator)
  
  nJ <- length(ftypeOfInterest)
  allJ <- sort(unique(ftype[ftype != 0]))
  ofInterestJ <- sort(ftypeOfInterest)

  # calculate number of groups
  ntrt <- length(trtOfInterest)
  uniqtrt <- sort(trtOfInterest)

  # TO DO: Move to check inputs
  if(!is.null(mediator)){
    stopifnot(length(uniqtrt) == 1)
  }

  # estimate trt probabilities
  trtOut <- estimateTreatment(
    dat = dat,
    ntrt = ntrt,
    uniqtrt = uniqtrt,
    adjustVars = adjustVars,
    SL.trt = SL.trt,
    glm.trt = glm.trt,
    cvControl = cvControl,
    returnModels = returnModels,
    gtol = gtol,
    trtOfInterest = trtOfInterest
  )
  dat <- trtOut$dat
  trtMod <- trtOut$trtMod

  if(!is.null(mediator)){
    trtMedOut <- estimateTreatment(
      dat = dat,
      ntrt = ntrt,
      uniqtrt = uniqtrt,
      adjustVars = adjustVars,
      mediator = mediator,
      mediatorSampWt = mediatorSampWt,
      SL.trt = SL.trtMediator,
      glm.trt = glm.trtMediator,
      cvControl = cvControl,
      returnModels = returnModels,
      gtol = gtol,
      trtOfInterest = trtOfInterest
    )
    dat <- trtMedOut$dat
    trtMedMod <- trtMedOut$trtMod
  }
  ftimeMissingOut <- estimateFtimeMissing(
    dat = dat,
    adjustVars = adjustVars,
    glm.ftimeMissing = glm.ftimeMissing,
    SL.ftimeMissing = SL.ftimeMissing,
    cvControl = cvControl,
    returnModels = returnModels,
    verbose = verbose,
    gtol = gtol,
    trtOfInterest = trtOfInterest,
    mediator = mediator,
    mediatorSampWt = mediatorSampWt,
    mediatorInCensMod = mediatorInCensMod
  )
  dat <- ftimeMissingOut$dat
  ftimeMissingMod <- ftimeMissingOut$ftimeMissingMod

  # TO DO: Move to internal estimateCensoring
  # hacky way to make long format data
  if(any(is.na(dat$ftime))){
    na_idx <- dat$id[is.na(dat$ftime)]
    dat$ftime[is.na(dat$ftime)] <- min(dat$ftime, na.rm = TRUE)
  }else{
    na_idx <- NULL
  }

  # make long version of data sets needed for estimation and prediction
  dataList <- makeDataList(
    dat = dat, J = allJ, ntrt = ntrt, uniqtrt = uniqtrt,
    t0 = t0, bounds = bounds
  )
  
  # TO DO: Move to internal estimateCensoring
  dataList[[1]] <- dataList[[1]][!(dataList[[1]]$id %in% na_idx), ]
  dat$ftime[dat$id %in% na_idx] <- NA

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
    cvControl = cvControl,
    returnModels = returnModels,
    gtol = gtolCens,
    mediator = mediator,
    mediatorSampWt = mediatorSampWt,
    mediatorInCensMod = mediatorInCensMod
  )
  dataList <- censOut$dataList
  ctimeMod <- censOut$ctimeMod

  # estimate cause specific hazards
  estOut <- estimateHazards(
    dataList = dataList,
    J = allJ,
    verbose = verbose,
    bounds = bounds,
    adjustVars = adjustVars,
    SL.ftime = SL.ftime,
    glm.ftime = glm.ftime,
    glm.family = glm.family,
    cvControl = cvControl,
    returnModels = returnModels,
    mediator = mediator,
    mediatorSampWt = mediatorSampWt,
    mediatorStratify.ftime = mediatorStratify.ftime,
    trtOfInterest = trtOfInterest
  )
  dataList <- estOut$dataList
  ftimeMod <- estOut$ftimeMod

  # calculate cum inc and clever covariates needed for fluctuations
  dataList <- updateVariables(
    dataList = dataList, allJ = allJ,
    ofInterestJ = ofInterestJ,
    nJ = nJ, uniqtrt = uniqtrt, ntrt = ntrt,
    t0 = t0, verbose = verbose, att = att,
    mediator = mediator, mediatorTrtVal = mediatorTrtVal,
    truncateH = truncateH
  )

  # TO DO: Note this is the full data influence function
  # calculate influence function
  dat <- getHazardInfluenceCurve(
    dataList = dataList, dat = dat,
    ofInterestJ = ofInterestJ, allJ = allJ,
    nJ = nJ, uniqtrt = uniqtrt, ntrt = ntrt,
    verbose = verbose, t0 = t0, att = att,
    mediator = mediator
  )
  infCurves <- dat[, grep("D.j", names(dat)), drop = FALSE]
  meanIC <- colMeans(infCurves)
  
  if(!is.null(mediator)){
    medOut <- estimateMedReg(
      dat = dat, 
      SL.mediator = SL.mediator,
      glm.mediator = glm.mediator,
      adjustVars = adjustVars,
      cvControl = cvControl,
      returnModels = returnModels,
      bounds = bounds,
      ofInterestJ = ofInterestJ,
      mediatorTrtVal = mediatorTrtVal,
      uniqtrt = uniqtrt,
      verbose = verbose
    )
    medMod <- medOut$mediatorMod
    dat <- medOut$dat

    dat <- getMediatorInfluenceCurve(
      dat = dat, uniqtrt = uniqtrt,
      mediatorTrtVal = mediatorTrtVal,
      ofInterestJ = ofInterestJ
    )

    eifOut <- estimateFullDataEIF(
      dat = dat, mediator = mediator,
      ofInterestJ = ofInterestJ, 
      uniqtrt = uniqtrt,
      adjustVars = adjustVars, 
      SL.eif = SL.eif, glm.eif = glm.eif, 
      cvControl = cvControl,
      verbose = verbose
    )
    eifMod <- eifOut$eifMod
    dat <- eifOut$dat

    samp_eif <- getSamplingInfluenceCurve(
      dat = dat, ofInterestJ = ofInterestJ, 
      uniqtrt = uniqtrt, R = R,
      mediatorSampProb = mediatorSampProb
    )
    targeted_mediatorSampProb <- mediatorSampProb

    ct <- 0
    while(max(abs(colMeans(samp_eif))) > tol & ct < maxIter){      
      ct <- ct + 1
      targeted_mediatorSampProb <- fluctuateSampWt(
        dat = dat, ofInterestJ = ofInterestJ, 
        uniqtrt = uniqtrt, R = R,
        mediatorSampProb = targeted_mediatorSampProb
      )
      samp_eif <- getSamplingInfluenceCurve(
        dat = dat, ofInterestJ = ofInterestJ, 
        uniqtrt = uniqtrt, R = R,
        mediatorSampProb = targeted_mediatorSampProb
      )
    }
    dat$sampWt <- R / targeted_mediatorSampProb
    # add sampling eif to dat
    dat <- data.frame(dat, samp_eif)
    # TO DO: more elegant way to do this
    dataList <- lapply(dataList, function(x, R, targeted_mediatorSampProb){
      for(i in seq_len(dim(dat)[1])){
        x[x$id == i, "sampWt"] <- R[i] / targeted_mediatorSampProb[i]
        return(x)
      }
    }, R = R, targeted_mediatorSampProb = targeted_mediatorSampProb)
  }else{
    medMod <- eifMod <- NULL
  }

  attr(dataList, "fluc") <- rep(Inf, ntrt * nJ^2)
  ct <- 0
  while (any(abs(meanIC) > tol) & ct < maxIter) {
    ct <- ct + 1
    dataList <- fluctuateHazards(
      dataList = dataList, ofInterestJ = ofInterestJ,
      allJ = allJ, nJ = nJ,
      uniqtrt = uniqtrt, ntrt = ntrt,
      verbose = verbose, t0 = t0, att = att,
      mediator = mediator, mediatorTrtVal = mediatorTrtVal,
      mediatorSampWt = mediatorSampWt, truncateH = truncateH
    )

    # calculate influence function
    dat <- getHazardInfluenceCurve(
      dataList = dataList, dat = dat,
      ofInterestJ = ofInterestJ, allJ = allJ,
      nJ = nJ, uniqtrt = uniqtrt, ntrt = ntrt,
      verbose = verbose, t0 = t0, att = att,
      mediator = mediator
    )
    infCurves <- dat[, grep("D.j", names(dat)), drop = FALSE]
    meanIC <- colMeans(infCurves)

    if (verbose) {
      # print(attr(dataList,"fluc"))
      cat("TMLE Iteration ", ct, "  : ", round(meanIC, 4), "\n")
    }
  }
  if (ct == maxIter + 1) {
    warning("maxIter reached. You may wish to check meanIC to ensure
             it is adequately small.")
  }

  if(!is.null(mediator)){
    medOut <- estimateMedReg(
      dat = dat, 
      SL.mediator = SL.mediator,
      glm.mediator = glm.mediator,
      adjustVars = adjustVars,
      cvControl = cvControl,
      returnModels = returnModels,
      bounds = bounds,
      ofInterestJ = ofInterestJ,
      mediatorTrtVal = mediatorTrtVal,
      uniqtrt = uniqtrt,
      verbose = verbose
    )
    medMod <- medOut$mediatorMod
    dat <- medOut$dat

    dat <- getMediatorInfluenceCurve(
      dat = dat, uniqtrt = uniqtrt,
      mediatorTrtVal = mediatorTrtVal,
      ofInterestJ = ofInterestJ
    )

    if(abs(mean(dat$Dmed.j1.z1)) > tol){
      dat <- fluctuateCumInc(
        dat = dat, ofInterestJ = ofInterestJ, 
        uniqtrt = uniqtrt, mediatorTrtVal = mediatorTrtVal
      )
    }

    dat <- getMediatorInfluenceCurve(
      dat = dat, uniqtrt = uniqtrt,
      mediatorTrtVal = mediatorTrtVal,
      ofInterestJ = ofInterestJ
    )
      # calculate point estimate
    est <- rowNames <- ic <- NULL 
    for (j in ofInterestJ) {
      for (z in uniqtrt) {
        est <- rbind(est, mean(dat[[paste0("margMedF", j, ".z", z, ".t0")]]))
        rowNames <- c(rowNames, paste(z, j))
        # browser()
        tmpD <- 
          # hazard component          
          dat[[paste0("D.j", j, ".z", z)]] + 
            # mediator component
            dat[[paste0("Dmed.j", j, ".z", z)]] -
              # sampling component
              dat[[paste0("Dsamp.j", j, ".z", z)]]
        if(att){
          # scale by P(trt = 1)
          tmpD <- tmpD / mean(dat$trt == 1) + 
            # adjust vars component
            as.numeric(dat$trt == 1) / mean(dat$trt == 1) * (
              dat[[paste0("margMedF", j, ".z", z, ".t0")]] - 
                mean(dat[[paste0("margMedF", j, ".z", z, ".t0")]])
            )
        }else{
          tmpD <- tmpD + 
            #adjustVars component
            dat[[paste0("margMedF", j, ".z", z, ".t0")]] - 
                mean(dat[[paste0("margMedF", j, ".z", z, ".t0")]])
        }
        ic <- cbind(ic, tmpD)
      }
    }
    # browser()
    row.names(est) <- rowNames
    var <- cov(ic) / n
    row.names(var) <- colnames(var) <- rowNames
  }else{
    # calculate point estimate
    est <- rowNames <- NULL
    for (j in ofInterestJ) {
      for (z in uniqtrt) {
        eval(parse(text = paste(
          "est <- rbind(est, dat$margF", j, ".z", z,
          ".t0[1])",
          sep = ""
        )))
        rowNames <- c(rowNames, paste(c(z, j), collapse = " "))
      }
    }
    row.names(est) <- rowNames

    # calculate standard error
    var <- t(as.matrix(infCurves)) %*% as.matrix(infCurves) / n^2
    row.names(var) <- colnames(var) <- rowNames
  }

  out <- list(
    est = est, var = var, meanIC = meanIC, ic = infCurves,
    trtMod = trtMod, ftimeMod = ftimeMod, ctimeMod = ctimeMod,
    # TO DO: add other Mods to output
    ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars
  )
  class(out) <- "survtmle"
  return(out)
}
