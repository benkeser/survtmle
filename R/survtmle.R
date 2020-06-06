#' Compute Targeted Minimum Loss-Based Estimators in Survival Analysis Settings
#'
#' @description This function estimates the marginal cumulative incidence for
#'  failures of specified types using targeted minimum loss-based estimation.
#'
#' @param ftime An integer-valued vector of failure times. Right-censored
#'  observations should have corresponding \code{ftype} set to 0.
#' @param ftype An integer-valued vector indicating the type of failure.
#'  Observations with \code{ftype=0} are treated as being right-censored. Each
#'  unique value besides zero is treated as a separate type of failure.
#' @param trt A numeric vector indicating observed treatment assignment. Each
#'  unique value will be treated as a different type of treatment. Currently,
#'  only two unique values are supported.
#' @param adjustVars A data.frame of adjustment variables that will be used in
#'  estimating the conditional treatment, censoring, and failure (hazard or
#'  conditional mean) probabilities.
#' @param t0 The time at which to return cumulative incidence estimates. By
#'  default this is set to \code{max(ftime[ftype > 0])}.
#' @param SL.ftime A character vector or list specification to be passed to the
#'  \code{SL.library} in the call to \code{\link[SuperLearner]{SuperLearner}}
#'  for the outcome regression (either cause-specific hazards or iterated
#'  mean). See the documentation of \code{\link[SuperLearner]{SuperLearner}}
#'  for more information on how to specify valid \code{SuperLearner} libraries.
#'  It is expected that the wrappers used in the library will play nicely with
#'  the input variables, which will be called \code{"trt"},
#'  \code{names(adjustVars)}, and \code{"t"} (if \code{method="hazard"}).
#' @param SL.ctime A character vector or list specification to be passed to the
#'  \code{SL.library} in the call to \code{\link[SuperLearner]{SuperLearner}}
#'  for the estimate of the conditional hazard for censoring. It is expected
#'  that the wrappers used in the library will play nicely with the input
#'  variables, which will be called \code{"trt"} and \code{names(adjustVars)}.
#' @param SL.trt A character vector or list specification to be passed to the
#'  \code{SL.library} in the call to \code{\link[SuperLearner]{SuperLearner}}
#'  for the estimate of the conditional probability of treatment. It is
#'  expected that the wrappers used in the library will play nicely with the
#'  input variables, which will be \code{names(adjustVars)}.
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
#'  interest. The default value computes estimates for all of the values in
#'  \code{unique(trt)}. Can alternatively be set to a vector of values found in
#'  \code{trt}.
#' @param cvControl A \code{list} providing control options to be fed directly
#'  into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
#'  the contents of \code{SuperLearner.CV.control} exactly. For details,
#'  consult the documentation of the \pkg{SuperLearner} package.
#' @param method A character specification of how the targeted minimum
#'  loss-based estimators should be computed, either \code{"mean"} or
#'  \code{"hazard"}. The \code{"mean"} specification uses a closed-form
#'  targeted minimum loss-based estimation based on the G-computation formula
#'  of Bang and Robins (2005). The \code{"hazard"} specification uses an
#'  iterative algorithm based on cause-specific hazard functions. The latter
#'  specification has no guarantee of convergence in finite samples. The
#'  convergence can be influenced by the stopping criteria specified in the
#'  \code{tol}. Future versions may implement a closed-form version of this
#'  hazard-based estimator.
#' @param bounds A \code{data.frame} of bounds on the conditional hazard
#'  function (if \code{method = "hazard"}) or on the iterated conditional means
#'  (if \code{method = "mean"}). The \code{data.frame} should have a column
#'  named \code{"t"} that includes values \code{seq_len(t0)}. The other columns
#'  should be names \code{paste0("l",j)} and \code{paste0("u",j)} for each
#'  unique failure type label j, denoting lower and upper bounds, respectively.
#'  See examples.
#' @param verbose A \code{logical} indicating whether the function should print
#'  messages to indicate progress. If \code{\link[SuperLearner]{SuperLearner}}
#'  is called internally, this option will be passed to it.
#' @param tol The stopping criteria when \code{method="hazard"}. The TMLE
#'  algorithm performs updates to the initial estimators until the empirical
#'  mean of the efficient influence function is smaller than \code{tol} or
#'  until \code{maxIter} iterations have been completed. The default
#'  (\code{1/length(ftime)}) is a sensible value. Larger values can be used in
#'  situations where convergence of the algorithm is an issue; however, this
#'  may result in large finite-sample bias.
#' @param Gcomp A \code{logical} indicating whether to compute the
#'  G-computation estimator (i.e., a substitution estimator with no targeting
#'  step). Theory does not support inference for the G-computation estimator if
#'  Super Learner is used to estimate failure and censoring distributions. The
#'  G-computation is only implemented for \code{method = "mean"}.
#' @param maxIter A maximum number of iterations for the algorithm when
#'  \code{method = "hazard"}. The algorithm will iterate until either the
#'  empirical mean of the efficient influence function is smaller than
#'  \code{tol} or until \code{maxIter} iterations have been completed.
#' @param gtol The truncation level of predicted censoring survival. Setting to
#'  larger values can help performance in data sets with practical positivity
#'  violations.
#' @param returnCall A \code{logical} specifying whether to return the function
#'  call via \code{match.call(expand.dots = TRUE)}. Set to \code{FALSE} to have
#'  the \code{call} slot return only \code{NA} when passing in pre-computed
#'  initial estimates to reduce memory inefficiency overhead. Defaults to
#'  \code{TRUE}.
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
#'   \item{ftimeMod}{If \code{returnModels=TRUE} the fit object(s) for the call
#'     to \code{\link[stats]{glm}} or \code{\link[SuperLearner]{SuperLearner}}
#'     for the outcome regression models. If \code{method="mean"} this will be
#'     a list of length \code{length(ftypeOfInterest)} each of length \code{t0}
#'     (one regression for each failure type and for each timepoint). If
#'     \code{method="hazard"} this will be a list of length
#'     \code{length(ftypeOfInterest)} with one fit corresponding to the hazard
#'     for each cause of failure. If \code{returnModels = FALSE}, this entry
#'     will be \code{NULL}.}
#'   \item{ctimeMod}{If \code{returnModels=TRUE} the fit object for the call to
#'     \code{\link[stats]{glm}} or \code{\link[SuperLearner]{SuperLearner}} for
#'     the pooled hazard regression model for the censoring distribution. If
#'     \code{returnModels=FALSE}, this entry will be \code{NULL}.}
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
#' @importFrom SuperLearner SuperLearner.CV.control
#'
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 200
#' trt <- rbinom(n, 1, 0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
#'
#' ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n, 0, 1))
#'
#' # Fit 1
#' # fit a survtmle object with glm estimators for treatment, censoring, and
#' # failure using the "mean" method
#' fit1 <- survtmle(
#'   ftime = ftime, ftype = ftype,
#'   trt = trt, adjustVars = adjustVars,
#'   glm.trt = "W1 + W2",
#'   glm.ftime = "trt + W1 + W2",
#'   glm.ctime = "trt + W1 + W2",
#'   method = "mean", t0 = 6
#' )
#' fit1
#'
#' # Fit 2
#' # fit an survtmle object with SuperLearner estimators for failure and
#' # censoring and empirical estimators for treatment using the "mean" method
#' fit2 <- survtmle(
#'   ftime = ftime, ftype = ftype,
#'   trt = trt, adjustVars = adjustVars,
#'   SL.ftime = c("SL.mean"),
#'   SL.ctime = c("SL.mean"),
#'   method = "mean", t0 = 6
#' )
#' fit2
#' @export
survtmle <- function(ftime, ftype, trt, adjustVars, t0 = max(ftime[ftype > 0]),
                     SL.ftime = NULL, SL.ctime = NULL, SL.trt = NULL,
                     glm.ftime = NULL, glm.ctime = NULL, glm.trt = NULL,
                     returnIC = TRUE, returnModels = TRUE,
                     ftypeOfInterest = unique(ftype[ftype != 0]),
                     trtOfInterest = unique(trt),
                     cvControl = list(
                       V = 10L, stratifyCV = FALSE,
                       shuffle = TRUE, validRows = NULL
                     ),
                     method = "hazard", bounds = NULL, verbose = FALSE,
                     tol = 1 / (sqrt(length(ftime))),
                     maxIter = 10, Gcomp = FALSE, gtol = 1e-3,
                     returnCall = TRUE) {
  # optionally catch function call
  if (returnCall) {
    call <- match.call(expand.dots = TRUE)
  } else {
    call <- NA
  }

  # check and clean inputs
  clean <- checkInputs(
    ftime = ftime, ftype = ftype, trt = trt,
    t0 = t0, adjustVars = adjustVars,
    SL.ftime = SL.ftime,
    SL.ctime = SL.ctime,
    SL.trt = SL.trt,
    glm.ftime = glm.ftime,
    glm.ctime = glm.ctime,
    glm.trt = glm.trt,
    returnIC = returnIC,
    returnModels = returnModels,
    ftypeOfInterest = ftypeOfInterest,
    trtOfInterest = trtOfInterest,
    bounds = bounds, verbose = verbose, tol = tol,
    Gcomp = Gcomp, method = method
  )

  # run cross-validation arguments through reformatting helper
  cvControl <- do.call(SuperLearner::SuperLearner.CV.control, cvControl)

  # hazard-based TMLE
  if (method == "hazard") {
    tmle.fit <- hazard_tmle(
      ftime = clean$ftime,
      ftype = clean$ftype,
      trt = clean$trt,
      t0 = t0,
      adjustVars = clean$adjustVars,
      SL.ftime = clean$SL.ftime,
      SL.ctime = clean$SL.ctime,
      SL.trt = clean$SL.trt,
      glm.ftime = clean$glm.ftime,
      glm.ctime = clean$glm.ctime,
      glm.trt = clean$glm.trt,
      returnIC = returnIC,
      returnModels = returnModels,
      ftypeOfInterest = ftypeOfInterest,
      trtOfInterest = trtOfInterest,
      cvControl = cvControl,
      bounds = bounds,
      verbose = verbose,
      tol = tol,
      maxIter = maxIter,
      gtol = gtol
    )
  } else if (method == "mean") {
    tmle.fit <- mean_tmle(
      ftime = clean$ftime,
      ftype = clean$ftype,
      trt = clean$trt,
      t0 = t0,
      adjustVars = clean$adjustVars,
      SL.ftime = clean$SL.ftime,
      SL.ctime = clean$SL.ctime,
      SL.trt = clean$SL.trt,
      glm.ftime = clean$glm.ftime,
      glm.ctime = clean$glm.ctime,
      glm.trt = clean$glm.trt,
      returnIC = returnIC,
      returnModels = returnModels,
      ftypeOfInterest = ftypeOfInterest,
      trtOfInterest = trtOfInterest,
      cvControl = cvControl,
      bounds = bounds,
      verbose = verbose,
      tol = tol,
      Gcomp = Gcomp,
      gtol = gtol
    )
  }

  out <- list(
    call = call, est = tmle.fit$est, var = tmle.fit$var,
    meanIC = tmle.fit$meanIC, ic = tmle.fit$ic,
    ftimeMod = tmle.fit$ftimeMod, ctimeMod = tmle.fit$ctimeMod,
    trtMod = tmle.fit$trtMod, t0 = t0, ftime = tmle.fit$ftime,
    ftype = tmle.fit$ftype, trt = tmle.fit$trt,
    adjustVars = tmle.fit$adjustVars
  )
  class(out) <- "survtmle"
  return(out)
}
