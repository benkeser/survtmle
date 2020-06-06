#' Estimate Treatment Mechanisms
#'
#' @description This function computes the conditional probability of having
#'  \code{trt} for each specified level either using \code{\link[stats]{glm}}
#'  or \code{\link[SuperLearner]{SuperLearner}}. Currently, only two unique
#'  values of treatment are acceptable. By default the function will compute
#'  estimates of the conditional probability of \code{trt == max(trt)} and
#'  compute the probability of \code{trt == min(trt)} as one minus this
#'  probability.
#'
#' @param dat An object of class \code{data.frame}. Must have named column
#'  \code{trt}.
#' @param adjustVars An object of class \code{data.frame} that will be used
#'  either as the \code{data} argument in a call to \code{\link[stats]{glm}} or
#'  as \code{X} in a call to \code{\link[SuperLearner]{SuperLearner}}.
#' @param glm.trt A character formula for the right-hand side of the
#'  \code{\link[stats]{formula}} in a call to \code{\link[stats]{glm}}. See the
#'  documentation of \code{\link{survtmle}} for information. Alternatively,
#'  this could be an object of class \code{\link[stats]{glm}} (as in calls to
#'  this function via \code{\link{timepoints}}), in which case predictions are
#'  obtained using this object with no new fitting.
#' @param SL.trt A specification of the \code{SL.library} option of a call to
#'  \code{\link[SuperLearner]{SuperLearner}}. See the documentation of
#'  \code{\link{survtmle}} for more information. Alternatively, this could be
#'  an object of class \code{SuperLearner} (as in calls to this function via
#'  \code{\link{timepoints}}), in which case predictions are obtained using
#'  this object with no new fitting.
#' @param cvControl A \code{list} providing control options to be fed directly
#'  into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
#'  the contents of \code{SuperLearner.CV.control} exactly. For details,
#'  consult the documentation of the \pkg{SuperLearner} package. This is passed
#'  in from \code{\link{mean_tmle}} or \code{\link{hazard_tmle}} via
#'  \code{\link{survtmle}}.
#' @param returnModels A \code{logical} indicating whether fitted model objects
#'  should be returned.
#' @param verbose A \code{logical} passed to the \code{verbose} option of the
#'  call to \code{\link[SuperLearner]{SuperLearner}}.
#' @param gtol The truncation level of predicted trt probabilities to handle
#'  positivity violations.
#' @param ... Other arguments. Not currently used
#'
#' @return dat The input \code{data.frame} object with two added columns
#'  corresponding with the conditional probability (given \code{adjustVars}) of
#'  \code{trt==max(trt)} and \code{trt==min(trt)}.
#' @return trtMod If \code{returnModels = TRUE}, the fitted \code{glm} or
#'  \code{SuperLearner} object. Otherwise, \code{NULL}
#'
#' @importFrom stats as.formula predict model.matrix optim glm
#' @importFrom SuperLearner SuperLearner SuperLearner.CV.control All SL.mean
#'  SL.glm SL.step
#'
#' @export
estimateTreatment <- function(dat,
                              adjustVars,
                              glm.trt = NULL,
                              SL.trt = NULL,
                              cvControl,
                              returnModels = FALSE,
                              verbose = FALSE,
                              gtol = 1e-3,
                              ...) {
  if (length(unique(dat$trt)) == 1) {
    eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
  } else {
    # binarize the outcome
    thisY <- as.numeric(dat$trt == max(dat$trt))

    # fit Super Learner
    if (!is.null(SL.trt)) {
      if (class(SL.trt) != "SuperLearner") {
        trtMod <- SuperLearner::SuperLearner(
          Y = thisY, X = adjustVars,
          newX = adjustVars,
          SL.library = SL.trt,
          id = dat$id, verbose = verbose,
          family = "binomial",
          cvControl = cvControl
        )
      } else {
        trtMod <- SL.trt
      }
      dat[[paste0("g_", max(dat$trt))]] <- trtMod$SL.predict
      dat[[paste0("g_", min(dat$trt))]] <- 1 - trtMod$SL.predict
    } else if (!is.null(glm.trt) & is.null(SL.trt)) {
      # set up model formula and data for the treatment regression
      trt_form <- paste("thisY", "~", glm.trt, sep = " ")
      trt_data_in <- as.data.frame(cbind(adjustVars, thisY))

      # fit GLM if Super Learner not requested
      if (!("glm" %in% class(glm.trt)) & !("speedglm" %in% class(glm.trt))) {
        # fit the treatment model
        trtMod <- fast_glm(
          reg_form = stats::as.formula(trt_form),
          data = trt_data_in,
          family = stats::binomial()
        )
      } else {
        trtMod <- glm.trt
      }
      suppressWarnings(
        pred <- predict(trtMod, newdata = trt_data_in, type = "response")
      )
      dat[[paste0("g_", max(dat$trt))]] <- pred
      dat[[paste0("g_", min(dat$trt))]] <- 1 - pred
    }
  }

  # truncate propensities
  eval(parse(text = paste0(
    "dat$g_", min(dat$trt), "[dat$g_", min(dat$trt),
    "< gtol]<- gtol"
  )))
  eval(parse(text = paste0(
    "dat$g_", max(dat$trt), "[dat$g_", max(dat$trt),
    "< gtol]<- gtol"
  )))
  out <- list()
  out$dat <- dat
  out$trtMod <- NULL
  if (returnModels) out$trtMod <- trtMod
  return(out)
}
