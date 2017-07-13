#' Estimate Treatment Mechanisms
#'
#' This function computes the conditional probability of having \code{trt} for
#' each specified level either using \code{glm} or \code{SuperLearner}.
#' Currently only two unique values of treatment are acceptable. By default the
#' function will compute estimates of the conditional probability of
#' \code{trt == max(trt)} and compute the probability of \code{trt == min(trt)}
#' as one minus this probability.
#'
#' @param dat An object of class \code{data.frame}. Must have named column
#'        \code{trt}.
#' @param adjustVars An object of class \code{data.frame} that will be used
#'        either as the \code{data} argument in a call to \code{glm} or as the
#'        \code{X} object in a call to \code{SuperLearner}.
#' @param glm.trt A character formula for the right-hand side of \code{formula}
#'        in a call to \code{glm}. See \code{?survtmle} for more documentation.
#'        Alternatively, this could be an object of class \code{glm} (as in
#'        calls to this function via \code{timepoints}), in which case
#'        predictions are obtained using this object with no new fitting.
#' @param SL.trt A valid specification of the \code{SL.library} option of a call
#'        to \code{SuperLearner}. See \code{?survtmle} for more documentation.
#'        Alternatively, this could be an object of class \code{SuperLearner}
#'        (as in calls to this function via \code{timepoints}), in which case
#'        predictions are obtained using this object with no new fitting.
#' @param returnModels A boolean indicating whether fitted model objects should
#'        be returned.
#' @param verbose A boolean passed to the \code{verbose} option of the call to
#'        \code{SuperLearner}.
#' @param gtol The truncation level of predicted trt probabilities to handle
#'        positivity violations.
#' @param ... Other arguments. Not currently used
#'
#' @return dat The input \code{data.frame} object with two added columns
#'         corresponding with the conditional probability (given
#'         \code{adjustVars}) of \code{trt==max(trt)} and \code{trt==min(trt)}.
#' @return trtMod If \code{returnModels = TRUE}, the fitted \code{glm} or
#'         \code{SuperLearner} object. Otherwise, \code{NULL}
#'
#' @importFrom stats as.formula predict model.matrix optim glm
#' @importFrom SuperLearner SuperLearner SuperLearner.CV.control All SL.mean SL.glm SL.step
#'

estimateTreatment <- function(dat, adjustVars, glm.trt = NULL, SL.trt = NULL,
                              returnModels = FALSE, verbose = FALSE,
                              gtol = 1e-3, ...) {
  if(length(unique(dat$trt)) == 1) {
    eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
  } else {
    if(!is.null(SL.trt)) {
      if(class(SL.trt) != "SuperLearner") {
        thisY <- as.numeric(dat$trt == max(dat$trt))
        trtMod <- SuperLearner::SuperLearner(Y = thisY, X = adjustVars,
                                             newX = adjustVars,
                                             SL.library = SL.trt,
                                             id = dat$id, verbose = verbose,
                                             family = "binomial")
      } else {
        trtMod <- SL.trt
      }
      dat[[paste0("g_",max(dat$trt))]] <- trtMod$SL.predict
      dat[[paste0("g_",min(dat$trt))]] <- 1-trtMod$SL.predict

    }else if(!is.null(glm.trt) & is.null(SL.trt)) {
      if(!("glm" %in% class(glm.trt))) {
        thisY <- as.numeric(dat$trt == max(dat$trt))
        trtMod <- stats::glm(stats::as.formula(paste0("thisY ~ ", glm.trt)),
                             data = adjustVars, family = "binomial")
      } else {
        trtMod <- glm.trt
      }
      suppressWarnings(
        pred <- predict(trtMod, type = "response")
      )
      dat[[paste0("g_",max(dat$trt))]] <- pred
      dat[[paste0("g_",min(dat$trt))]] <- 1-pred
    }
  }

  # truncate propensities
  eval(parse(text = paste0("dat$g_", min(dat$trt), "[dat$g_", min(dat$trt),
                           "< gtol]<- gtol")))
  eval(parse(text=paste0("dat$g_", max(dat$trt), "[dat$g_", max(dat$trt),
                         "< gtol]<- gtol")))

  out <- list(dat = dat,
              trtMod = if(returnModels & length(unique(dat$trt)) > 1)
                trtMod
              else
                NULL)
  out
}
