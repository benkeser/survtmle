#' Estimate Censoring Mechanisms
#'
#' @description Computes an estimate of the hazard for censoring using either
#'  \code{\link[stats]{glm}} or \code{\link[SuperLearner]{SuperLearner}} based
#'  on log-likelihood loss. The function then computes the censoring survival
#'  distribution based on these estimates. The structure of the function is
#'  specific to how it is called within \code{\link{survtmle}}. In particular,
#'  \code{dataList} must have a very specific structure for this function to
#'  run properly. The list should consist of \code{data.frame} objects. The
#'  first will have the number of rows for each observation equal to the
#'  \code{ftime} corresponding to that observation. Subsequent entries will
#'  have \code{t0} rows for each observation and will set \code{trt} column
#'  equal to each value of \code{trtOfInterest} in turn. One of these columns
#'  must be named \code{C} that is a counting process for the right-censoring
#'  variable. The function will fit a regression with \code{C} as the outcome
#'  and functions of \code{trt} and \code{names(adjustVars)} as specified by
#'  \code{glm.ctime} or \code{SL.ctime} as predictors.
#'
#' @param dataList A list of \code{data.frame} objects as described in the
#'  documentation of \code{\link{makeDataList}}.
#' @param adjustVars Object of class \code{data.frame} that contains the
#'  variables to adjust for in the regression.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#'  Needed only because the naming convention for the regression if
#'  \code{t == t0} is different than if \code{t != t0}.
#' @param SL.ctime A character vector or list specification to be passed to the
#'  \code{SL.library} argument of \code{\link[SuperLearner]{SuperLearner}} for
#'  the outcome regression (either cause-specific hazards or conditional mean).
#'  See the documentation of \code{\link[SuperLearner]{SuperLearner}} for more
#'  information on how to specify valid \code{SuperLearner} libraries. It is
#'  expected that the wrappers used in the library will play nicely with the
#'  input variables, which will be called \code{"trt"} and
#'  \code{names(adjustVars)}.
#' @param glm.ctime A character specification of the right-hand side of the
#'  equation passed to the \code{\link[stats]{formula}} option of a call to
#'  \code{\link[stats]{glm}} for the outcome regression (either cause-specific
#'  hazards or conditional mean). Ignored if \code{SL.ctime != NULL}. Use
#'  \code{"trt"} to specify the treatment in this formula (see examples). The
#'  The formula can additionally include any variables found in
#'  \code{names(adjustVars)}.
#' @param glm.family The type of regression to be performed if fitting GLMs in
#'  the estimation and fluctuation procedures. The default is "binomial" for
#'  logistic regression. Only change this from the default if there are
#'  justifications that are well understood. This is inherited from the calling
#'  function (either \code{\link{mean_tmle}} or \code{\link{hazard_tmle}}).
#' @param cvControl A \code{list} providing control options to be fed directly
#'  into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
#'  the contents of \code{SuperLearner.CV.control} exactly. For details,
#'  consult the documentation of the \pkg{SuperLearner} package. This is passed
#'  in from \code{\link{mean_tmle}} or \code{\link{hazard_tmle}} via
#'  \code{\link{survtmle}}.
#' @param returnModels A \code{logical} indicating whether to return the
#'  \code{glm} or \code{SuperLearner} objects used to estimate the nuisance
#'  parameters. Must be set to \code{TRUE} to make downstream calls to
#'  \code{\link{timepoints}} for obtaining estimates at times other than
#'  \code{t0}. See documentation of \code{\link{timepoints}} for more
#'  information.
#' @param gtol The truncation level of predicted censoring survival to handle
#'  positivity violations.
#' @param verbose A \code{logical} indicating whether the function should print
#'  messages to indicate progress.
#' @param stratify If \code{TRUE}, then the censoring hazard model is estimated using only
#'  the observations with \code{trt == trtOfInterest}. Only works if 
#'  \code{length(trtOfInterest) == 1}. If \code{stratify = TRUE} then \code{glm.ftime}
#'  cannot include \code{trt} in the model formula and any learners in \code{SL.ftime}
#'  should not assume a variable named \code{trt} will be included in the candidate 
#'  super learner estimators.
#' @param trtOfInterest An input specifying which levels of \code{trt} are of
#'  interest. The default value computes estimates for all values in
#'  \code{unique(trt)}. Can alternatively be set to a vector of values found in
#'  \code{trt}. Ignored unless \code{stratify == TRUE}.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats as.formula predict glm
#' @importFrom SuperLearner SuperLearner
#'
#' @return The function returns a list that is exactly the same as the input
#'  \code{dataList}, but with a column named \code{G_dC} added to it, which is
#'  the estimated conditional survival distribution for the censoring variable
#'  evaluated at the each of the rows of each \code{data.frame} in
#'  \code{dataList}.
#'
#' @export
estimateCensoring <- function(dataList,
                              adjustVars,
                              t0,
                              SL.ctime = NULL,
                              glm.ctime = NULL,
                              glm.family,
                              cvControl,
                              returnModels = FALSE,
                              verbose = TRUE,
                              gtol = 1e-3,
                              stratify = FALSE,
                              trtOfInterest,
                              ...) {
  if(length(trtOfInterest) > 1 & stratify){
    stratify <- FALSE
    warning("stratify option only supported if there is only a single trtOfInterest.")
  }
  if(stratify){
    stratify_include <- dataList[[1]]$trt == trtOfInterest[1]
  }else{
    stratify_include <- rep(TRUE, length(dataList[[1]][,1]))
  }
  include <- !(dataList[[1]]$t == dataList[[1]]$ftime & dataList[[1]]$C != 1 &
    dataList[[1]]$t < t0) & !(dataList[[1]]$t == dataList[[1]]$ftime & 
    dataList[[1]]$C == 1 & dataList[[1]]$t == t0) & 
    stratify_include

  ## determine whether to use linear or logistic regression in GLM fit
  if (!is.null(glm.family)) {
    glm_family <- parse(text = paste0("stats::", glm.family, "()"))
  }

  # if no SL library is specified, the code defaults to the specific GLM form
  if (is.null(SL.ctime)) {
    if (!(any(c("glm", "speedglm") %in% class(glm.ctime)))) {
      if (!all(dataList[[1]]$C == 0)) {
        ctimeForm <- stats::as.formula(sprintf("%s ~ %s", "C", glm.ctime))
        ctimeMod <- fast_glm(
          reg_form = ctimeForm,
          data = dataList[[1]][include, ],
          family = eval(glm_family)
        )
        if (unique(class(ctimeMod) %in% c("glm", "lm"))) {
          ctimeMod <- cleanglm(ctimeMod)
        }
      } else {
        dataList <- lapply(dataList, function(x) {
          x$G_dC <- 1
          x
        })
        ctimeMod <- "No censoring observed"
        class(ctimeMod) <- "noCens"
      }
    } else {
      ctimeMod <- glm.ctime
    }
    # as long as there are some observed censoring events,
    # get predictions from ctimeMod
    if (all(class(ctimeMod) != "noCens")) {
      dataList <- lapply(dataList, function(x) {
        g_dC <- rep(1, length(x[, 1]))
        if (t0 != 1) {
          # temporarily replace time with t-1
          # NOTE: this will fail if t enters model as a factor
          x$t <- x$t - 1

          suppressWarnings(
            g_dC <- 1 - predict(ctimeMod, newdata = x, type = "response")
          )

          # put time back to normal
          x$t <- x$t + 1
          # replace any observations with t = 1
          # to avoid extrapolation
          g_dC[x$t == 1] <- 1
        }
        x$G_dC <- as.numeric(unlist(by(g_dC, x$id, FUN = cumprod)))
        x
      })
      # if no observed censoring events, everybody gets 1
    } else {
      dataList <- lapply(dataList, function(x) {
        x$G_dC <- 1
        x
      })
    }
  } else {
    if (class(SL.ctime) != "SuperLearner") {
      if (!all(dataList[[1]]$C == 0)) {
        predictor_variables <- c(
            "t",
            names(adjustVars)
          )
        if(!stratify){
          predictor_variables <- c(
            predictor_variables, "trt"
          )
        }
        ctimeMod <- SuperLearner::SuperLearner(
          Y = dataList[[1]]$C[include],
          X = dataList[[1]][include, predictor_variables, drop = FALSE],
          id = dataList[[1]]$id[include],
          family = "binomial",
          SL.library = SL.ctime,
          verbose = verbose,
          cvControl = cvControl
        )
      } else {
        dataList <- lapply(dataList, function(x) {
          x$G_dC <- 1
          x
        })
        ctimeMod <- "No censoring observed"
        class(ctimeMod) <- "noCens"
      }
    } else { # if input SLlibrary.time is Super Learner object, just use that
      ctimeMod <- SL.ctime
    }
    if (class(ctimeMod) != "noCens") {
      predictor_variables <- c(
            "t",
            names(adjustVars)
          )
      if(!stratify){
        predictor_variables <- c(
          predictor_variables, "trt"
        )
      }
      dataList <- lapply(dataList, function(x) {
        g_dC <- rep(1, nrow(x))
        if (t0 != 1) {
          # temporarily replace time with t-1
          # NOTE: this will fail if t enters model as a factor
          x$t <- x$t - 1
          g_dC <-
            suppressWarnings(
              1 - predict(
                ctimeMod,
                newdata = x[, predictor_variables, drop = FALSE],
                onlySL = TRUE
              )[[1]]
            )

          # put time back to normal
          x$t <- x$t + 1
          # replace any observations with t = 1
          # to avoid extrapolation at t = 0
          g_dC[x$t == 1] <- 1
        }
        x$G_dC <- as.numeric(unlist(by(g_dC, x$id, FUN = cumprod)))
        x
      })
    } else {
      dataList <- lapply(dataList, function(x) {
        x$G_dC <- 1
        x
      })
    }
  }
  # truncate small propensities at gtol
  dataList <- lapply(dataList, function(x) {
    x$G_dC[x$G_dC < gtol] <- gtol
    x
  })

  out <- list(
    dataList = dataList,
    ctimeMod = if (returnModels) {
      ctimeMod
    } else {
      NULL
    }
  )
  return(out)
}
