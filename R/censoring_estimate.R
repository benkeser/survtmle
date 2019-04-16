#' Estimate Censoring Mechanisms
#'
#' Computes an estimate of the hazard for censoring using either \code{glm} or
#' \code{SuperLearner} based on log-likelihood loss. The function then computes
#' the censoring survival distribution based on these estimates. The structure
#' of the function is specific to how it is called within \code{survtmle}. In
#' particular, \code{dataList} must have a very specific structure for this
#' function to run properly. The list should consist of \code{data.frame}
#' objects. The first will have the number of rows for each observation equal to
#' the \code{ftime} corresponding to that observation. The subsequent entries
#' will have \code{t0} rows for each observation and will set \code{trt} column
#' equal to each value of \code{trtOfInterest} in turn. One of these columns
#' must be named \code{C} that is a counting process for the right-censoring
#' variable. The function will fit a regression with \code{C} as the outcome and
#' functions of \code{trt} and \code{names(adjustVars)} as specified by
#' \code{glm.ctime} or \code{SL.ctime} as predictors.
#'
#' @param dataList A list of \code{data.frame} objects as described in
#'        \code{?makeDataList}.
#' @param adjustVars Object of class \code{data.frame} that contains the
#'        variables to adjust for in the regression.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#'        Needed only because the naming convention for the regression if
#'        \code{t == t0} is different than if \code{t != t0}.
#' @param SL.ctime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        outcome regression (either cause-specific hazards or conditional mean).
#'        See \code{?SuperLearner} for more information on how to specify valid
#'        \code{SuperLearner} libraries. It is expected that the wrappers used
#'        in the library will play nicely with the input variables, which will
#'        be called \code{"trt"} and \code{names(adjustVars)}.
#' @param glm.ctime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the outcome regression (either cause-specific hazards or
#'        conditional mean). Ignored if \code{SL.ctime != NULL}. Use \code{"trt"}
#'        to specify the treatment in this formula (see examples). The formula
#'        can additionally include any variables found in
#'        \code{names(adjustVars)}.
#' @param glm.family The type of regression to be performed if fitting GLMs in
#'        the estimation and fluctuation procedures. The default is "binomial"
#'        for logistic regression. Only change this from the default if there
#'        are justifications that are well understood. This is inherited from
#'        the calling function (either \code{mean_tmle} or \code{hazard_tmle}).
#' @param returnModels A boolean indicating whether to return the
#'        \code{SuperLearner} or \code{glm} objects used to estimate the
#'        nuisance parameters. Must be set to \code{TRUE} if the user plans to
#'        use calls to \code{timepoints} to obtain estimates at times other than
#'        \code{t0}. See \code{?timepoints} for more information.
#' @param gtol The truncation level of predicted censoring survival to handle
#'        positivity violations.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats as.formula predict glm
#' @importFrom SuperLearner SuperLearner
#'
#' @return The function returns a list that is exactly the same as the input
#'         \code{dataList}, but with a column named \code{G_dC} added to it,
#'         which is the estimated conditional survival distribution for the
#'         censoring variable evaluated at the each of the rows of each
#'         \code{data.frame} in \code{dataList}.
#' @export

estimateCensoring <- function(dataList,
                              adjustVars,
                              t0,
                              SL.ctime = NULL,
                              glm.ctime = NULL,
                              glm.family,
                              returnModels = FALSE,
                              verbose = TRUE,
                              gtol = 1e-3,
                              ...) {
  include <- !(dataList[[1]]$t == dataList[[1]]$ftime & dataList[[1]]$C != 1 &
    dataList[[1]]$t < t0) & !(dataList[[1]]$t == dataList[[1]]$ftime
  & dataList[[1]]$C == 1 & dataList[[1]]$t == t0)

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
        ctimeMod <- SuperLearner::SuperLearner(
          Y = dataList[[1]]$C[include],
          X = dataList[[1]][include, c(
            "t", "trt",
            names(adjustVars)
          )],
          id = dataList[[1]]$id[include],
          family = "binomial",
          SL.library = SL.ctime,
          verbose = verbose
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
                newdata = x[, c(
                  "t", "trt",
                  names(adjustVars)
                )],
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
