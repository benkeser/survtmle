#' Super Learner prediction of hazards
#'
#' @param ftime A numeric vector of failure times. Right-censored observations
#'        should have corresponding \code{ftype} set to 0.
#' @param ftype A numeric vector indicating the type of failure. Observations
#'        with \code{ftype=0} are treated as a right-censored observation. Each
#'        unique value besides zero is treated as a separate type of failure.
#' @param adjustVars A data.frame of adjustment variables that will be used in
#'        estimating the conditional treatment, censoring, and failure (hazard
#'        or conditional mean) probabilities.
#' @param t0 The time at which to return cumulative incidence estimates. By
#'        default this is set to \code{max(ftime[ftype > 0])}.
#' @param SL.ftime A character vector or list specification to be passed to the
#'        \code{SL.library} option in the call to \code{SuperLearner} for the
#'        cause-specific hazards. See \code{?SuperLearner} for more information
#'        on how to specify valid \code{SuperLearner} libraries. It is expected
#'        that the wrappers used in the library will play nicely with the input
#'        variables, which will be called \code{"trt"},
#'        \code{names(adjustVars)}, and \code{"t"} if \code{method = "hazard"}.
#' @param SL.ctime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        estimate of the conditional hazard for censoring. It is expected that
#'        the wrappers used in the library will play nicely with the input
#'        variables, which will be called \code{"trt"} and
#'        \code{names(adjustVars)}.
#' @param glm.ctime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the estimate of the conditional hazard for censoring. Ignored if
#'        \code{SL.ctime} is not equal to \code{NULL}. Use \code{"trt"} to
#'        specify the treatment in this formula (see examples). The formula can
#'        additionally include any variables found in \code{names(adjustVars)}.

survSuperLearner <- function(ftime, ftype, adjustVars, 
                             t0, SL.ftime, SL.ctime, glm.ctime, ...){
  # assemble data frame of necessary variables
  n <- length(ftime)
  id <- seq_len(n)
  dat <- data.frame(id = id, ftime = ftime, ftype = ftype, trt = 1)
  if (is.null(adjustVars)) stop("adjustVars needed for SuperLearner prediction")

  nJ <- length(ftypeOfInterest)
  allJ <- sort(unique(ftype[ftype != 0]))

  # make long version of data sets needed for estimation and prediction
  longData <- makeDataList(
    dat = dat, J = allJ, ntrt = 1, uniqtrt = 1, t0 = t0, bounds = NULL
  )[[1]]

  # get rid of arbitrary trt column
  longData <- longData[,-which(colnames(longData) == "trt")]

  # estimate censoring distribution
  # TO DO: see how this handles when no trt is input
  censOut <- estimateCensoring(
    dataList = dataList,
    ntrt = 0,
    uniqtrt = NULL,
    t0 = t0,
    verbose = verbose,
    adjustVars = adjustVars,
    SL.ctime = SL.ctime,
    glm.ctime = glm.ctime,
    glm.family = glm.family,
    returnModels = returnModels,
    gtol = gtol
  )


}