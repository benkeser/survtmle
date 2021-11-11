#' Estimate Missing Failure Time mechanism
#'
#' @description This function computes the conditional probability of having
#'  non-missing \code{ftime} for each level of treatment using \code{\link[stats]{glm}}
#'  or \code{\link[SuperLearner]{SuperLearner}}. 
#'
#' @param dat An object of class \code{data.frame}. Must have named column
#'  \code{trt}.
#' @param adjustVars An object of class \code{data.frame} that will be used
#'  either as the \code{data} argument in a call to \code{\link[stats]{glm}} or
#'  as \code{X} in a call to \code{\link[SuperLearner]{SuperLearner}}.
#' @param glm.ftimeMissing A character specification of the right-hand side of the
#'  equation passed to the \code{\link[stats]{formula}} option of a call to
#'  \code{\link[stats]{glm}} for the missingness failure times model. Ignored if
#'  \code{SL.ftimeMissing} is not equal to \code{NULL}. Use \code{"trt"} to specify
#'  the treatment in this formula (see examples). The formula can additionally
#'  include any variables found in \code{names(adjustVars)}.
#' @param SL.ftimeMissing A character vector or list specification to be passed to the
#'  \code{SL.library} in the call to \code{\link[SuperLearner]{SuperLearner}}
#'  for the missingness model for failure time.
#'  It is expected that the wrappers used in the library will play nicely with
#'  the input variables, which will be called \code{"trt"} and 
#'  \code{names(adjustVars)}.
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
#' @param trtOfInterest An input specifying which levels of \code{trt} are of
#'  interest. The default value computes estimates for all of the values in
#'  \code{unique(trt)}. Can alternatively be set to a vector of values found in
#'  \code{trt}.
#' @param gtol The truncation level of predicted missingness probabilities to handle
#'  positivity violations.
#' @param ... Other arguments. Not currently used
#'
#' @return dat The input \code{data.frame} object with two added columns
#'  corresponding with the conditional probability (given \code{adjustVars}) of
#'  \code{trt==max(trt)} and \code{trt==min(trt)}.
#' @return ftimeMissingMod If \code{returnModels = TRUE}, the fitted \code{glm} or
#'  \code{SuperLearner} object. Otherwise, \code{NULL}
#'
#' @importFrom stats as.formula predict model.matrix optim glm
#' @importFrom SuperLearner SuperLearner SuperLearner.CV.control All SL.mean
#'  SL.glm SL.step
#'
#' @export
estimateFtimeMissing <- function(dat,
                              adjustVars,
                              glm.ftimeMissing = NULL,
                              SL.ftimeMissing = NULL,
                              cvControl,
                              returnModels = FALSE,
                              verbose = FALSE,
                              gtol = 1e-3,
                              trtOfInterest,
                              ...) {

  if(any(is.na(dat$ftime))){
  # binarize the outcome
  thisY <- as.numeric(!is.na(dat$ftime))

  # fit Super Learner
  if (!is.null(SL.ftimeMissing)) {
    if (class(SL.ftimeMissing) != "SuperLearner") {
      ftimeMissingMod <- SuperLearner::SuperLearner(
        Y = thisY, X = data.frame(trt = dat$trt, adjustVars),
        SL.library = SL.ftimeMissing,
        id = dat$id, verbose = verbose,
        family = "binomial",
        cvControl = cvControl
      )
    } else {
      ftimeMissingMod <- SL.ftimeMissing
    }
    # predict under each level of treatment!
    for(a in trtOfInterest){
      new_data <- data.frame(trt = a, adjustVars)
      dat[[paste0("g_", a)]] <- dat[[paste0("g_", a)]] * as.numeric(predict(
        ftimeMissingMod, newdata = new_data
      )[[1]])
    } 

  } else if (!is.null(glm.ftimeMissing) & is.null(SL.ftimeMissing)) {
    # set up model formula and data for the treatment regression
    ftimeMissing_form <- paste("thisY", "~", glm.ftimeMissing)
    fit_data <- data.frame(thisY = thisY, trt = dat$trt, adjustVars)
    ftimeMissingMod <- stats::glm(as.formula(ftimeMissing_form), 
                                  family = binomial(),
                                  data = fit_data)

    for(a in trtOfInterest){
      new_data <- data.frame(trt = a, adjustVars)
      dat[[paste0("g_", a)]] <- dat[[paste0("g_", a)]] * as.numeric(predict(
        ftimeMissingMod, newdata = new_data, type = "response"
      ))
    } 
  }

  # truncate propensities
  for(a in trtOfInterest){
    eval(parse(text = paste0(
      "dat$g_", a, "[dat$g_", a,
      "< gtol]<- gtol"
    )))
  }
  }else{
    ftimeMissingMod <- "No missing failure times"    
  }
  out <- list()
  out$dat <- dat
  out$ftimeMissingMod <- NULL
  if (returnModels) out$ftimeMissingMod <- ftimeMissingMod
  return(out)
}
