#' Evaluate Results over Time Points of Interest
#'
#' Wrapper function for \code{survtmle} that takes a fitted \code{survtmle}
#' object and computes the TMLE estimated incidence for all times specified in
#' the \code{times} argument. For this function to work, the original call to
#' \code{survtmle} should have been executed with \code{returnModels = TRUE}.
#' This allows the function to be more efficient than repeated calls to
#' \code{survtmle} in that \code{timepoints} will use fitted censoring (and
#' hazard if \code{method="hazard"} was used in the original call) estimates. It
#' is therefore advisable that the vector \code{times} used in the call to
#' \code{timepoints} not include times beyond the time specified in \code{t0} in
#' the original call to \code{survtmle}. This can be ensured be making the
#' original call to \code{survtmle} with \code{t0 = max(ftime)}.
#'
#' @param object A fitted \code{survtmle} object with \code{returnModels = TRUE}
#' @param times The times to evaluate incidence.
#' @param returnModels Should the function return fitted GLM or Super Learner
#'        models at each timepoint. If set to \code{TRUE}, memory issues could
#'        arise depending on the number of timepoints specified and the size of
#'        the Super Learner library.
#'
#' @return An object of class \code{tp.survtmle} with number of entries equal to
#'         \code{length(times)}. Each entry is named "tX", where X denotes a
#'         single value of \code{times}.
#'
#' @export
#'
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' ftime <- round(runif(n, 1, 4))
#' ftype <- round(runif(n, 0, 2))
#' trt <- rbinom(n, 1, 0.5)
#' adjustVars <- data.frame(W1 = rnorm(n), W2 = rnorm(n))
#'
#' # fit an initial survtmle object with t0=max(ftime)
#' fm <- survtmle(ftime = ftime, ftype = ftype,
#'                trt = trt, adjustVars = adjustVars,
#'                glm.trt = "1", glm.ftime = "trt + W1 + W2",
#'                glm.ctime = "trt + W1 + W2", method="mean",
#'                returnModels = TRUE)
#'
#' # call timepoints to get cumulative incidence estimates at each timepoint
#' allTimes <- timepoints(object = fm, times = 1:4, returnModels = FALSE)
#'
#' # look at results for time 1
#' class(allTimes$t1)
#' allTimes$t1
#' # look at results for time 2
#' allTimes$t2
#'

timepoints <- function(object, times, returnModels = FALSE) {
  if(is.null(object$trtMod))
    stop("object must have returnModels = TRUE")
 
  callList <- as.list(object$call)[-1]
  cglm <- any(class(object$ctimeMod) %in% c("glm", "speedglm")) |
    any(class(object$ctimeMod) == "noCens")

  tglm <- any(class(object$trtMod) %in% c("glm", "speedglm"))
  ftglm <- ifelse(callList$method == "hazard",
                  any(class(object$ftimeMod[[1]]) %in% c("glm",
                                                         "speedglm")), FALSE)

  myOpts <- c("t0", "returnModels",
              ifelse(cglm, "glm.ctime", "SL.ctime"),
              ifelse(tglm, "glm.trt", "SL.trt"))
  if(callList$method == "hazard") {
    myOpts <- c(myOpts, ifelse(ftglm, "glm.ftime", "SL.ftime"))
  }
  funOpts <- callList[-which(names(callList) %in% myOpts)]

  funOpts$returnModels <- returnModels
  # used glm for censoring?
  if(cglm) {
    funOpts$glm.ctime <- object$ctimeMod
    funOpts$SL.ctime <- NULL
  } else {
    funOpts$SL.ctime <- object$ctimeMod
  }
  # used glm for trt?
  if(tglm) {
    funOpts$glm.trt <- object$trtMod
  } else {
    funOpts$SL.trt <- object$trtMod
  }
  # used glm for ftime
  if(ftglm & callList$method == "hazard") {
    funOpts$glm.ftime <- object$ftimeMod
  } else if(!ftglm & callList$method == "hazard") {
    funOpts$SL.ftime <- object$ftimeMod
  }
  # add in failure times, types, trt, and adjust
  funOpts$ftime <- object$ftime
  funOpts$ftype <- object$ftype
  funOpts$trt <- object$trt
  funOpts$adjustVars <- object$adjustVars

  outList <- vector(mode = "list", length = length(times))
  ct <- 0
  for(i in times) {
    ct <- ct + 1
    funOpts$t0 <- i
    if(all(object$ftime[object$ftype > 0] > i)) {
      outList[[ct]] <- list(est = rep(0, length(object$est)),
                            var = matrix(NA, nrow = length(object$est),
                                         ncol = length(object$est)))
    } else {
      if(i != object$t0) {
        outList[[ct]] <- do.call("survtmle", args = funOpts)
      } else {
        outList[[ct]] <- object
      }
    }
  }
  names(outList) <- paste0("t", times)
  class(outList) <- "tp.survtmle"
  return(outList)
}
