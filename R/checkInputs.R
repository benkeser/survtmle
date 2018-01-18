#' Check Function Inputs
#'
#' Check the input values of function parameters for errors.
#'
#' @param ftime A numeric vector of failure times. Right-censored observations
#'        should have corresponding \code{ftype} set to 0.
#' @param ftype A numeric vector indicating the type of failure. Observations
#'        with \code{ftype == 0} are treated as right-censored. Each unique
#'        value besides zero is treated as a separate type of failure.
#' @param trt A numeric vector indicating observed treatment assignment. Each
#'        unique value will be treated as an (unordered) separate type of
#'        treatment. Currently, only two unique values of \code{trt} are
#'        supported.
#' @param adjustVars A \code{data.frame} of adjustment variables that will be
#'        used in estimating the conditional treatment, censoring, and failure
#'        (hazard or conditional mean) probabilities.
#' @param t0 The time at which to return cumulative incidence estimates. By
#'        default this is set to \code{max(ftime)}.
#' @param SL.ftime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        outcome regression (either cause-specific hazards or conditional mean).
#'        See \code{?SuperLearner} for more information on how to specify valid 
#'        \code{SuperLearner} libraries. It is expected that the wrappers used
#'        in the library will play nicely with the input variables, which will
#'        be called \code{"trt"} and \code{names(adjustVars)}.
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
#'        for the outcome regression (either cause-specific hazards or
#'        conditional mean). Ignored if \code{SL.ftime != NULL}. Use \code{"trt"}
#'        to specify the treatment in this formula (see examples). The formula
#'        can additionally include any variables found in
#'        \code{names(adjustVars)}.
#' @param glm.ctime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the estimate of the conditional hazard for censoring. Ignored if
#'        \code{SL.ctime != NULL}. Use \code{"trt"} to specify the treatment in
#'        this formula (see examples). The formula can additionally include any
#'        variables found in \code{names(adjustVars)}.
#' @param glm.trt A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the estimate of the conditional probability of treatment. Ignored
#'        if \code{SL.trt != NULL}. By default set to "1", corresponding to
#'        using empirical estimates of each value of \code{trt}. The formula can
#'        include any variables found in \code{names(adjustVars)}.
#' @param returnIC A boolean indicating whether to return vectors of influence
#'        curve estimates. These are needed for some post-hoc comparisons, so it
#'        is recommended to leave as \code{TRUE} (the default) unless the user
#'        is sure these estimates will not be needed later.
#' @param returnModels A boolean indicating whether to return the
#'        \code{SuperLearner} or \code{glm} objects used to estimate the
#'        nuisance parameters. Must be set to \code{TRUE} if the user plans to
#'        use calls to \code{timepoints} to obtain estimates at times other
#'        than \code{t0}. See \code{?timepoints} for more information.
#' @param ftypeOfInterest An input specifying what failure types to compute
#'        estimates of incidence for. The default value computes estimates for
#'        values \code{unique(ftype)}. Can alternatively be set to a vector of
#'        values found in \code{ftype}.
#' @param trtOfInterest An input specifying which levels of \code{trt} are of
#'        interest. The default value computes estimates for values
#'        \code{unique(trt)}. Can alternatively be set to a vector of values
#'        found in \code{trt}.
#' @param method A character specification of how the targeted minimum
#'       loss-based estimators should be computer, either \code{"mean"} or
#'       \code{"hazard"}. The \code{"mean"} specification uses a closed-form
#'       targeted minimum loss-based estimation based on the G-computation
#'       formula of Bang and Robins (2005). The \code{"hazard"} specification
#'       uses an iteratively algorithm based on cause-specific hazard functions.
#'       The latter specification has no guarantee of convergence in finite
#'       samples. The convergence can be influenced by the stopping criteria
#'       specified in the \code{tol}. Future versions may implement a closed
#'       form version of this hazard-based estimator.
#' @param bounds A list of bounds.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param tol The stopping criteria when \code{method = "hazard"}. The algorithm
#'        will continue performing targeting updates to the initial estimators
#'        until the empirical mean of the efficient influence function is
#'        smaller than \code{tol}. The default (\code{1/length(ftime)}) is a
#'        sensible value. Larger values can be used in situations where
#'        convergence of the algorithm is an issue; however, this may result in 
#'        large finite-sample bias.
#' @param Gcomp A boolean indicating whether to compute the G-computation
#'        estimator (i.e., a substitution estimator with no targeting step).
#'        Note, theory does not support inference for the G-computation
#'        estimator if Super Learner is used to estimate failure and censoring
#'        mechanisms. Only implemented for \code{method = "mean"}.
#' @param maxIter A maximum number of iterations for the algorithm when
#'        \code{method = "hazard"}. The algorithm will iterate until either the
#'        empirical mean of the efficient influence function is smaller than
#'        \code{tol} or until \code{maxIter} iterations have been completed.
#'
#' @return Options to be passed to \code{mean_tmle} or \code{hazard_tmle}.
#'

checkInputs <- function(ftime,
	                ftype,
	                trt,
                        adjustVars,
                        t0 = max(ftime[ftype > 0]),
                        SL.ftime = NULL,
                        SL.ctime = NULL,
                        SL.trt = NULL,
                        glm.ftime = NULL,
                        glm.ctime = NULL,
                        glm.trt = "1",
                        returnIC = TRUE,
                        returnModels = TRUE,
                        ftypeOfInterest = unique(ftype[ftype != 0]),
                        trtOfInterest=unique(trt),
                        method = "hazard",
                        bounds = NULL,
                        verbose = FALSE,
                        tol = 1 / (length(ftime)),
                        maxIter = 100,
                        Gcomp = FALSE) {

  # check for NULL values
  if(sum(is.null(ftime)) > 0 |
     sum(is.null(ftype)) > 0 |
     sum(is.null(trt)) > 0) {
    stop("NULL values in ftime, ftype, trt, and adjustVars not allowed")
  }

  # check for missing values
  if(sum(is.na(ftime)) > 0 | sum(is.na(ftype)) > 0 | sum(is.na(trt)) > 0) {
    stop("Missing values in ftime, ftype, trt, and adjustVars not supported.")
  }

  # check that trt is vector
  if(!is.vector(trt)) {
    stop("trt must be a vector.")
  }

  # check that not all failures are < t0
  if(all(ftime[ftype > 0] > t0)) {
    stop("No observed events by t0.")
  }

  # if(length(unique(trt)) > 2) {
  #   stop("trt with more than 2 unique values not yet supported.")
  # }

  # check for reserved names in columns of adjustVars
  if(any(colnames(adjustVars) == "t")) {
    stop("t is a reserved name. Please rename that column from adjustVars.")
  }

  if(!is.null(adjustVars)) {
    if(sum(is.na(adjustVars)) > 0) {
      stop("Missing values in ftime, ftype, trt, and adjustVars not supported.")
    }
  }

  # check for G-comp for hazard
  if(method == "hazard" & Gcomp) {
    warning("G-computation estimator not implemented for method='hazard'.
	    Proceeding with TMLE.")
  }

  # check for ftime with 0
  if(any(ftime <= 0)) {
    stop("Some failure times less than or equal zero. Remove these observations
	 and try again")
  }

  # check if no events in each trt/type combo
  for(j in ftypeOfInterest) {
    for(z in trtOfInterest) {
      if(sum(ftype == j & trt == z) == 0) {
        stop(paste0("No endpoints of type ftype = ", j, " in group trt = ",
	            z, ". Adjust ftypeOfInterest and trtOfInterest accordingly."
	           )
	    )
      }
    }
  }

  # haven't figured out how to fix updateVariables yet when ftypeOfInterest is
  # not equal to unique(trt).
  if(!all(unique(ftypeOfInterest) == unique(ftypeOfInterest)) &
     method == "hazard") {
    stop("Hazard implementation is not yet functional when ftypeOfInterest does
	 not include all unique values of trt")
  }

  # check that all trt of interest are observed
  if(!(all(trtOfInterest %in% trt))) {
    stop("At least one trtOfInterest not observed. Remove from trtOfInterest and
	 try again.")
  }

  # check if adjustVars is data.frame
  if(!is.data.frame(adjustVars) & !is.null(adjustVars)) {
    stop("adjustVars should be a data.frame or NULL.")
  }

  # check if method is known
  if(!(method %in% c("hazard", "mean"))) {
    stop("method should be either 'hazard' or 'mean'.")
  }

  # warn if tol is too large
  if(tol > 1 / sqrt(length(ftime))) {
    warning("tol is larger than 1/sqrt(n), consider decreasing.")
  }
 
  # stop if negative ftypes
  if(any(ftype < 0)) {
    stop("Some ftype less than 0. ftype should use 0 to denote censoring and
	 positive numeric values to denote other endpoints.")
  }

  # check if time enters as a factor in glm.ftime or glm.ctime
  if(!is.null(glm.ftime) & !is.list(glm.ftime) & method == "hazard") {
    if(grepl("factor(t)", glm.ftime)) {
      stop("Time can only be modeled as a factor in hazard implementation if
	   there are observed endpoints at every time 1:t0.")
    }
  }
  if(!is.null(glm.ctime) & !any(class(glm.ctime) %in% c("speedglm", "glm"))) {
    if(grepl("factor(t)", glm.ctime)) {
      stop("Time can only be modeled as a factor in hazard implementation if
	   there are observed endpoints at every time 1:t0.")
    }
  }

  # warn if no events at t0
  if(method == "hazard") {
    for(j in ftypeOfInterest) {
      for(z in trtOfInterest) {
	if(t0 > max(ftime[ftype == j & trt == z])) {
	  warning(paste0("t0 larger than last observed endpoint of ftype = ", j,
	                 " and trt = ", z, ". Hazard TMLE may be extrapolating
			 to estimate incidence."))
	}
      }
    }
  }

  if(t0 > max(ftime[ftype > 0]) & method == "mean") {
    warning("t0 larger than last observed endpoint. Mean-based TMLE assumes
	    constant incidence between last observed failure type and this
	    time."
	   )
  }

  # check if both glm and SL are specified
  if(!(is.null(glm.trt)) & !(is.null(SL.trt))) {
    warning("glm.trt and SL.trt specified. Proceeding with SL.trt")
    glm.trt <- NULL
  }
  if(!(is.null(glm.ftime)) & !(is.null(SL.ftime))) {
    warning("glm.ftime and SL.ftime specified. Proceeding with SL.ftime")
    glm.ftime <- NULL
  }
  if(!(is.null(glm.ctime)) & !(is.null(SL.ctime))) {
    warning("glm.ctime and SL.ctime specified. Proceeding with SL.ctime")
    glm.ctime <- NULL
  }
  # check glm formulas
  if(!is.null(glm.trt)) {
    tryCatch({
      tmp <- as.formula(paste0("trt ~", glm.trt))
    }, error = function(e) {
      stop("glm.trt formula appears to be invalid.")
    })
  }
  if(!is.null(glm.ctime)) {
    if(all(glm.ctime != "No censoring observed")) {
      tryCatch({
	tmp <- as.formula(paste0("C ~", glm.ctime))
      }, error = function(e) {
	stop("glm.ctime formula appears to be invalid.")
      })
    }
  }
  if(!is.null(glm.ftime)) {
    tryCatch({
      tmp <- as.formula(paste0("N ~", glm.ftime))
    }, error = function(e) {
      stop("glm.ftime formula appears to be invalid.")
    })
  }

  # check that one of glm.trt or SL.trt is specified
  if(is.null(glm.trt) & is.null(SL.trt)) {
    warning("glm.trt and SL.trt not specified. Proceeding with glm.trt = '1'")
    glm.trt <- "1"
  }
  # check that one of glm.ftime or SL.ftime is specified
  if(is.null(glm.ftime) & is.null(SL.ftime)) {
    warning("glm.ftime and SL.ftime not specified. Computing empirical estimates.")
    if(method == "hazard") {
      glm.ftime <- paste0("-1 + ", paste0("I(t == ", unique(ftime[ftype > 0]),
					  ")", collapse = "+"), "+",
			  paste0("I(trt*t == ", unique(ftime[ftype > 0]), ")",
				 collapse = "+"))
    } else {
      glm.ftime <- "trt"
    }
  }
  # check that one of glm.ctime or SL.ctime is specified
  if(is.null(glm.ctime) & is.null(SL.ctime)) {
    warning("glm.ctime and SL.ctime not specified. Computing Kaplan-Meier estimates.")
    glm.ctime <- paste0("-1 + ", paste0("I(t == ", unique(ftime[ftype == 0]),
					")", collapse = "+"), "+",
			paste0("I(trt*t == ", unique(ftime[ftype == 0]), ")",
			       collapse="+"))
  }

  # if covariates are null
  if(is.null(adjustVars)) {
    warning("adjustVars = NULL. Computing unadjusted estimates.")
    if(is.null(glm.trt)) {
      glm.trt <- "1"
    }
    if(is.null(glm.ctime)) {
      glm.ctime <- paste0("-1 + ", paste0("I(t == ", unique(ftime[ftype == 0]),
					  ")", collapse = "+"), "+",
			  paste0("I(trt*t == ", unique(ftime[ftype == 0]), ")",
				 collapse = "+"))
    }
    if(is.null(glm.ftime)) {
      glm.ftime <- paste0("-1 + ", paste0("I(t==", unique(ftime[ftype > 0]),
					  ")",collapse = "+"), "+",
			  paste0("I(trt*t == ", unique(ftime[ftype > 0]), ")",
				 collapse = "+"))
    }
    SL.trt <- SL.ctime <- SL.ftime <- NULL
    # add in dummy adjustVars so nothing else complains about NULL
    adjustVars <- data.frame(dummy = rep(1,length(ftime)))
    } else {
    if(all(apply(adjustVars, 2, function(x){length(unique(x)) == 1}))) {
      warning("Columns of adjustVars are constantly valued. Computing unadjusted estimates.")
      if(is.null(glm.trt)) {
	glm.trt <- "1"
      }
      if(is.null(glm.ctime)) {
	glm.ctime <- paste0("-1 + ", paste0("I(t == ",
					    unique(ftime[ftype == 0]), ")",
					    collapse = "+"), "+",
			    paste0("I(trt*t == ", unique(ftime[ftype == 0]),
				   ")", collapse = "+"))
      }
      if(is.null(glm.ftime)) {
	glm.ftime <- paste0("-1 + ", paste0("I(t == ", unique(ftime[ftype > 0]),
					    ")", collapse = "+"), "+",
			    paste0("I(trt*t == ", unique(ftime[ftype > 0]), ")",
				   collapse = "+"))
      }
      SL.trt <- SL.ctime <- SL.ftime <- NULL
      # add in dummy adjustVars so nothing else complains about NULL
    }
  }

  # check format of bound inputs
  if(!(is.null(bounds))) {
    if(!any(colnames(bounds) == "t")) {
      stop("bounds should have a column named 't' containing a row for each
	   value 1:t0.")
    }
    # if(!(all(1:t0 %in% bounds$t))){
    # stop("bounds should have one row for each time 1:t0.")
    # }
  }

  # return clean variables
  return(list(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	      glm.trt = glm.trt, SL.trt = SL.trt, glm.ftime = glm.ftime,
	      glm.ctime = glm.ctime, SL.ftime = SL.ftime, SL.ctime = SL.ctime))
}
