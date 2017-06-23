#' survtmle
#' 
#' This function estimates the marginal cumulative incidence for failures of
#' specified types using targeted minimum loss-based estimation. More description to come...
#' 
#' @param ftime A numeric vector of failure times. Missing values are not supported. Right 
#' censored observations should have a non-missing value for \code{ftime} while 
#' \code{ftype} should be set to 0. 
#' @param ftype A numeric vector indicating the type of failure with type 0 reserved for 
#' right censored observations. Each unique value will be treated as an (unordered) separate
#' type of failure.
#' @param trt A numeric vector indicating observed treatment assignment. Each unique value will 
#' be treated as an (unordered) separate type of treatment. Currently, only two unique values of 
#' \code{trt} are supported.
#' @param adjustVars A data.frame of adjustment variables that will be used in estimating the 
#' conditional treatment, censoring, and failure (hazard or conditional mean) probabilities. 
#' @param t0 The time at which to return cumulative incidence estimates. By default this is set
#' to \code{max(ftime)}.
#' @param incidence If \code{TRUE} the function return estimates of cumulative incidence. If 
#' \code{FALSE} the function returns esimtates of survival. 
#' @param SL.ftime A character vector or list specification to be passed to the \code{SL.library} argument 
#' in the call to \code{SuperLearner} for the outcome regression (either cause-specific hazards or 
#' condtional mean). See \code{?SuperLearner} for more information on how to specify valid 
#' \code{SuperLearner} libraries. It is expected that the wrappers used in the library will play nicely
#' with the input variables, which will be called \code{"trt"} and \code{names(adjustVars)}. 
#' @param SL.ctime A character vector or list specification to be passed to the \code{SL.library} argument 
#' in the call to \code{SuperLearner} for the estimate of the conditional hazard for censoring. It is expected 
#' that the wrappers used in the library will play nicely with the input variables, which will be 
#' called \code{"trt"} and \code{names(adjustVars)}. 
#' @param SL.trt A character vector or list specification to be passed to the \code{SL.library} argument 
#' in the call to \code{SuperLearner} for the estimate of the conditional probabiltiy of treatment.
#' It is expected that the wrappers used in the library will play nicely
#' with the input variables, which will be \code{names(adjustVars)}. 
#' @param glm.ftime A character specification of the right-hand side of the equation passed to the
#' \code{formula} option of a call to \code{glm} for the outcome regression (either cause-specific hazards or 
#' condtional mean). Ignored if \code{SL.ftime != NULL}. Use \code{"trt"} to specify the treatment 
#' in this formula (see examples). The formula can additionally include any variables found in 
#' \code{names(adjustVars)}. 
#' @param glm.ctime A character specification of the right-hand side of the equation passed to the
#' \code{formula} option of a call to \code{glm} for the estimate of the conditional hazard for 
#' censoring. Ignored if \code{SL.ctime != NULL}. Use \code{"trt"} to specify the treatment 
#' in this formula (see examples). The formula can additionally include any variables found in 
#' \code{names(adjustVars)}. 
#' @param glm.trt A character specification of the right-hand side of the equation passed to the
#' \code{formula} option of a call to \code{glm} for the estimate of the conditional probabiltiy of 
#' treatment. Ignored if \code{SL.trt != NULL}. By default set to "1", corresponding to using empirical
#' estimates of each value of \code{trt}. The formula can include any variables found in 
#' \code{names(adjustVars)}. 
#' @param returnIC A boolean indicating whether to return vectors of influence curve estimates. These are
#' needed for some post-hoc comarisons, so it is recommended to leave as \code{TRUE} (the default)
#' unless the user is sure these estimates will not be needed later. 
#' @param returnModels A boolean indicating whether to return the \code{SuperLearner} or \code{glm} 
#' objects used to estimate the nuisance parameters. Must be set to \code{TRUE} if the user plans to 
#' use calls to \code{timepoints} to obtain estimates at times other than \code{t0}. See \code{?timepoints}
#' for more information. 
#' @param ftypeOfInterest An input specifying what failure types to compute estimates of incidence for. 
#' The default value is \code{"all"}, which computes estimates for values \code{unique(ftype)}. Can alternatively
#' be set to a vector of values found in \code{ftype}.
#' @param trtOfInterest An input specifying which levels of \code{trt} are of interest. The default value
#' is \code{"all"}, which computes estimates for values \code{unique(trt)}. Can alternatively be set to a 
#' vector of values found in \code{trt}.
#' @param method A character specification of how the targeted minimum loss-based estimators should be 
#' computer, either \code{"mean"} or \code{"hazard"}. The \code{"mean"} specification uses a closed-form
#' targeted minimum loss-based estimation based on the G-computation formula of Bang and Robins (2005).
#' The \code{"hazard"} specification uses an iteratively algorithm based on cause-specific hazard functions.
#' The latter specification has no guarantee of convergence in finite samples. The convergence can be 
#' influenced by the stopping criteria specified in the \code{tol}. Future versions may implement a closed 
#' form version of this hazard-based estimator. 
#' @param bounds A list of bounds... XXX NEED MORE DESCRIPTION HERE XXX
#' @param verbose A boolean indicating whether the function should print messages to indicate progress.
#' @param tol The stopping criteria when \code{method="hazard"}. The algorithm will continue performing
#' targeting updates to the initial estimators until the empirical mean of the efficient influence function
#' is smaller than \code{tol}. The default (\code{1/length(ftime)}) is a sensible value. Larger values can
#' be used in situations where convergence of the algorithm is an issue; however, this may result in 
#' large finite-sample bias. 
#' @param Gcomp A boolean indicating whether to compute the G-computation estimator (i.e., a substitution
#' estimator with no targeting step). Note, theory does not support inference for the Gcomp estimator if using 
#' super learner is used to estimate failure and censoring mechanisms. Only implemented if \code{method="mean"}.
#' @param maxIter A maximum number of iterations for the algorithm when \code{method="hazard"}. The 
#' algorithm will iterate until either the empirical mean of the efficient influence function
#' is smaller than \code{tol} or until \code{maxIter} iterations have been completed. 
#' 
#' 
#' @return An object of class \code{survtmle}.
#' @return call The call to \code{survtmle}.
#' @return est A numeric vector of point estimates -- one for each combination of \code{ftypeOfInterest}
#' and \code{trtOfInterest}.
#' @return var A covariance matrix of all the point estimates
#' @return meanIC The empirical mean of the efficient influence function at the estimated, fluctuated
#' nuisance parameters. If all goes well, each value should be small. If
#' \code{method="hazard"} then each value should be less than \code{tol}. If code{method="mean"}
#' then should be on the order of \code{1e-7}. 
#' @return ic The efficient influence function at the estimated, fluctuated nuisance parameters,
#' evaluated on each of the observations (summed over all times). These may be used to perform
#' post-hoc comparisons. More details coming. 
#' @return ftimeMod If \code{returnModels=TRUE} the fit object(s) for the call to \code{glm} or 
#' \code{SuperLearner} for the outcome regression models. If \code{method="mean"} this will be a list
#' of length \code{length(ftypeOfInterest)} each of length \code{t0} (one regression for each
#' failure type and for each time point). If \code{method="hazard"} this will be a list of length
#' \code{length(ftypeOfInterest)} with one model corresponding with the (pseudo-)hazard for each
#' cause of failure. If \code{returnModels=FALSE}, this will equal \code{NULL}.
#' @return ctimeMod If \code{returnModels=TRUE} the fit object for the call to \code{glm} or 
#' \code{SuperLearner} for the censoring hazard regression model.  
#' If \code{returnModels=FALSE}, this will equal \code{NULL}.
#' @return trtMod If \code{returnModels=TRUE} the fit object for the call to \code{glm} or 
#' \code{SuperLearner} for the conditioanl probability of \code{trt} regression model. 
#' If \code{returnModels=FALSE}, this will equal \code{NULL}.
#' @return t0 The timepoint the function was evaluated at. 
#' 
#' @export
#' 
#' @examples
#' 
#' ## Single failure type examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' trt <- rbinom(n,1,0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
#' 
#' ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n,0,1))
#' 
#' # Fit 1
#' # fit a survtmle object with glm estimators for treatment, censoring, and failure
#' # using the "mean" method
#' fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", returnModels = TRUE)
#' # fit1
#' 
#' # Fit 2
#' # fit an survtmle object with SuperLearner estimators for failure and censoring
#' # and empirical estimators for treatment (default) using the "mean" method
#' fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="mean", returnModels = TRUE)
#' # fit2
#' 
#' # Fit 3
#' # repeat Fit 1 using the "hazard" method 
#' fit3 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="hazard", returnModels = TRUE)
#' # fit3
#' 
#' # Fit 4
#' # repeat Fit 2 using the "hazard" method
#' fit4 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="hazard", returnModels = TRUE)
#' # fit4
#' 
#' ## Multiple failure type examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' trt <- rbinom(n,1,0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
#' 
#' ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n,0,2))
#' 
#' # Fit 5
#' # fit a survtmle object with glm estimators for treatment, censoring, and failure
#' # using the "mean" method
#' fit5 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", returnModels = TRUE)
#' # fit5
#' 
#' # Fit 6
#' # fit an survtmle object with SuperLearner estimators for failure and censoring
#' # and empirical estimators for treatment (default) using the "mean" method
#' fit6 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="mean", returnModels = TRUE)
#' # fit6
#' 
#' # Fit 7
#' # repeat Fit 5 using the "hazard" method 
#' fit7 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="hazard", returnModels = TRUE)
#' # fit7
#' 
#' # Fit 8
#' # repeat Fit 6 using the "hazard" method
#' fit8 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="hazard", returnModels = TRUE)
#' # fit8 
#' 
#' # Fit 9
#' # repeat Fit 5, but only return estimates for ftype = 1
#' fit9 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", returnModels = TRUE, ftypeOfInterest = 1)
#' # fit9
#' 
#' # TODO: Add bounds example

survtmle <- function(
  ftime, 
  ftype,
  trt,
  adjustVars,
  t0=max(ftime[ftype > 0]),
  incidence=TRUE,
  SL.ftime=NULL,
  SL.ctime=NULL,
  SL.trt=NULL,
  glm.ftime=NULL,
  glm.ctime=NULL,
  glm.trt="1",
  returnIC=TRUE,
  returnModels=TRUE,
  ftypeOfInterest="all",
  trtOfInterest="all",
  method="hazard",
  bounds=NULL,
  verbose=FALSE,
  tol=1/(length(ftime)),
  maxIter=100,
  Gcomp = FALSE
){
  
  call <- match.call(expand.dots = TRUE)
  
  # check for missing values
  if(sum(is.na(ftime))>0 | sum(is.na(ftype))>0 | sum(is.na(trt))>0 | sum(is.na(adjustVars))>0){
    stop("Missing values in ftime, ftype, trt, or adjustVars not currently supported.")
  }
  
  # check for G-comp for hazard
  if(method=="hazard" & Gcomp){
    warning("G-computation estimator not implemented for method='hazard'. Proceeding with TMLE.")
  }

  # check for ftime with 0
  if(any(ftime<=0)){
    warning("Some failure times less than or equal zero. Dropping these observations")
    ind <- which(ftime>0)
    ftime <- ftime[ind]; ftype <- ftype[ind]
    adjustVars <- adjustVars[ind,,drop=FALSE],
    trt <- trt[ind]
  }
  # number of failure types
  nJ <- length(unique(ftype))-1
  #if(nJ >= 2) print(paste("ftype has ", nJ, " unique failure types. Calculating cumulative incidence estimates."))
  
  # hazard-based TMLE
  if(method=="hazard"){
    tmle.fit <- hazard_tmle(ftime=ftime, 
                            ftype=ftype,
                            trt=trt,
                            t0=t0,
                            incidence=incidence,
                            adjustVars=adjustVars,
                            SL.ftime=SL.ftime,
                            SL.ctime=SL.ctime,
                            SL.trt = SL.trt,
                            glm.ftime=glm.ftime,
                            glm.ctime=glm.ctime,
                            glm.trt=glm.trt,
                            returnIC=returnIC,
                            returnModels=returnModels,
                            ftypeOfInterest=ftypeOfInterest,
                            trtOfInterest=trtOfInterest,
                            bounds=bounds,
                            verbose=verbose, 
                            tol=tol, 
                            maxIter=maxIter)
  }else if(method=="mean"){
    tmle.fit <- mean_tmle(ftime=ftime, 
                          ftype=ftype,
                          trt=trt,
                          t0=t0,
                          incidence=incidence,
                          adjustVars=adjustVars,
                          SL.ftime=SL.ftime,
                          SL.ctime=SL.ctime,
                          SL.trt = SL.trt,
                          glm.ftime=glm.ftime,
                          glm.ctime=glm.ctime,
                          glm.trt=glm.trt,
                          returnIC=returnIC,
                          returnModels=returnModels,
                          ftypeOfInterest=ftypeOfInterest,
                          trtOfInterest=trtOfInterest,
                          bounds=bounds,
                          verbose=verbose, 
                          tol=tol,
                          Gcomp=Gcomp
    )
  }
  
  out <- list(
    call=call,
    est=tmle.fit$est, var=tmle.fit$var,
    meanIC=tmle.fit$meanIC, ic=tmle.fit$ic,
    ftimeMod=tmle.fit$ftimeMod, ctimeMod=tmle.fit$ctimeMod,
    trtMod=tmle.fit$trtMod, t0=t0
  )
  class(out) <- "survtmle"
  out
}
