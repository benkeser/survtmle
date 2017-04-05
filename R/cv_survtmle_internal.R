#' cv_survtmle_internal
#' 
#' This function estimates the marginal cumulative incidence for failures of
#' specified types using cross-validated targeted minimum loss-based estimation. More description to come...
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
#' @param survtmleFits A list of objects of class \code{survtmle} fit on train samples of data
#' @param validFolds A \code{list} of \code{vectors} indicating which observations were in the 
#' validation fold for each of the \code{survtmleFits}. 
#' @param returnIC A boolean indicating whether to return vectors of influence curve estimates. These are
#' needed for some post-hoc comarisons, so it is recommended to leave as \code{TRUE} (the default)
#' unless the user is sure these estimates will not be needed later. 
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
#' @param maxIter A maximum number of iterations for the algorithm when \code{method="hazard"}. The 
#' algorithm will iterate until either the empirical mean of the efficient influence function
#' is smaller than \code{tol} or until \code{maxIter} iterations have been completed. 
#' @param returnModels A boolean indicating whether to return the \code{SuperLearner} or \code{glm} 
#' objects used to estimate the nuisance parameters. Must be set to \code{TRUE} if the user plans to 
#' use calls to \code{timepoints} to obtain estimates at times other than \code{t0}. See \code{?timepoints}
#' for more information. 
#' 
#' 
#' @return An object of class \code{cv_survtmle}.
#' @return call The call to \code{cv_survtmle}.
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

cv_survtmle_internal <- function(
  ftime, 
  ftype,
  trt,
  adjustVars,
  t0=max(ftime),
  incidence=TRUE,
  survtmleFits,
  validFolds,
  returnIC=TRUE,
  returnModels=TRUE,
  ftypeOfInterest="all",
  trtOfInterest="all",
  method="hazard",
  bounds=NULL,
  verbose=TRUE,
  tol=1/(length(ftime)),
  maxIter=100
){
  
  call <- match.call(expand.dots = TRUE)
  
  # check for missing values
  if(sum(is.na(ftime))>0 | sum(is.na(ftype))>0 | sum(is.na(trt))>0 | sum(is.na(adjustVars))>0){
    stop("Missing values in ftime, ftype, trt, or adjustVars not currently supported.")
  }
  
  # number of failure types
  nJ <- length(unique(ftype))-1
  #if(nJ >= 2) print(paste("ftype has ", nJ, " unique failure types. Calculating cumulative incidence estimates."))
  
  # hazard-based TMLE
  if(method=="hazard"){
    cv_tmle.fit <- cv_hazard.tmle(ftime=ftime, 
                            ftype=ftype,
                            trt=trt,
                            t0=t0,
                            incidence=incidence,
                            adjustVars=adjustVars,
                            validFolds = validFolds, 
                            survtmleFits = survtmleFits,
                            returnIC=returnIC,
                            returnModels=returnModels,
                            ftypeOfInterest=ftypeOfInterest,
                            trtOfInterest=trtOfInterest,
                            bounds=bounds,
                            verbose=verbose, 
                            tol=tol,
                            maxIter=maxIter)
  }else if(method=="mean"){
    stop("Not written for method='mean'")
  }
  
  out <- list(
    call=call,
    est=cv_tmle.fit$est, var=cv_tmle.fit$var,
    meanIC=cv_tmle.fit$meanIC, ic=cv_tmle.fit$ic,
    ftimeMod=cv_tmle.fit$ftimeMod, ctimeMod=cv_tmle.fit$ctimeMod,
    trtMod=cv_tmle.fit$trtMod, t0=t0
  )
  class(out) <- "survtmle"
  out
}
