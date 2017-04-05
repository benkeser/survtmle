#' cv_hazard.tmle
#' 
#' This function estimates the marginal cumulative incidence for failures of
#' specified types using targeted minimum loss-based estimation based on the 
#' initial estimates of the cause-specific hazard functions for failures of each type. 
#' The function is called by \code{survtmle} whenever \code{method="hazard"} is specified. 
#' However, power users could, in theory, make calls directly to this function. 
#' 
#' @param ftime A numeric vector of failure times. Missing values or zero values
#' are not supported. Right censored observations should have a non-missing value for \code{ftime} while 
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
#' @param bounds A list of bounds... XXX NEED MORE DESCRIPTION HERE XXX
#' @param verbose A boolean indicating whether the function should print messages to indicate progress.
#' @param tol The stopping criteria for the fluctuation step. The algorithm will continue performing
#' targeting updates to the initial estimators until the empirical mean of the efficient influence function
#' is smaller than \code{tol}. The default (\code{1/length(ftime)}) is a sensible value. Larger values can
#' be used in situations where convergence of the algorithm is an issue; however, this may result in 
#' finite-sample bias. 
#' @param maxIter A maximum number of iterations for the algorithm when \code{method="hazard"}. The 
#' algorithm will iterate until either the empirical mean of the efficient influence function
#' is smaller than \code{tol} or until \code{maxIter} iterations have been completed. 
#' @param ... Other arguments. Not currently used. 
#' 
#' 
#' @return est A numeric vector of point estimates -- one for each combination of \code{ftypeOfInterest}
#' and \code{trtOfInterest}.
#' @return var A covariance matrix of all the point estimates
#' @return meanIC The empirical mean of the efficient influence function at the estimated, fluctuated
#' nuisance parameters. If all goes well, each value should be less than \code{tol}. 
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
#'   
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
#' #' # Fit 1
#' # fit mean.tmle object with glm estimators for treatment, censoring, and failure
#' fit1 <- mean.tmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="hazard", returnModels = TRUE)
#' 

cv_hazard.tmle <- function(
  ftime, 
  ftype,
  trt,
  t0,
  incidence=TRUE, 
  adjustVars=NULL,
  survtmleFits,
  validFolds,
  returnIC=TRUE,
  returnModels=FALSE,
  ftypeOfInterest="all",
  trtOfInterest="all",
  method="hazard",
  bounds=NULL, 
  verbose,
  tol=1/(length(ftime)),
  maxIter=100,
  ...
){
  
  # assemble data frame of necessary variables
  n <- length(ftime)
  id <- 1:n
  dat <- data.frame(id=id, ftime=ftime,ftype=ftype,trt=trt)
  if(!is.null(adjustVars)) dat <- cbind(dat, adjustVars)
  
  # calculate number of failure types
  if(all(ftypeOfInterest=="all")){
    nJ <- length(unique(ftype[ftype!=0]))
    allJ <- ofInterestJ <- sort(unique(ftype[ftype!=0]))
  }else{
    nJ <- length(ftypeOfInterest)
    allJ <- sort(unique(ftype[ftype!=0]))
    ofInterestJ <- sort(ftypeOfInterest)
  }
  
  # calculate number of groups
  if(all(trtOfInterest=="all")){
    ntrt <- length(unique(trt))
    uniqtrt <- sort(unique(trt))
  }else{
    ntrt <- length(trtOfInterest)
    uniqtrt <- sort(trtOfInterest)
  }

  # number of folds
  nfolds <- length(validFolds)
  #------------------------------
  # Work flow for a given split
  # 1. Get treatment probabilities in validation sample from fitted survtmle object
  # 2. Make long version of validation data set 
  # 3. Get censoring probabilities in validation sample from fitted survtmle object
  # 4. Get estimated hazards in validation sample from fitted survtmle object
  # 5. Fluctuate using validation data until convergence
  # 6. Save estimates and influence functions
  # 7. Work on combining on back end 
  #
  #
  cv_targeting_out <- lapply(split(1:nfolds,1:nfolds), function(ind){
    #browser()
    #--------------------------------------------------------------
    # get predicted treatment probabilities in validation sample
    #--------------------------------------------------------------
    # first figure out if glm or superlearner was used
    glm.trt <- NULL; SL.trt <- NULL
    is.glm.trtMod <- "glm" %in% class(survtmleFits[[ind]]$trtMod)
    if(is.glm.trtMod){
      glm.trt <- survtmleFits[[ind]]$trtMod
    }else{
      SL.trt <- survtmleFits[[ind]]$trtMod
    }

    trtOut <- estimateTreatment(dat=dat[validFolds[[ind]],], ntrt=ntrt, uniqtrt=uniqtrt, 
                                adjustVars=adjustVars[validFolds[[ind]],], 
                                SL.trt = SL.trt, glm.trt = glm.trt, returnModels = returnModels)
    dat.valid <- trtOut$dat
    trtMod <- trtOut$trtMod

    #--------------------------------------------------------------
    # make long version of validation data
    #--------------------------------------------------------------
    dataList <- makeDataList(dat=dat.valid, J=allJ, ntrt=ntrt, uniqtrt=uniqtrt, 
                             t0=t0, bounds=bounds)
    
    #--------------------------------------------------------------
    # get predicted censoring in validation
    #--------------------------------------------------------------
    # first figure out if glm or superlearner was used
    glm.ctime <- NULL; SL.ctime <- NULL
    is.glm.ctimeMod <- "glm" %in% class(survtmleFits[[ind]]$ctimeMod)
    if(is.glm.ctimeMod){
      glm.ctime <- survtmleFits[[ind]]$ctimeMod
    }else{
      SL.ctime <- survtmleFits[[ind]]$ctimeMod
    }

    censOut <- estimateCensoring(dataList=dataList, ntrt=ntrt, uniqtrt=uniqtrt, t0=t0,verbose=verbose,
                                  adjustVars=adjustVars[1,], ## NOTE: The estimateCensoring function is kind of stupid
                                  # right now and only looks at the names of this. Should change that to be the option instead 
                                  # so that the whole data frame need not be passed into this function. 
                                  SL.ctime=SL.ctime, glm.ctime=glm.ctime,
                                  returnModels=returnModels)
    dataList <- censOut$dataList
    ctimeMod <- censOut$ctimeMod
    
    #---------------------------------------------------------------
    # get predicted cause specific hazards
    #---------------------------------------------------------------    
    # first figure out if glm or superlearner was used
    glm.ftime <- NULL; SL.ftime <- NULL
    is.glm.ftimeMod <- "glm" %in% class(survtmleFits[[ind]]$ftimeMod[[1]])
    if(is.glm.ftimeMod){
      glm.ftime <- survtmleFits[[ind]]$ftimeMod
    }else{
      SL.ftime <- survtmleFits[[ind]]$ftimeMod
    }

    estOut <- estimateHazards(dataList=dataList, J=allJ, 
                                verbose=verbose,bounds=bounds,
                                adjustVars=adjustVars[1,], # same comment as with estimateCensoring applies here!
                                SL.ftime=SL.ftime, glm.ftime=glm.ftime,
                                returnModels=returnModels)
    dataList <- estOut$dataList
    ftimeMod <- estOut$ftimeMod

    # check for convergence
    suppressWarnings(
      if(all(dataList[[1]]=="convergence failure")){
        return("estimation convergence failure")
      })
      
      return(list(dat.valid, dataList))
    })
    # browser()
    # now let's put things together
    c.dat <- Reduce(rbind, lapply(cv_targeting_out, function(x){x[[1]]}))
    c.dataList <- apply(matrix(1:3), 1, function(x){
      tmp <- lapply(cv_targeting_out, function(y){ y[[2]][[x]] })
      Reduce(rbind, tmp)
    })


    # lapply(cv_targeting_out, function(x){ 
    #   x[[2]] 
    # }),    #FUN = function(...){ Reduce(rbind,list(...))}, 
    #                         lapply(cv_targeting_out, function(x){ x[[2]] }), 
    #                         SIMPLIFY = FALSE)

    # now we're off and targeting in validation sample...

    # calculate cum inc and clever covariates needed for fluctuations
    c.dataList <- updateVariables(dataList=c.dataList, allJ=allJ, ofInterestJ=ofInterestJ, 
                                nJ=nJ, uniqtrt=uniqtrt, ntrt=ntrt, t0=t0, verbose=verbose, cvSieve = TRUE)
    
    # calculate influence function for sample split average parameter
    # this will now have the proper "clever covariates" as a result of the updated
    # updateVariables function.   
    dat.valid <- getHazardInfluenceCurve(dataList=c.dataList, dat=c.dat, ofInterestJ=ofInterestJ, 
                                   allJ=allJ, nJ=nJ, uniqtrt=uniqtrt, ntrt=ntrt, verbose=verbose, t0=t0)
    infCurves <- dat.valid[,grep("D.j",names(dat.valid))]
    meanIC <- colMeans(infCurves)
    
    attr(c.dataList, "fluc") <- rep(Inf,ntrt*nJ^2); ct <- 0
    while(any(abs(meanIC) > tol) & ct <= maxIter){
      ct <- ct + 1
      c.dataList <- fluctuateHazards(dataList=c.dataList, ofInterestJ=ofInterestJ, tol=tol,
                                   allJ=allJ, nJ=nJ, uniqtrt=uniqtrt, ntrt=ntrt, verbose=verbose, t0=t0,
                                   cvSieve = TRUE)
      suppressWarnings(
        if(all(c.dataList[[1]]=="convergence failure")){
          return("fluctuation convergence failure")
        } )
      
      # calculate influence function
      dat.valid <- getHazardInfluenceCurve(dataList=c.dataList, dat=dat.valid, ofInterestJ=ofInterestJ, 
                                     allJ=allJ, nJ=nJ, uniqtrt=uniqtrt, ntrt=ntrt, verbose=verbose, t0=t0)
      infCurves <- dat.valid[,grep("D.j",names(dat.valid))]
      meanIC <- colMeans(infCurves)
      
      if(verbose){
        # print(attr(dataList,"fluc"))
        cat("TMLE Iteration ", ct, "  : ", round(meanIC, 4), "\n")
      }
    }
    if(ct == maxIter+1){
      warning("TMLE fluctuations did not converge. Check that meanIC is adequately small and proceed with caution.")
    }
    
    # calculate point estimate
    est <- rowNames <- NULL
    for(j in ofInterestJ){
      for(z in uniqtrt){
        eval(parse(text=paste("est <- rbind(est, dat.valid$margF",j,".z",z,".t0[1])",sep="")))
        rowNames <- c(rowNames, paste(c(z,j),collapse=" "))
      }
    }
    row.names(est) <- rowNames
    
    # calculate standard error
    var <- t(as.matrix(infCurves))%*%as.matrix(infCurves)/n^2
    row.names(var) <- colnames(var) <- rowNames
    
    out <- list(est=est,var=var, meanIC=meanIC, ic=infCurves, trtMod=NULL, ftimeMod=NULL, ctimeMod=NULL)
    # out

  return(out)
}
