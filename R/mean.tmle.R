#' mean_tmle
#' 
#' This function estimates the marginal cumulative incidence for failures of
#' specified types using targeted minimum loss-based estimation based on the 
#' G-computation representation of cumulative incidence. The function is called
#' by \code{survtmle} whenever \code{method="mean"} is specified. However, power 
#' users could, in theory, make calls directly to this function. 
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
#' @param bounds A list of bounds... XXX NEED MORE DESCRIPTION HERE XXX
#' @param verbose A boolean indicating whether the function should print messages to indicate progress.
#' @param ... Other arguments. Not currently used. 
#' 
#' 
#' @return call The call to \code{survtmle}.
#' @return est A numeric vector of point estimates -- one for each combination of \code{ftypeOfInterest}
#' and \code{trtOfInterest}.
#' @return var A covariance matrix of all the point estimates
#' @return meanIC The empirical mean of the efficient influence function at the estimated, fluctuated
#' nuisance parameters. If all goes well, each value should be small. This should be confirmed, 
#' particularly if \code{bounds} were specified. 
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
#' # fit mean_tmle object with glm estimators for treatment, censoring, and failure
#' fit1 <- mean_tmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2")
#' 


mean_tmle <- function(
  ftime, 
  ftype,
  trt,
  t0=max(ftime[ftype > 0]),
  incidence=TRUE,
  adjustVars=NULL,
  SL.ftime=NULL,
  SL.ctime=NULL,
  SL.trt=NULL,
  glm.ftime=NULL,
  glm.ctime=NULL,
  glm.trt="1",
  returnIC=TRUE,
  returnModels=FALSE,
  ftypeOfInterest="all",
  trtOfInterest="all",
  bounds=NULL, 
  verbose = FALSE,
  ...
){
  # assemble data frame of necessary variables
  n <- length(ftime)
  id <- 1:n
  dat <- data.frame(id=id, ftime=ftime,ftype=ftype, trt=trt)
  
  if(!is.null(adjustVars)){
    dat <- cbind(dat, adjustVars)
  } 
  
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
  
  # estimate trt probabilities
  trtOut <- estimateTreatment(dat=dat, ntrt=ntrt, uniqtrt=uniqtrt, adjustVars=adjustVars, 
                           SL.trt=SL.trt, glm.trt = glm.trt,
                           returnModels = returnModels)
  dat <- trtOut$dat
  trtMod <- trtOut$trtMod
  
  # make long version of data sets needed for estimation of censoring
  dataList <- makeDataList(dat=dat, J=allJ, ntrt=ntrt, uniqtrt=uniqtrt, t0=t0, bounds=bounds)
  
  # estimate censoring
  censOut <- estimateCensoring(dataList=dataList, ntrt=ntrt, uniqtrt=uniqtrt, t0=t0,verbose=verbose,
                                adjustVars=adjustVars,
                                SL.ctime=SL.ctime,glm.ctime=glm.ctime,
                                returnModels=returnModels)
  dataList <- censOut$dataList
  ctimeMod <- censOut$ctimeMod
  
  wideDataList <- makeWideDataList(dat=dat, dataList=dataList,adjustVars=adjustVars, t0=t0, allJ=allJ, ntrt=ntrt, uniqtrt=uniqtrt)
  
  # estimate/fluctuate iterated means
  timeAndType <- expand.grid(t0:1,ofInterestJ)
  
  # empty list for Qmod if returnModels
  ftimeMod <- vector(mode="list",length=length(ofInterestJ))
  names(ftimeMod) <- paste0("J",ofInterestJ)
  for(j in 1:length(ofInterestJ)){
    ftimeMod[[j]] <- vector(mode="list",length=t0)
    names(ftimeMod[[j]]) <- paste0("t",1:t0)
  }

  for(i in 1:nrow(timeAndType)){
    estOut <- estimateIteratedMean(wideDataList=wideDataList, t=timeAndType[i,1],whichJ=timeAndType[i,2],ntrt=ntrt,uniqtrt=uniqtrt,allJ=allJ,t0=t0,
                                         SL.ftime=SL.ftime, adjustVars=adjustVars,
                                         glm.ftime=glm.ftime,verbose=verbose,
                                         returnModels=returnModels, bounds=bounds)
    wideDataList <- estOut$wideDataList
    eval(parse(text=paste0("ftimeMod$J",timeAndType[i,2],"$t",timeAndType[i,1],"<-estOut$ftimeMod")))
    wideDataList <- fluctuateIteratedMean(wideDataList=wideDataList,t=timeAndType[i,1],whichJ=timeAndType[i,2],ntrt=ntrt,uniqtrt=uniqtrt,allJ=allJ,t0=t0,
                                          SL.ftime=SL.ftime, 
                                          glm.ftime=glm.ftime,
                                          returnModels=returnModels, bounds=bounds)
  }
  
  
  # get point estimate
  est <- rowNames <- NULL
  for(j in ofInterestJ){
    for(z in 1:length(uniqtrt)){
      thisEst <- eval(parse(text=paste("mean(wideDataList[[",z+1,"]]$Q",j,"star.1)",sep="")))
      est <- rbind(est,thisEst)
      rowNames <- c(rowNames, paste(c(uniqtrt[z],j),collapse=" "))
      eval(parse(text=paste("wideDataList[[1]]$Q",j,"star.0.Z",uniqtrt[z]," <- rep(thisEst,n)",sep="")))
      eval(parse(text=paste("wideDataList[[1]]$Q",j,"star.1.Z",uniqtrt[z]," <- wideDataList[[(z+1)]]$Q",j,"star.1",sep="")))
    }
  }
  row.names(est) <- rowNames
  
  # calculate influence function
  for(j in ofInterestJ){
    for(z in 1:length(uniqtrt)){
      for(t in t0:1){
        outcomeName <- ifelse(t==t0, paste("N",j,".",t0,sep=""), paste("Q",j,"star.",t+1,sep=""))
        eval(parse(text=paste(
          "wideDataList[[1]]$D.Z",uniqtrt[z],".",j,"star.",t,
          " <- wideDataList[[1]]$H",uniqtrt[z],".",t,"*(wideDataList[[1]][,outcomeName] - wideDataList[[1]]$Q",j,"star.",t,")",sep="")))
      }
      eval(parse(text=paste("wideDataList[[1]]$D.Z",uniqtrt[z],".",j,"star.0 <- wideDataList[[1]]$Q",j,"star.1.Z",uniqtrt[z]," - wideDataList[[1]]$Q",j,"star.0.Z",uniqtrt[z],sep="")))
      ind <- eval(parse(text=paste("grep('D.Z",uniqtrt[z],".",j,"star', names(wideDataList[[1]]))",sep="")))
      eval(parse(text=paste("wideDataList[[1]]$IC",j,"star.Z",uniqtrt[z],
                            " <- rowSums(cbind(rep(0, nrow(wideDataList[[1]])),wideDataList[[1]][,ind]))",sep="")))
    }
  }
  
  # calculate standard error
  infCurves <- wideDataList[[1]][,grep("IC",names(wideDataList[[1]]))]
  #   browser()
  meanIC <- apply(infCurves, MARGIN=2, FUN=mean)
  var <- t(as.matrix(infCurves))%*%as.matrix(infCurves)/n^2
  row.names(var) <- colnames(var) <- rowNames
  
  return(list(est=est,var=var, meanIC=meanIC, ic=infCurves, trtMod=trtMod,
              ftimeMod=ftimeMod,ctimeMod=ctimeMod)) 
}
