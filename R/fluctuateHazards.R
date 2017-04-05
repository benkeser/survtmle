#' fluctuateHazards
#' 
#' This function performs a fluctuation of an initial estimate of the cause-specific
#' hazard functions using a call to \code{glm} (i.e., a logistic submodel) or a call to
#' \code{optim} (to ensure fluctuations stay within model space).
#' The structure of the function is specific to how it is 
#' called within \code{hazard.tmle}. In particular, \code{dataList} must have a very specific structure for this 
#' function to run properly. The list should consist of \code{data.frame} objects. 
#' The first will have the number of rows for each observation
#' equal to the \code{ftime} corresponding to that observation. The subsequent entries will
#' have \code{t0} rows for each observation and will set \code{trt} column equal to each value
#' of \code{trtOfInterest} in turn. The function will fit a logistic regression with 
#' (a scaled version of) \code{Nj} as outcome, the logit of the current (pseudo-) hazard 
#' estimate as offset and the targeted minimum loss-based estimation "clever covariates". The function
#' then obtains predictions based on this fit on each of the \code{data.frame} objects in 
#' \code{dataList}. 
#' 
#' @param dataList A list of \code{data.frame} objects. 
#' @param allJ Numeric vector indicating the labels of all causes of failure. 
#' @param ofInterestJ Numeric vector indicating \code{ftypeOfInterest} that was passed to 
#' \code{hazard.tmle}. 
#' @param nJ The number of unique failure types. 
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean.tmle}.
#' @param ntrt The number of \code{trt} values of interest. 
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate. 
#' @param verbose A boolean indicating whether the function should print messages to indicate progress.
#' @param ... Other arguments. Not currently used. 
#' 
#' @return The function returns a list that is exactly the same as the input \code{dataList}, 
#' but with updated columns corresponding with estimated cumulative incidence at each time
#' and estimated "clever covariates" at each time. 
#' 
#' @export


fluctuateHazards <- function(
  dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt, t0, verbose, cvSieve = FALSE, ...
){
  eps <- NULL
  for(z in uniqtrt){
    for(j in allJ){
      # clever covariates
      cleverCovariatesNotSelf <- NULL
      if(length(ofInterestJ[ofInterestJ!=j])>0){
        cleverCovariatesNotSelf <- c(cleverCovariatesNotSelf,paste0("H",ofInterestJ[ofInterestJ!=j],".jNotSelf.z",z))
      }
      if(j %in% ofInterestJ){
        cleverCovariatesSelf <- paste0("H",j,".jSelf.z",z)
      }else{
        cleverCovariatesSelf <- NULL
      }
    
    # calculate offset term and outcome
    dataList <- lapply(dataList, function(x,j,allJ){
      eval(parse(text=paste("x$thisScale <- pmin(x$u",j,",1-x$hazNot",j,") - x$l",j,sep="")))
      eval(parse(text=paste("x$thisOffset <- qlogis( pmin((x$Q",j,"Haz - x$l",j,") /x$thisScale ,1-.Machine$double.neg.eps))",sep="")))
      eval(parse(text=paste("x$thisOutcome <- (x$N",j,"- x$l",j,") / x$thisScale",sep="")))
      x
    }, j=j,allJ=allJ)

#    if(length(c(cleverCovariatesNotSelf,cleverCovariatesSelf))>1){
      fluc.mod <- optim(par=rep(0,length(c(cleverCovariatesNotSelf,cleverCovariatesSelf))), 
                        fn=LogLikelihood.offset, 
                        Y=dataList[[1]]$thisOutcome, 
                        H=suppressWarnings(
                          as.matrix(Diagonal(x=dataList[[1]]$thisScale)%*%as.matrix(dataList[[1]][,c(cleverCovariatesNotSelf,cleverCovariatesSelf)]))
                          ),
                        offset=dataList[[1]]$thisOffset,    
                        method="BFGS",gr=grad.offset,
                        control=list(reltol=1e-7,maxit=50000))
#    }else{
    #   fluc.mod <- optim(par=rep(0,length(c(cleverCovariatesNotSelf,cleverCovariatesSelf))), 
    #                     fn=LogLikelihood.offset, 
    #                     Y=dataList[[1]]$thisOutcome, 
    #                     H=as.matrix(Diagonal(x=dataList[[1]]$thisScale)%*%as.matrix(dataList[[1]][,c(cleverCovariatesNotSelf,cleverCovariatesSelf)])),
    #                     offset=dataList[[1]]$thisOffset,    
    #                     method="Brent",lower=-1000,upper=1000,
    #                     control=list(reltol=1e-12,maxit=50000))
    # }
    
    if(fluc.mod$convergence!=0){
      warning("Fluctuation convergence failure. Using with initial estimates. Proceed with caution")
      beta <- rep(0, length(fluc.mod$par))
    }else{
      beta <- fluc.mod$par
    }
    eps <- c(eps, beta)
      
    dataList <- lapply(dataList, function(x,j){
      eval(parse(text=paste("x$Q",j,"PseudoHaz[x$trt==",z,"] <- plogis(x$thisOffset[x$trt==",z,"] + suppressWarnings(as.matrix(Diagonal(x=x$thisScale[x$trt==",z,"])%*%as.matrix(x[x$trt==",z,",c(cleverCovariatesNotSelf,cleverCovariatesSelf)]))%*%as.matrix(beta)))",sep="")))
      eval(parse(text=paste("x$Q",j,"Haz[x$trt==",z,"] <- x$Q",j,"PseudoHaz[x$trt==",z,"] * x$thisScale[x$trt==",z,"] + x$l",j,"[x$trt==",z,"]",sep="")))
      x 
    },j=j)
      
    # update variables based on new haz
    dataList <- updateVariables(dataList=dataList, allJ=allJ, ofInterestJ=ofInterestJ, 
                                nJ=nJ, uniqtrt=uniqtrt, ntrt=ntrt, verbose=verbose, t0=t0,
                                cvSieve = cvSieve)
    }
  }
  attr(dataList,"fluc") <- eps
  dataList
}
