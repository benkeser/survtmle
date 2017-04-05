#' updateVariables
#' 
#' A helper function that maps hazard estimates into estimates of cumulative 
#' incidence and updates the "clever covariates" used by the targeted minimum 
#' loss-based estimation fluctuation step. 
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


updateVariables <- function(
  dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt, t0, verbose, ...
){
  
  dataList <- lapply(dataList, function(x, allJ){
    # total hazard
    Q_dot <- rowSums(cbind(rep(0,nrow(x)), x[,paste0("Q",allJ,"Haz")]))
    # survival at t
    x$S.t <- unlist(by((1-Q_dot), x$id, FUN=cumprod))
    x$S.t[x$S.t==0] <- .Machine$double.neg.eps
    # survival at t-1
    S.tminus1 <- c(1, x$S.t[1:(length(x$S.t)-1)])
    S.tminus1[x$t==1] <- 1
    
    for(j in ofInterestJ){
      # calculate CIF at time t
      eval(parse(text=paste("x$F",j,".t <- unlist(by(x[,","paste0('Q',j,'Haz')","]*S.tminus1, x$id, FUN=cumsum))",sep="")))
    }
    x
  }, allJ=allJ)
  
  # calculate CIF at time t0
  for(j in ofInterestJ){
    Fj.t0.allZ <- NULL
    for(i in 1:ntrt){
      Fj.t0.allZ <- cbind(Fj.t0.allZ, eval(parse(text=paste("dataList[[i+1]]$F",j,".t[dataList[[i+1]]$t==t0]",sep=""))))
    }
    
    dataList <- lapply(dataList, function(x,j,uniqtrt,Fj.t0.allZ){
      for(i in 1:length(uniqtrt)){
        ind <- tapply(X=x$id,INDEX=x$id,FUN=NULL)
        eval(parse(text=paste("x$F",j,".z",uniqtrt[i],".t0 <- Fj.t0.allZ[ind,",i,"]",sep="")))
      }
      x
    },j=j,uniqtrt=uniqtrt,Fj.t0.allZ=Fj.t0.allZ)
  }
  
  dataList <- lapply(dataList, function(x,allJ){
    for(j in allJ){
      if(length(allJ)>1){
        eval(parse(text=paste("x$hazNot",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',allJ[allJ != j],'Haz')]))",sep="")))
        eval(parse(text=paste("x$hazNot",j,"[x$hazNot",j,"==1] <- 1-.Machine$double.neg.eps",sep="")))
      }else{
        eval(parse(text=paste("x$hazNot",j," <- 0",sep="")))
      }
    }
    x
  },allJ=allJ)
  
  
  # set up clever covariates needed for fluctuation
  dataList <- lapply(dataList, function(x,ofInterestJ,uniqtrt){
    for(z in uniqtrt){
      for(j in ofInterestJ){
        eval(parse(text=paste("x$H",j,".jSelf.z",z,"<- (x$ftime>=x$t & x$trt==",z,")/(x$g_",z,"*x$G_dC) * (1-x$hazNot",j,") * ((x$t<t0)*(1-(x$F",j,".z",z,".t0 - x$F",j,".t)/c(x$S.t)) + (x$t==t0))",sep="")))
        eval(parse(text=paste("x$H",j,".jNotSelf.z",z,"<- -(x$ftime>=x$t & x$trt==",z,")/(x$g_",z,"*x$G_dC) *(1-x$hazNot",j,") * ((x$t<t0)*(x$F",j,".z",z,".t0 - x$F",j,".t)/c(x$S.t))",sep="")))
      }
    }
    x
  },ofInterestJ=ofInterestJ,uniqtrt=uniqtrt)
  
  dataList
}
