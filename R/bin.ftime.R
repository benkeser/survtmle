#' bin_ftime
#' 
#' This is a helper function that takes a vector of continuously valued failure times and
#' discretizes into \code{nKnots} number of discrete times. By default the bins
#' will be evenly spaced, but users can also input their own cutpoints through the
#' \code{bins} argument. 
#' 
#' @param ftime A vector of failure times
#' @param bins How should the failure times be discretized? Can be set to 
#' \code{"default"} in which case cut points will be evenly over the range of \code{ftime}. 
#' Or can be set to a vector containing the desired cutpoints for the failure time bins.
#' @param nKnots If \code{bins="default"}, how many bins of failure times should there be?
#' The default value is 10 and the option is ignored if \code{bins != "default"}
#' 
#' @return A vector of binned failure times
#' 
#' @examples
#' ftime <- runif(100,0,10)
#' # first do default binning
#' binned.ftime <- bin_ftime(ftime = ftime)
#' table(binned.ftime)
#' # now specify own cutpoints
#' binned.ftime2 <- bin_ftime(ftime = ftime, bins = c(2,5,8))
#' table(binned.ftime)

bin_ftime <- function(ftime, bins="default", nKnots=10){
  n <- length(ftime)
  if(bins=="default"){
    # check for proper input
    if(sum(is.na(ftime))>0) stop("Missing data in ftime not accepted")
    if(!is.numeric(nKnots)) stop("nKnots must be numeric")
    
    cutPoints <- seq(0, max(ftime), length=nKnots)
  }else{
    if(!is.numeric(ftime)) stop("bins must be numeric vector or set to default")
    cutPoints <- bins
  }
  
  # generate new failure time variable
  tTilde <- rep(NA,length(ftime))
  for(i in 1:(length(cutPoints)-1)){
    ind <- which(ftime > cutPoints[i] & ftime <= cutPoints[i+1])
    tTilde[ind] <- i
  }
  
  attr(tTilde, "cutPoints") <- cutPoints
  
  return(tTilde)
}