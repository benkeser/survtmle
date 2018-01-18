#' Compute the MSM projection parameter
#' 
#' @param dat
#' @param dataList
#' @param msm.formula 
#' @param msm.family
#' @param msmWeightList
#' @param ...
#' 
getHazardMSMEstimate <- function(dat, dataList, msm.formula, msm.family, 
                                 msmWeightList, adjustVarNames, ofInterestJ, ...){
	min.t <- min(dataList[[1]]$t)
	outcomeName <- paste0("F",ofInterestJ,".z.t0")
	stackedData <- Reduce("rbind", lapply(dataList[2:length(dataList)], function(d){
		d[d$t == min.t, c("trt", adjustVarNames, outcomeName)]
	}))
	stackedWeights <- Reduce("c", mapply(dl = dataList[2:length(dataList)],
	                                     wl = msmWeightList[3:length(msmWeightList)],
	                                     FUN = function(dl, wl){
	                                     	wl[dl$t == min.t]
	                                     }))
	suppressWarnings(
		msm.fit <- glm(as.formula(paste0(outcomeName, " ~ ", msm.formula)),
	    	           data = stackedData, family = msm.family,
	    	           weights = stackedWeights)
	)
	return(msm.fit$coefficients)
}