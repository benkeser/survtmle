#' Determine weights for MSM
#' 
#' TO DO: Add more documentation here. 
#' 
#' @param dat
#' @param dataList
#' @param ntrt
#' @param uniqtrt
#' @param adjustVars
#' @param msm.formula
#' @param msm.weights
#' @param returnModels
#' @param ...
#' 
#' @return A list
#' 

estimateMSMWeights <- function(dat, dataList = NULL, ntrt, uniqtrt, adjustVars,
                               msm.formula, msm.weights,
                               returnModels, t0 = NULL){
	long <- !is.null(dataList)
	n <- n.ii <- n.iv <- length(dat[,1])
	if(long){
		n.ii <- length(dataList[[1]][dataList[[1]]$t <= t0, 1])
		n.iv <- length(dataList[[2]][,1])
	}
	n <- length(dat[,1])
	if(msm.weights == "marginal"){
		msmWeightList <- vector(mode = "list", length = ntrt + 1 + long)
		msmWeightList[[1]] <- rep(NA, n)
		ct <- 0
		for(z in uniqtrt){
			ct <- ct + 1
			margz <- mean(dat$trt == z)
			if(!long){
				trt_ind <- dat$trt == z
				msmWeightList[[1]][trt_ind] <- margz
				msmWeightList[[ct + 1]] <- rep(margz, n)	
			}else{
				trt_ind <- dat$trt == z
				msmWeightList[[1]][trt_ind] <- margz
				trt_ind <- dataList[[1]]$trt[dataList[[1]]$t <= t0] == z
				msmWeightList[[2]][trt_ind] <- margz
				msmWeightList[[ct + 2]] <- rep(margz, n.iv)
			}
		}
	}else if(msm.weights == "equal"){
		msmWeightList <- vector(mode = "list", length = ntrt + ifelse(long, 1,2))
		ct <- 0
		for(i in 1:length(msmWeightList)){
			this_n <- n
			if(i == 2) this_n <- n.ii
			if(i > 2) this_n <- n.iv
			msmWeightList[[i]] <- rep(1, this_n)
		}
	}else{
		stop("Only msm.weights = 'marginal' or 'equal' supported at this time.")
	}

	return(msmWeightList)
}