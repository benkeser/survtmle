getMediatorInfluenceCurve <- function(
	dat, uniqtrt, mediatorTrtVal, ofInterestJ
){
	for(j in ofInterestJ){
		ic_label <- paste0("Dmed.j", j, ".z", uniqtrt)
		g_a <- dat[[paste0("g_", mediatorTrtVal)]]
		ind_a <- as.numeric(dat$trt == mediatorTrtVal)
		Fmed_j <- dat[[paste0("F", j, ".z", uniqtrt, ".t0")]]
		F_j <- dat[[paste0("margMedF", j, ".z", uniqtrt, ".t0")]]
		
		dat[,ic_label] <- 
		  dat$sampWt * ind_a / g_a * (Fmed_j - F_j)
	}
	return(dat)
}