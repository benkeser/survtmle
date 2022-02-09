# TO DO: document
getSamplingInfluenceCurve <- function(
  dat, ofInterestJ, uniqtrt, mediatorSampProb, R
){
	eif_mat <- matrix(NA, nrow = dim(dat)[1],
	                  ncol = length(ofInterestJ))
	ct <- 0
	for(j in ofInterestJ){
		ct <- ct + 1
		ED <- dat[[paste0("EDfull.j", j, ".z", uniqtrt)]]
		eif_mat[,ct] <- ED / mediatorSampProb * (R - mediatorSampProb)
	}
	eif_dat <- data.frame(eif_mat)
	colnames(eif_dat) <- paste0("Dsamp.j", ofInterestJ, ".z", uniqtrt)
	return(eif_dat)
}