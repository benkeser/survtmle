# TO DO: document
estimateFullDataEIF <- function(
	dat, mediator, ofInterestJ, uniqtrt,
	adjustVars, SL.eif, glm.eif, cvControl,
	verbose
){
	j <- ofInterestJ # TODO: generalize for competing risks
	eifMod <- vector(mode = 'list', length = length(ofInterestJ))
	names(eifMod) <- paste0("J", j)
	for(j in ofInterestJ){
		eif_name <- paste0("Dfull.j", j, ".z", uniqtrt)
		estimated_eif_name <- paste0("EDfull.j", j, ".z", uniqtrt)

		fit_idx <- which(complete.cases(mediator))
		ic1_label <- paste0("D.j", j, ".z", uniqtrt)
		ic2_label <- paste0("Dmed.j", j, ".z", uniqtrt)
		# ic_label variables contain observed data influence function
		# which equals full data EIF * sampWt, so we divide by sampWt
		outcome_j <- (dat[fit_idx, ic1_label] + dat[fit_idx, ic2_label]) / dat$sampWt[fit_idx]
		covariates_j <- dat[fit_idx , c(colnames(adjustVars), "trt", "ftime", "ftype")]

		outcome_j_full <- dat[,ic1_label] # doesn't matter
		covariates_j_full <- dat[ , c(colnames(adjustVars), "trt", "ftime", "ftype")]
		# browser()
		if(is.null(SL.eif)){
			fit_data <- data.frame(outcome_j, covariates_j)
			pred_data <- data.frame(outcome_j_full, covariates_j_full)
			colnames(fit_data)[1] <- eif_name
			str_formula <- paste0(eif_name, " ~ ", glm.eif)
			fit_j <- glm(str_formula, data = fit_data, family = gaussian())
			dat[, estimated_eif_name] <- 
				predict(fit_j, newdata = pred_data, type = "response")
		}else{
			fit_j <- SuperLearner(
        Y = outcome_j,
        X = covariates_j,
       	newX = covariates_j_full,
       	SL.library = SL.eif,
       	family = gaussian(),
       	cvControl = cvControl,
       	verbose = verbose
      )
      dat[, estimated_eif_name] <- fit_j$SL.pred
		}
		eifMod[[paste0("J", j)]] <- fit_j
	}
	out <- list(dat = dat, eifMod = eifMod)
	return(out)
}