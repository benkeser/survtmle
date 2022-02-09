# TO DO: document
estimateMedReg <- function(
  dat, 
  SL.mediator,
  glm.mediator,
  adjustVars,
  cvControl,
  returnModels,
  bounds,
  ofInterestJ,
  mediatorTrtVal,
  uniqtrt
){
	if(!is.null(bounds)){
		stop("bounds not yet supported for mediation analysis")
	}
	mediatorMod <- vector(mode = 'list', length = length(ofInterestJ))
	names(mediatorMod) <- paste0("J", ofInterestJ)
	for(j in ofInterestJ){
		col_label_outcome <- paste0("F", j, ".z", uniqtrt, ".t0")
		subset_idx <- which(dat$trt == mediatorTrtVal)
		outcome_j <- dat[subset_idx, col_label_outcome]
		covariates_j <- dat[subset_idx, colnames(adjustVars), drop = FALSE]
		samp_weights_j <- dat$sampWt[subset_idx]

		outcome_j_full <- dat[, col_label_outcome]
		covariates_j_full <- dat[, colnames(adjustVars), drop = FALSE]
		samp_weights_j_full <- dat$sampWt

		new_col_label <- paste0("margMedF", j, ".z", uniqtrt, ".t0")
		
		if(is.null(SL.mediator)){
			fit_data_j <- data.frame(outcome_j, covariates_j, samp_weights_j)
			pred_data_j <- data.frame(outcome_j_full, covariates_j_full, samp_weights_j_full)
			colnames(fit_data_j)[1] <- col_label_outcome
			str_formula_j <- paste0(
			   col_label_outcome, " ~ ", glm.mediator
			)
			suppressWarnings(
				fit_j <- glm(str_formula_j, data = fit_data_j, 
				             family = stats::binomial(),
				             weights = samp_weights_j)
			)
			
			dat[,new_col_label] <- predict(
			   fit_j, newdata = pred_data_j, type = "response"
			)
		}else{
			fit_j <- SuperLearner(
			  Y = outcome_j, X = covariates_j, 
			  newX = covariates_j_full,
			  SL.library = SL.mediator,
			  family = binomial(), # should this be gaussian()?
			  cvControl = cvControl,
			  obsWeights = samp_weights_j,
			  verbose = verbose
      )
			dat[,new_col_label] <- fit_j$SL.pred
		}
		mediatorMod[[paste0("J", j)]] <- fit_j
	}
	out <- list(
    dat = dat,
    mediatorMod = mediatorMod
  )
  return(out)
}