# TO DO: document
fluctuateCumInc <- function(
  dat, ofInterestJ, uniqtrt,
  mediatorTrtVal
){
  for(j in ofInterestJ){
  	outcome_j <- paste0("F", j, ".z", uniqtrt, ".t0")
  	offset_label_j <- paste0("margMedF", j, ".z", uniqtrt, ".t0")
    dat[,offset_label_j][
      data[,offset_label_j] < .Machine$double.neg.eps
    ] <- .Machine$double.neg.eps

    dat[,offset_label_j][
      data[,offset_label_j] > 1 - .Machine$double.neg.eps
    ] <- 1 - .Machine$double.neg.eps
  	
    g_label <- paste0("g_", uniqtrt)
  	fit_data <- data.frame(
      outcome = dat[, outcome_j],
      logit_F = qlogis(dat[,offset_label_j]),
      covar = as.numeric(dat$trt == mediatorTrtVal) / dat[,g_label],
      in_weight = dat$sampWt
  	)
  	pred_data <- data.frame(
      outcome = dat[, outcome_j],
      logit_F = qlogis(dat[,offset_label_j]),
      covar = 1 / dat[,g_label],
      in_weight = dat$sampWt
  	)
  	suppressWarnings(
	  	fit_j <- glm("outcome ~ -1 + offset(logit_F) + covar", 
	  	             data = fit_data, weights = in_weight,
	  	             family = stats::binomial(), 
	  	             start = 0)
  	)
  	dat[,offset_label_j] <-
  	  predict(fit_j, newdata = pred_data, type = "response") 
  }
  return(dat)
}