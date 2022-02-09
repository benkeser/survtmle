fluctuateSampWt <- function(
  dat, ofInterestJ, uniqtrt, R,
  mediatorSampProb
){
	covariate_names <- paste0("EDfull.j", ofInterestJ, ".z", uniqtrt)
	mediatorSampProb[mediatorSampProb > 1 - .Machine$double.neg.eps] <- 1 - .Machine$double.neg.eps
  mediatorSampProb[mediatorSampProb < .Machine$double.neg.eps] <- .Machine$double.neg.eps

	fit_data <- data.frame(
    R = R, logit_sampProb = qlogis(mediatorSampProb),
    sampWt_in = 1 / mediatorSampProb,
    dat[ , covariate_names, drop = FALSE]
  )
  start_vals <- rep(0, length(covariate_names))
  str_form <- paste0(
    "R ~ -1 + offset(logit_sampProb) + ",
    paste0(covariate_names, collapse = "+")
  )
  # browser()
  fit <- glm(str_form, family = stats::binomial(), 
             data = fit_data, start = start_vals,
             weights = sampWt_in)
  targeted_mediatorSampProb <- fit$fitted.values
  return(targeted_mediatorSampProb)
}