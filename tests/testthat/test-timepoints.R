library(survtmle)
library(survival)
library(cmprsk)
context("Testing timepoints function")


test_that("hazard_tmle and mean_tmle timepoints equal kaplan-meier with no covariates", {
	set.seed(1234)
	n <- 500
	trt <- rbinom(n,1,0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
	
	ftime <- round(runif(n,0,4)) + 1
	ftype <- round(runif(n,0,1))
	
	# hazard fit
	suppressWarnings(
	fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "1", 
	glm.ftime = paste0("trt + ",paste0("I(t==",1:max(ftime),")",collapse="+"),"+",paste0("I(trt*t==",1:max(ftime),")",collapse="+")),
	glm.ctime = paste0("trt + ",paste0("I(t==",1:max(ftime),")",collapse="+"),"+",paste0("I(trt*t==",1:max(ftime),")",collapse="+")), 
	method="hazard", t0=5, returnModels = TRUE)
	)
	suppressWarnings(
		tp.fit1 <- survtmle::timepoints(fit1, times = 1:5)
	)
	est.fit1 <- Reduce(cbind,lapply(tp.fit1,"[[","est"))
	# mean fit
	fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "1", 
	glm.ftime = "trt", 
	glm.ctime = "trt", 
	method="mean", t0=5)
	suppressWarnings(
		tp.fit2 <- timepoints(fit2, times = 1:5)
	)
	est.fit2 <- Reduce(cbind,lapply(tp.fit2,"[[","est"))
	# compare to kaplan meier
	km <- matrix(1-summary(survfit(Surv(ftime,ftype)~trt))$surv,nrow=2,byrow=TRUE)

	expect_equal(as.numeric(km),as.numeric(est.fit1))
	expect_equal(as.numeric(km),as.numeric(est.fit2))
})

test_that("hazard_tmle and mean_tmle equal aalen-johansen with no covariates", {
	set.seed(1234)
	n <- 500
	trt <- rbinom(n,1,0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
	
	ftime <- round(runif(n,0,4)) + 1
	ftype <- round(runif(n,0,1)) + round(runif(n,0,1))
	
	# hazard fit
	suppressWarnings(
	fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "1", 
	glm.ftime = paste0("trt + ",paste0("I(t==",1:max(ftime),")",collapse="+"),"+",paste0("I(trt*t==",1:max(ftime),")",collapse="+")),
	glm.ctime = paste0("trt + ",paste0("I(t==",1:max(ftime),")",collapse="+"),"+",paste0("I(trt*t==",1:max(ftime),")",collapse="+")), 
	method="hazard", t0=5, returnModels = TRUE)
	)
	suppressWarnings(
		tp.fit1 <- survtmle::timepoints(fit1, times = 1:5)
	)
	est.fit1 <- Reduce(cbind,lapply(tp.fit1,"[[","est"))
	# mean fit
	fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "1", 
	glm.ftime = "trt", 
	glm.ctime = "trt", 
	method="mean", t0=5)
	suppressWarnings(
		tp.fit2 <- timepoints(fit2, times = 1:5)
	)
	est.fit2 <- Reduce(cbind,lapply(tp.fit2,"[[","est"))

	# compare to kaplan meier
	aj <- cuminc(ftime = ftime, fstatus = ftype, group = trt)
	fit.aj <- cmprsk::timepoints(aj,5)

	expect_equal(as.numeric(fit.aj$est),as.numeric(fit1$est))
	expect_equal(as.numeric(fit.aj$est),as.numeric(fit2$est))
})
