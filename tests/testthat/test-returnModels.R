library(survtmle)
library(survival)
context("Checking model types under differing runs with returnModels")

test_that("Objects in fit are of correct classes for differing returnModels", {
	# simulation parameters
	set.seed(1234)
	n <- 200

	# generate data
	trt <- rbinom(n, 1, 0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
	ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
	ftype <- round(runif(n, 0, 1))
	
	# fit iterative mean method with no bounds and return models
	fit1 <- survtmle(ftime = ftime, ftype = ftype,
					 trt = trt, adjustVars = adjustVars,
					 glm.trt = "W1 + W2",
					 glm.ftime = "trt + W1 + W2",
					 glm.ctime = "trt + W1 + W2",
					 method = "mean", t0 = 6,
					 returnModels = TRUE)
	
	# extract objects in model fit
	ftime_mod_1 <- lapply(lapply(fit1$ftimeMod, "[[", 1), class)
	names(ftime_mod_1) <- NULL
	ftime_mod_1 <- unlist(ftime_mod_1)
	ctime_mod_1 <- class(fit1$ctimeMod)
	trt_mod_1 <- class(fit1$trtMod)
	
	# should have roughly same point estimates
	#expect_equivalent(class(ftime_mod_1)

	# fit iterative mean method with no bounds and do NOT return models
	fit2 <- survtmle(ftime = ftime, ftype = ftype,
					 trt = trt, adjustVars = adjustVars,
					 glm.trt = "W1 + W2",
					 glm.ftime = "trt + W1 + W2",
					 glm.ctime = "trt + W1 + W2",
					 method = "mean", t0 = 6,
					 returnModels = FALSE)

	# extract objects in model fit
	ftime_mod_2 <- unlist(fit2$ftimeMod)
	ctime_mod_2 <- fit2$ctimeMod
	trt_mod_2 <- fit2$trtMod

	# should have roughly same point estimates
	#expect_equivalent(ftime_mod_2, ctime_mod_2)
	#expect_equivalent(ftime_mod_2, trt_mod_2)
	#expect_equivalent(trt_mod_2, NULL)
})
