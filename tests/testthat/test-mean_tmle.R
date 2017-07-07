library(survtmle)
library(survival)
context("Testing mean_tmle function")

test_that("mean_tmle with bounds of (0,1) gives same results as unbounded with one failure type", {
	set.seed(1234)
	n <- 200
	trt <- rbinom(n,1,0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
	
	ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
	ftype <- round(runif(n,0,1))
	
	# fit with no bounds
	fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "W1 + W2", 
	glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
	method="mean", t0=6)
	# fit with bounds
	bf <- data.frame(t=1:6,l1 = rep(0,6),u1 = rep(1,6))
	fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "W1 + W2", 
	glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
	method="mean", t0=6, bounds = bf)

	# should have roughly same point estimates
	expect_true(all(abs(fit1$est - fit2$est) < 1e-4))
})


test_that("mean_tmle with glm and super learner with only glm give same answers (one failure type)", {
	set.seed(1234)
	n <- 200
	trt <- rbinom(n,1,0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
	
	ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
	ftype <- round(runif(n,0,1))
	
	# fit with super learner
	fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	SL.trt = "SL.glm", 
	SL.ftime = "SL.glm", SL.ctime = "SL.glm",
	method="mean", t0=6)
	

	# fit with glm
	fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "W1 + W2", 
	glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2 + t", 
	method="mean", t0=6)

	# should have roughly same point estimates
	expect_equal(fit1$est,fit2$est)
})

test_that("mean_tmle with glm and super learner with only glm give same answers (two failure types)", {
	set.seed(1234)
	n <- 200
	trt <- rbinom(n,1,0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
	
	ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
	ftype <- round(runif(n,0,1)) + round(runif(n,0,1))
	
	# fit with super learner
	fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	SL.trt = "SL.glm", 
	SL.ftime = "SL.glm", SL.ctime = "SL.glm",
	method="mean", t0=6)


	# fit with glm
	fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "W1 + W2", 
	glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2 + t", 
	method="mean", t0=6)

	# should have roughly same point estimates
	expect_equal(fit1$est,fit2$est)
})

test_that("mean_tmle with one ftypeOfInterest and one trtOfInterest gives same answer as multiple.", {
	set.seed(1234)
	n <- 200
	trt <- rbinom(n,1,0.5)
	adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
	
	ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
	ftype <- round(runif(n,0,1)) + round(runif(n,0,1))
	
	# fit with super learner
	fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "1", ftypeOfInterest = 1, trtOfInterest = 0,
	glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2",
	method="mean", t0=3)


	# fit with glm
	fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
	glm.trt = "1", 
	glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
	method="mean", t0=3)

	# should have roughly same point estimates
	expect_true(abs(fit1$est[1]-fit2$est[1]) < 1e-4)
})
