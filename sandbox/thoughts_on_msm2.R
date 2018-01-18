install.packages("~/Dropbox/R/survtmle", type = "source", repos = NULL)
q("no")

library(survtmle)

makeData <- function(n){
	adjustVars <- data.frame(w1 = rnorm(n), w2 = rbinom(n,1,0.5))
	trt <- rbinom(n, 3, 0.25)
	ftype <- rbinom(n, 1, 0.5) + 1 
	ctime <- 1 + rgeom(n, plogis(-4)) # random censoring
	ftime <- 1 + rgeom(n, plogis(-1 + 0.5*adjustVars$w1 - adjustVars$w2*trt))
	time <- pmin(ctime, ftime)
	ftype[ctime < ftime] <- 0
	return(list(adjustVars = adjustVars, trt = trt, ftime = time, ftype = ftype))
}

dat <- makeData(n = 1000)
debug(hazard_tmle)
# debug(survtmle:::updateVariables)
debug(survtmle:::getHazardInfluenceCurve)
fit <- survtmle(ftime = dat$ftime, ftype = dat$ftype, adjustVars = dat$adjustVars,
	                t0 = 6,
	                trt = dat$trt, msm.formula = "trt*w2", ftypeOfInterest = 1,
	                glm.trt = "1",
	                glm.ctime = "t",
	                glm.ftime = "w1 + trt*w2",
	                msm.family = "binomial", method = "hazard")
