# some thoughts on MSM implementation

# will need to make sure that making data lists handles trtOfInterest
# appropriately. 

# restrict only to one ftypeOfInterest for now...

# put h(a,v) into wideDataList

# need input of h(a,v) 
# default could be to set = 1 
# an alternative that would seemingly always avoid positivity
# would be to fit an additional super learner of g_z onto V
# to generate a compatible estimate of P(Z = z | V)

# how to specify which ftype you want an msm for?

#
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
fit <- survtmle(ftime = dat$ftime, ftype = dat$ftype, adjustVars = dat$adjustVars,
	                t0 = 6,
	                trt = dat$trt, msm.formula = "trt*w2", ftypeOfInterest = 1,
	                glm.trt = "1",
	                glm.ctime = "t",
	                glm.ftime = "w1 + trt*w2",
	                msm.family = "binomial", method = "hazard")

		


getTrueCumInc <- function(w1, w2, trt, t0 = 6){
	haz1 <- plogis(-1 + 0.5*w1 - w2*trt)/2
	haz2 <- plogis(-1 + 0.5*w1 - w2*trt)/2
	ci <- haz1
	for(t in 2:t0){
		ci <- ci + haz1 * (1 - haz1 - haz2)^(t-1)
	}
	return(ci)
}

getTruth <- function(n = 1e5, t0 = 6, msm.formula = "w2*trt"){
	adjustVars <- data.frame(w1 = rnorm(n), w2 = rbinom(n,1,0.5))
	ftype <- rbinom(n, 1, 0.5) + 1

	ci_by_trt <- sapply(0:3, getTrueCumInc, w1 = adjustVars$w1, w2 = adjustVars$w2, t0 = 6,
	                    simplify = FALSE)
	weight_by_trt <- sapply(0:3, function(x){
		rep(dbinom(x, 3, 0.25), n)
	}, simplify = FALSE)

	stackedOutcome <- Reduce("c", ci_by_trt)
	stackedWeights <- Reduce("c", weight_by_trt)
	stackedTrt <- sort(rep(0:3, n))
	stackedAdjustVars <- adjustVars[rep(seq_len(n), 4),]
	stackedData <- data.frame(ci = stackedOutcome, 
	                          weights = stackedWeights,
	                          trt = stackedTrt, stackedAdjustVars)
	fit <- glm(as.formula(paste0("ci ~", msm.formula)), family = binomial(),
	           weights = stackedData$weights, data = stackedData)

	return(fit$coefficients)
}

truth <- getTruth()

do.one <- function(n){
	ct <<- ct + 1
	if(ct%%50 == 0){
		cat("On iteration:", ct," \n")
	}
	dat <- makeData(n = n)
	# debug(mean_tmle)
	# undebug(survtmle:::estimateTreatment)
	# debug(survtmle:::fluctuateIteratedMean)
	fit <- tryCatch({survtmle(ftime = dat$ftime, ftype = dat$ftype, adjustVars = dat$adjustVars,
	                t0 = 6,
	                trt = dat$trt, msm.formula = "trt*w2", ftypeOfInterest = 1,
	                glm.trt = "1",
	                glm.ctime = "t",
	                glm.ftime = "w1 + trt*w2",
	                msm.family = "binomial", method = "mean")
		}, error = function(e){ rep(NA, 8) })
	if(class(fit) == "survtmle"){
		est <- fit$est
		se <- sqrt(diag(fit$var))
		cil <- est - 1.96 * se
		ciu <- est + 1.96 * se
		coverage <- truth < ciu & truth > cil
		bias <- est - truth
		return(c(bias,coverage))
	}else{
		return(fit)
	}	
}

set.seed(1234)
ct <- 0
rslt <- replicate(500, do.one(1000))
