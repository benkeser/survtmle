make_data <- function(n){
	W1 <- rbinom(n, 1, 0.5)
	W2 <- rbinom(n, 1, 0.5)
	gA <- plogis(-2 + W1 + W2)
	A <- rbinom(n, 2, gA)

	haz <- plogis(-2 + W1 + W2 + (A == 1) + (A == 2) / 2)
	ftime <- rgeom(n, haz) + 1
	ftype <- rep(1, n)

	gftimeMissing <- plogis(1 + W1 + W2 - (A == 1))
	ftimeMissing <- rbinom(n, 1, gftimeMissing)

	ftime[ftimeMissing == 0] <- NA
	ftype[ftimeMissing == 0] <- 0

	out <- list(adjustVars = data.frame(W1, W2),
	            ftime = ftime,
	            ftype = ftype,
	            trt = A)
	return(out)
}

get_truth <- function(n = 1e6, t0 = 5){
	W1 <- rbinom(n, 1, 0.5)
	W2 <- rbinom(n, 1, 0.5)
	truth <- NULL
	for(a in 0:2){
		haz <- plogis(-2 + W1 + W2 + (a == 1) + (a == 2) / 2)
		ftime <- rgeom(n, haz) + 1
		truth <- c(truth, mean(ftime <= 5))
	}
	return(truth)
}

do_one_sim <- function(n, t0 = 5){
	data <- make_data(n = n)
	fit1 <- survtmle(
	  method = "hazard",     
	  t0 = 5,
	  ftime = data$ftime, ftype = data$ftype,
	  trt = data$trt, adjustVars = data$adjustVars,
	  glm.trt = "W1 + W2",
	  glm.ftimeMissing = "factor(trt) + W1 + W2",
	  # SL.trt = c("SL.glm", "SL.mean"),
	  glm.ftime = "factor(trt) + W1 + W2",
	  glm.ctime = "1",
	  returnModels = FALSE
	)
	est <- fit1$est
	ci <- confint(fit1)
	out <- c(est[1], ci[1,],
	         est[2], ci[2,],
	         est[3], ci[3,])
	return(out)
}


devtools::load_all("~/Dropbox/R/survtmle")
library(future.apply)
plan(multicore, workers = 25)

rslt_1000 <- future_sapply(1:500, do_one_sim, n = 1000)
rslt_2000 <- future_sapply(1:500, do_one_sim, n = 2000)
rslt_4000 <- future_sapply(1:500, do_one_sim, n = 4000)

truth <- get_truth()
bias_1000 <- rowMeans(rslt_1000[c(1, 4, 7),]) - truth
bias_2000 <- rowMeans(rslt_2000[c(1, 4, 7),]) - truth
bias_4000 <- rowMeans(rslt_4000[c(1, 4, 7),]) - truth

cover_1000 <- c(
  mean(rslt_1000[2,] < truth[1] & rslt_1000[3,] > truth[1]),
  mean(rslt_1000[5,] < truth[2] & rslt_1000[6,] > truth[2]),
  mean(rslt_1000[8,] < truth[3] & rslt_1000[9,] > truth[3])
)

cover_2000 <- c(
  mean(rslt_2000[2,] < truth[1] & rslt_2000[3,] > truth[1]),
  mean(rslt_2000[5,] < truth[2] & rslt_2000[6,] > truth[2]),
  mean(rslt_2000[8,] < truth[3] & rslt_2000[9,] > truth[3])
)

cover_4000 <- c(
  mean(rslt_4000[2,] < truth[1] & rslt_4000[3,] > truth[1]),
  mean(rslt_4000[5,] < truth[2] & rslt_4000[6,] > truth[2]),
  mean(rslt_4000[8,] < truth[3] & rslt_4000[9,] > truth[3])
)







