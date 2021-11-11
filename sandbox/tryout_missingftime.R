devtools::load_all("~/Dropbox/R/survtmle")
# simulate data
set.seed(1234)
n <- 1000
trt <- rbinom(n, 2, 0.25)
adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
ftime <- round(1 + runif(n, 1, 4) - trt / 2 + adjustVars$W1 + adjustVars$W2)
ftime[1:20] <- NA
ftype <- round(runif(n, 0, 1))
# Fit 1 - fit hazard_tmle object with GLMs for treatment, censoring, failure
# debug(hazard_tmle)
debug(survtmle:::getHazardInfluenceCurve)
fit1 <- survtmle(
  method = "hazard",     
  t0 = 5,
  ftime = ftime, ftype = ftype,
  trt = trt, adjustVars = adjustVars,
  glm.trt = "W1 + W2",
  glm.ftimeMissing = "trt + W1 + W2",
  # SL.trt = c("SL.glm", "SL.mean"),
  glm.ftime = "trt + W1 + W2",
  glm.ctime = "trt + W1 + W2",
  returnModels = TRUE
)

timepoints(fit1, times = 2:5)

