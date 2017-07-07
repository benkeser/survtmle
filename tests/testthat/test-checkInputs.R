library(survtmle)
library(survival)
library(cmprsk)
context("Testing checkInputs function")

test_that("checkInputs errors with bad inputs", {
  ftime <- rep(10,3)
  ftype <- rep(10,3)
  trt <- 1:3
  # too many trt levels
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = NULL, glm.trt = "1",
                     glm.ftime = "1",
                     glm.ctime = "1"))
  # wrong trt format
  trt <- cbind(1,2,2)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                   adjustVars = NULL, glm.trt = "1",
                   glm.ftime = "1",
                   glm.ctime = "1"))
  # no events observed prior to t0
  trt <- c(1,1,0)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                 adjustVars = NULL, glm.trt = "1",
                 glm.ftime = "1",
                 glm.ctime = "1",t0=3)) 
  # adjustVars has column 't'
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
               adjustVars = data.frame(t=1:3), glm.trt = "1",
               glm.ftime = "1",
               glm.ctime = "1"))   

})

test_that("NAs/NULLs throw errors", {
  ftime <- rep(10,3)
  ftime[3] <- NA
  ftype <- rep(10,3)
  trt <- c(0,0,1)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = NULL, glm.trt = "1",
                     glm.ftime = "1",
                     glm.ctime = "1"))
  # add a NULL
  ftime <- NULL
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                   adjustVars = NULL, glm.trt = "1",
                   glm.ftime = "1",
                   glm.ctime = "1"))
  # add missing in adjustVars
  ftime <- rep(10,3)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = data.frame(grbg = c(NA,1,1)), glm.trt = "1",
                     glm.ftime = "1",
                     glm.ctime = "1"))
  # no observed events in trt arm
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                 adjustVars = NULL, glm.trt = "1",
                 glm.ftime = "1",
                 glm.ctime = "1", ftypeOfInterest = 1))

})

test_that("ftime < 0 throw errors", {
  ftime <- rep(10,3)
  ftime[3] <- -10
  ftype <- rep(10,3)
  trt <- c(0,0,1)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = NULL, glm.trt = "1",
                     glm.ftime = "1",
                     glm.ctime = "1"))
})

test_that("adjustVars = NULL works as expected", {
  set.seed(1234)
  n <- 200
  trt <- rbinom(n,1,0.5)

  ftime <- round(runif(n,0,4)) + 1
  ftype <- round(runif(n,0,1)) +  round(runif(n,0,1))
  
  expect_warning(fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = NULL))
  expect_true(is.numeric(fit1$est))
})

test_that("adjustVars constantly valued works as expected", {
  set.seed(1234)
  n <- 200
  trt <- rbinom(n,1,0.5)
  adjustVars <- data.frame(W1=rep(1,n),W2=rep(0,n))

  ftime <- round(runif(n,0,4)) + 1
  ftype <- round(runif(n,0,1)) +  round(runif(n,0,1))
  
  expect_warning(fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = adjustVars))
  expect_true(is.numeric(fit1$est))
  # this one should go ahead fitting ctime and ftime glm's because they don't
  # use W
  expect_warning(fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = adjustVars,glm.trt = "W1",glm.ctime="trt",
                     glm.ftime="trt", t0 =3))
  expect_true(is.numeric(fit2$est))

  # this one should remove SL and fit unadjusted estimates
  expect_warning(fit3 <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = adjustVars,glm.trt = "W1",SL.ctime="SL.glm",
                     SL.ftime="SL.glm", t0 =3))
  expect_true(is.numeric(fit3$est))
})

test_that("when glm and super learner are specified, things behave.", {
  set.seed(1234)
  n <- 100
  trt <- rbinom(n,1,0.5)
  adjustVars <- data.frame(W1=rnorm(n),W2=rnorm(n))

  ftime <- round(runif(n,0,4)) + 1
  ftype <- round(runif(n,0,1)) +  round(runif(n,0,1))
  
  expect_warning(fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = adjustVars,glm.trt = "1",SL.trt = "SL.glm",
                     glm.ctime = "t + W1", SL.ctime = "SL.glm", glm.ftime = "trt + t",
                     SL.ftime = "SL.glm", t0 = 3))
  expect_true(is.numeric(fit1$est))
})

test_that("bad formulas cause errors.", {
  set.seed(1234)
  n <- 100
  trt <- rbinom(n,1,0.5)
  adjustVars <- data.frame(W1=rnorm(n),W2=rnorm(n))

  ftime <- round(runif(n,0,4)) + 1
  ftype <- round(runif(n,0,1)) +  round(runif(n,0,1))
  
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,t0=3,
                     adjustVars = adjustVars,glm.trt = "1245$%",glm.ftime = "1",glm.ctime = "1"))
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,t0=3,
                     adjustVars = adjustVars,glm.ftime = "1245$%",glm.trt = "1",glm.ctime = "1"))
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,t0=3,
                     adjustVars = adjustVars,glm.ctime = "1245$%",glm.trt = "1",glm.ftime = "1"))

})
