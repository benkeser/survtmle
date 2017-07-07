library(survtmle)
library(survival)
library(cmprsk)
context("Testing checkInputs function")

test_that("checkInputs errors with bad treatment inputs", {
  ftime <- rep(10,3)
  ftype <- rep(10,3)
  trt <- 1:3
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = NULL, glm.trt = "1",
                     glm.ftime = "1",
                     glm.ctime = "1"))
  trt <- cbind(1,2,2)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                   adjustVars = NULL, glm.trt = "1",
                   glm.ftime = "1",
                   glm.ctime = "1"))
})

test_that("NAs throw errors", {
  ftime <- rep(10,3)
  ftime[3] <- NA
  ftype <- rep(10,3)
  trt <- c(0,0,1)
  expect_error(survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = NULL, glm.trt = "1",
                     glm.ftime = "1",
                     glm.ctime = "1"))
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
  n <- 100
  trt <- rbinom(n,1,0.5)
  adjustVars <- data.frame(W1=rep(1,n),W2=rep(0,n))

  ftime <- round(runif(n,0,4)) + 1
  ftype <- round(runif(n,0,1)) +  round(runif(n,0,1))
  
  expect_warning(fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt,
                     adjustVars = adjustVars))
  expect_true(is.numeric(fit1$est))
})


