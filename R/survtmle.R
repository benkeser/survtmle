#' survtmle
#' 
#' This function estimates the marginal cumulative incidence for failures of
#' specified types using targeted minimum loss-based estimation. More description to come...
#' 


#' @export
#' 
#' @examples
#' 
#' ## Single failure type examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' trt <- rbinom(n,1,0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
#' 
#' ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n,0,1))
#' 
#' # Fit 1
#' # fit a survtmle object with glm estimators for treatment, censoring, and failure
#' # using the "mean" method
#' fit1 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", returnModels = TRUE)
#' # fit1
#' 
#' # Fit 2
#' # fit an survtmle object with SuperLearner estimators for failure and censoring
#' # and empirical estimators for treatment (default) using the "mean" method
#' fit2 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="mean", returnModels = TRUE)
#' # fit2
#' 
#' # Fit 3
#' # repeat Fit 1 using the "hazard" method 
#' fit3 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="hazard", returnModels = TRUE)
#' # fit3
#' 
#' # Fit 4
#' # repeat Fit 2 using the "hazard" method
#' fit4 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="hazard", returnModels = TRUE)
#' # fit4
#' ## Examples with bounded TMLE
#' # Fit 5
#' # repeat Fit 1, but now specifying bounds on the iterated conditional means
#' fit5 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", 
#' # one type of failure so bounds should have columns 't','l1',and 'u1'
#' bounds = data.frame(t=1:8,l1=rep(0.01,8),u1=rep(0.9,8)))
#' # fit5
#' # Fit 6
#' # repeat Fit 10 using the bounded hazard method
#' 
#' # make a data.frame of bounds in proper format
#' # one type of failure so bounds should have columns 't','l1',and 'u1'
#' # and contain t=1:t0 
#' bf1 <- data.frame(t=1:5,l1=rep(0.01,5),u1=rep(0.9,5))
#' fit6 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", 
#' bounds = bf1)
#' # fit6
#' 
#' ## Multiple failure type examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' trt <- rbinom(n,1,0.5)
#' adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))
#' 
#' ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
#' ftype <- round(runif(n,0,2))
#' 
#' # Fit 7
#' # fit a survtmle object with glm estimators for treatment, censoring, and failure
#' # using the "mean" method
#' fit7 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean")
#' # fit7
#' 
#' # Fit 8
#' # fit an survtmle object with SuperLearner estimators for failure and censoring
#' # and empirical estimators for treatment (default) using the "mean" method
#' fit8 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="mean")
#' # fit8
#' 
#' # Fit 9
#' # repeat Fit 7 using the "hazard" method 
#' fit9 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="hazard", returnModels = TRUE)
#' # fit9
#' 
#' # Fit 10
#' # repeat Fit 8 using the "hazard" method
#' fit10 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' SL.ftime = c("SL.glm","SL.mean","SL.step"), 
#' SL.ctime = c("SL.glm","SL.mean","SL.step"), 
#' method="hazard", returnModels = TRUE)
#' # fit10 
#' 
#' # Fit 11
#' # repeat Fit 7, but only return estimates for ftype = 1
#' fit11 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", returnModels = TRUE, ftypeOfInterest = 1)
#' # fit11
#' 
#' # Fit 12 
#' # repeat Fit 7, but include bounds
#' 
#' # make a data.frame of bounds in proper format
#' # two types of failure that are labeled with ftype = 1
#' # and ftype = 2, so bounds should have columns 't','l1', 'u1',
#' # 'l2', and 'u2'.
#' bf2 <- data.frame(t=1:5,l1=rep(0.01,5),u1=rep(0.9,5),l2=rep(0.02,5),u2=rep(0.5,5))
#' 
#' fit12 <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#' glm.trt = "W1 + W2", 
#' glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
#' method="mean", returnModels = TRUE, 
#' # two failure types, so bounds should have columns 't','l1','u1','l2', and 'u2'
#' bounds = bf2)
#' # fit12


survtmle <- function(
  ftime, 
  ftype,
  trt,
  adjustVars,
  t0=max(ftime[ftype > 0]),
  SL.ftime=NULL,
  SL.ctime=NULL,
  SL.trt=NULL,
  glm.ftime=NULL,
  glm.ctime=NULL,
  glm.trt="1",
  returnIC=TRUE,
  returnModels=TRUE,
  ftypeOfInterest=unique(ftype[ftype!=0]),
  trtOfInterest=unique(trt),
  method="hazard",
  bounds=NULL,
  verbose=FALSE,
  tol=1/(length(ftime)),
  maxIter=100,
  Gcomp = FALSE,
  gtol = 1e-3
){
  
  call <- match.call(expand.dots = TRUE)
  
  # check and clean inputs
  clean <- checkInputs(ftime=ftime,ftype=ftype,trt=trt,
                       t0=t0, adjustVars=adjustVars,SL.ftime=SL.ftime,
                       SL.ctime=SL.ctime,SL.trt = SL.trt,
                       glm.ftime=glm.ftime,glm.ctime=glm.ctime,glm.trt=glm.trt,
                       returnIC=returnIC, returnModels=returnModels, 
                       ftypeOfInterest=ftypeOfInterest,trtOfInterest=trtOfInterest,
                       bounds=bounds,verbose=verbose,tol=tol,Gcomp=Gcomp,
                       method = method)

  # hazard-based TMLE
  if(method=="hazard"){
    tmle.fit <- hazard_tmle(ftime=clean$ftime, 
                            ftype=clean$ftype,
                            trt=clean$trt,
                            t0=t0,
                            adjustVars=clean$adjustVars,
                            SL.ftime=clean$SL.ftime,
                            SL.ctime=clean$SL.ctime,
                            SL.trt = clean$SL.trt,
                            glm.ftime=clean$glm.ftime,
                            glm.ctime=clean$glm.ctime,
                            glm.trt=clean$glm.trt,
                            returnIC=returnIC,
                            returnModels=returnModels,
                            ftypeOfInterest=ftypeOfInterest,
                            trtOfInterest=trtOfInterest,
                            bounds=bounds,
                            verbose=verbose, 
                            tol=tol, 
                            maxIter=maxIter, gtol = gtol)
  }else if(method=="mean"){
    tmle.fit <- mean_tmle(ftime=clean$ftime, 
                          ftype=clean$ftype,
                          trt=clean$trt,
                          t0=t0,
                          adjustVars=clean$adjustVars,
                          SL.ftime=clean$SL.ftime,
                          SL.ctime=clean$SL.ctime,
                          SL.trt = clean$SL.trt,
                          glm.ftime=clean$glm.ftime,
                          glm.ctime=clean$glm.ctime,
                          glm.trt=clean$glm.trt,
                          returnIC=returnIC,
                          returnModels=returnModels,
                          ftypeOfInterest=ftypeOfInterest,
                          trtOfInterest=trtOfInterest,
                          bounds=bounds,
                          verbose=verbose, 
                          tol=tol,
                          Gcomp=Gcomp, gtol = gtol
    )
  }
  
  out <- list(
    call=call,
    est=tmle.fit$est, var=tmle.fit$var,
    meanIC=tmle.fit$meanIC, ic=tmle.fit$ic,
    ftimeMod=tmle.fit$ftimeMod, ctimeMod=tmle.fit$ctimeMod,
    trtMod=tmle.fit$trtMod, t0=t0
  )
  class(out) <- "survtmle"
  out
}
