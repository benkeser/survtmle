#' SuperLearning for Multinomial Outcomes
#' 
#' This function is a bare-bones implementation of super learning for a
#' multinomial outcome. It generates compatible ensemble estimates of each class 
#' of an outcome by fitting a sequence of conditional binomial regressions.
#' For example, if there are three class of outcome, two prediction functions are fit:
#' the first regresses an indicator of the class with most observations onto 
#' the predictors, the second regresses an indicator of the class with second most 
#' observations onto the predictors in the subset of data not in class one. The second
#' class probabilities are generated as the probability of being in the second class given
#' not in the first times the probability of not being in the first. The third class
#' probabilities are generated as one minus the sum of the other two probabilities. 
#' 
#' Each candidate learner specified in \code{SL.library} is fit this way with a separate
#' call to \code{SuperLearner} with \code{family = binomial()}. The ensemble is generated
#' as the convex combination of learners that maximizes cross-validated log-likelihood for
#' the multinomial outcome. 
#' 
#' @param Y The multinomial outcome. Must be a numeric vector.
#' @param X The predictor variables, usually a data.frame.
#' @param verbose Print messages 
#' @param SL.library A character vector of prediction algorithms. 
#' See details below for examples on the structure. 
#' @param control A list of parameters to control the estimation process. 
#' Parameters include saveFitLibrary and trimLogit. 
#' See SuperLearner::SuperLearner.control for details.
#' @param cvControl A list of parameters to control the cross-validation process. 
#' Parameters include V, stratifyCV, shuffle and validRows. 
#' See SuperLearner::SuperLearner.CV.control for details.
#' 
#' @return An object of class \code{multiSuperLearner}. \describe{
#' \item{SuperLearnerFits}{The conditional binomial SuperLearner fits}
#' \item{coef}{The super learner ensemble weights}
#' \item{SL.predict}{A matrix with columns corresponding to super learner predictions
#' of each of the outcomes based on \code{X}.}
#' \item{library.predict}{A list with each entry corresponding to the candidate 
#' learner fits in each of the conditional binomial regressions.}
#' \item{cvRisk}{The cross-validated negative log-multinomial likelihood loss for 
#' each candidate learner.}
#' }
#' 
#' @export
#' 
#' @examples
#' # simulate data
#' set.seed(1234)
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' Y <- rbinom(100, 1, 0.2) + rbinom(100, 1, plogis(X$x1)) + 
#'          rbinom(100, 1, plogis(X$x2^2))
#' 
#' # fit super learner
#' fit <- multiSuperLearner(Y = Y, X = X, SL.library = c("SL.glm","SL.mean"))
#' 

multiSuperLearner <- function(Y, X, SL.library, 
    verbose = FALSE, control = list(), 
    cvControl = list(), ...){
  # Work flow
  # Order outcomes by prevalence
  # Starting with most prevalent outcome, a_1, fit SuperLearner for
  #     A = a_1 ~ W 
  # Move to next most prevalent outcome, a_2, fit superLearner for
  #     A = a_2 ~ W | A != a_1
  # Map to A = a_2, multiplying probabilities 

  # Order outcomes by prevalence
  tabY <- table(Y)
  uniqY <- unique(Y)
  Y_uniq <- sort(uniqY)
  orderY <- order(-tabY)
  rankY <- rank(-tabY)
  Y_ord <- Y_uniq[orderY]
  n_outcome <- length(Y_ord)
  n_algo <- length(SL.library)
  n <- length(Y)
  sl_list <- vector(mode = "list", length = n_outcome - 1)
  learner_pred_list <- vector(mode = "list", length = n_outcome)
  ct <- 0
  # loop over unique values of Y
  # other than the last, which gets 1 - the rest
  for(y in Y_ord[-length(Y_ord)]){
    ct <- ct + 1
    # who to include
    include <- rep(TRUE, n) 
    if(ct > 1){
      include[Y %in% Y_ord[1:(ct-1)]] <- FALSE
    }
    # fit SuperLearner for Y == y ~ X
    sl_list[[ct]] <- SuperLearner(Y = as.numeric(Y == y), 
                                  X = X,
                                  family = binomial(),
                                  SL.library = SL.library, 
                                  verbose = verbose, control = control,
                                  cvControl = cvControl)
    if(ct == 1){
      learner_pred_list[[ct]] <- sl_list[[ct]]$Z      
    }else{
      learner_pred_list[[ct]] <- sl_list[[ct]]$Z * 
        Reduce("*",lapply(learner_pred_list[1:(ct-1)], function(x){ 1 - x }))
    }
    names(learner_pred_list[[ct]]) <- SL.library
  }
  # fill in last values
  learner_pred_list[[n_outcome]] <- 1 - Reduce("+", learner_pred_list[1:(n_outcome-1)])

  # name sl_list
  names(sl_list) <- paste0("y = ", Y_ord[-length(Y_ord)])
  names(learner_pred_list) <- paste0("y = ", Y_ord[-length(Y_ord)])
  # now need to compute CV risk
  # stupid way to do this
  ind <- rep(n_outcome, n)
  for(j in 1:(n_outcome-1)){
    ind[Y == Y_ord[j]] <- j
  }

  # this will hold the fit predicted probability of each algorithm
  # for the outcome that was observed
  obs_outcome_pred <- matrix(NA, nrow = n, ncol = n_algo)
  cvRisk <- rep(NA, n_algo)
  for(i in 1:n_algo){
    # get this algorithm from the fit for each outcome and
    # take the trimmed logit
    logit_thisZ <- Reduce("cbind",lapply(learner_pred_list, function(x){ 
      SuperLearner::trimLogit(x[, i]) }))
    # now get a single vector, where we take the proper column
    # for the observed value of the outcome
    obs_outcome_pred[,i] <- logit_thisZ[cbind(seq_along(ind),ind)]
    cvRisk[i] <- -sum(plogis(obs_outcome_pred[,i], log.p = TRUE))
  }

  # now need to compute super learner
  # we want x to come in already linked up to the appropriate y
  obj_and_grad <- function(x) {
    function(beta) {
        xB <- x %*% cbind(beta)
        obj <- -sum(plogis(xB, log.p = TRUE))
        p <- plogis(xB)
        grad <- crossprod(x, cbind(-1 / (1 + exp(xB))))
        list(objective = obj, gradient = grad)
    }
  }
  lower_bounds = rep(0, n_algo)
  upper_bounds = rep(1, n_algo)
  if (anyNA(cvRisk)) {
      upper_bounds[is.na(cvRisk)] = 0
  }
  r <- nloptr::nloptr(x0 = rep(1/n_algo, n_algo), 
                      eval_f = obj_and_grad(obs_outcome_pred), 
                      lb = lower_bounds, ub = upper_bounds, 
                      eval_g_eq = function(beta) (sum(beta) - 1), 
                      eval_jac_g_eq = function(beta) rep(1, length(beta)), 
                      opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
  coef <- r$solution
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] <- 0
        }
  coef[coef < 1e-04] <- 0
  coef <- coef/sum(coef)

  # compute super learner predictions
  # stupid way to do it 
  pred <- vector(mode = "list", length = n_outcome)
  for(i in 1:(n_outcome-1)){
    if(i == 1){
      learner_mat <- sl_list[[i]]$library.predict
    }else{
      learner_mat <- sl_list[[i]]$library.predict * 
        Reduce("*",lapply(sl_list[1:(i-1)], function(x){ 1 - x$library.predict }))
    }
    pred[[i]] <- plogis(SuperLearner::trimLogit(learner_mat) %*% matrix(coef))
  }
  pred[[n_outcome]] <- 1 - Reduce("+", pred[1:(n_outcome - 1)])

  # reorder columns
  pred <- pred[rankY]
  names(pred) <- paste0("y = ", Y_uniq)
  # predY <- apply(pred, 1, function(x){
  #   Y_ord[which.max(x)]
  # })

  # format output
  out <- list()
  out$SuperLearnerFits <- sl_list
  out$coef <- coef
  out$SL.predict <- pred
  out$library.predict <- learner_pred_list
  out$cvRisk <- cvRisk

  class(out) <- "multiSuperLearner"
  return(out)
}

#' Estimate Treatment Mechanisms
#'
#' This function computes the conditional probability of having \code{trt} for
#' each specified level either using \code{glm} or \code{SuperLearner}.
#' Currently only two unique values of treatment are acceptable. By default the
#' function will compute estimates of the conditional probability of
#' \code{trt == max(trt)} and compute the probability of \code{trt == min(trt)}
#' as one minus this probability.
#'
#' @param dat An object of class \code{data.frame}. Must have named column
#'        \code{trt}.
#' @param adjustVars An object of class \code{data.frame} that will be used
#'        either as the \code{data} argument in a call to \code{glm} or as the
#'        \code{X} object in a call to \code{SuperLearner}.
#' @param glm.trt A character formula for the right-hand side of \code{formula}
#'        in a call to \code{glm}. See \code{?survtmle} for more documentation.
#'        Alternatively, this could be an object of class \code{glm} (as in
#'        calls to this function via \code{timepoints}), in which case
#'        predictions are obtained using this object with no new fitting.
#' @param SL.trt A valid specification of the \code{SL.library} option of a call
#'        to \code{SuperLearner}. See \code{?survtmle} for more documentation.
#'        Alternatively, this could be an object of class \code{SuperLearner}
#'        (as in calls to this function via \code{timepoints}), in which case
#'        predictions are obtained using this object with no new fitting.
#' @param returnModels A boolean indicating whether fitted model objects should
#'        be returned.
#' @param verbose A boolean passed to the \code{verbose} option of the call to
#'        \code{SuperLearner}.
#' @param gtol The truncation level of predicted trt probabilities to handle
#'        positivity violations.
#' @param ... Other arguments. Not currently used
#'
#' @return dat The input \code{data.frame} object with two added columns
#'         corresponding with the conditional probability (given
#'         \code{adjustVars}) of \code{trt==max(trt)} and \code{trt==min(trt)}.
#' @return trtMod If \code{returnModels = TRUE}, the fitted \code{glm} or
#'         \code{SuperLearner} object. Otherwise, \code{NULL}
#'
#' @importFrom stats as.formula predict model.matrix optim glm
#' @importFrom SuperLearner SuperLearner SuperLearner.CV.control All SL.mean SL.glm SL.step
#'

estimateTreatment <- function(dat, adjustVars, glm.trt = NULL, SL.trt = NULL,
                              returnModels = FALSE, verbose = FALSE,
                              gtol = 1e-3, ...) {
  n <- length(dat[,1])
  uniqtrt <- unique(dat$trt)
  ntrt <- length(uniqtrt)

  if(length(unique(dat$trt)) == 1) {
    eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
  } else if(length(unique(dat$trt)) == 2) {
    # binarize the outcome
    thisY <- as.numeric(dat$trt == max(dat$trt))

    # fit Super Learner
    if(!is.null(SL.trt)) {
      if(class(SL.trt) != "SuperLearner") {
        trtMod <- SuperLearner::SuperLearner(Y = thisY, X = adjustVars,
                                             newX = adjustVars,
                                             SL.library = SL.trt,
                                             id = dat$id, verbose = verbose,
                                             family = "binomial")
      } else {
        trtMod <- SL.trt
      }
      dat[[paste0("g_",max(dat$trt))]] <- trtMod$SL.predict
      dat[[paste0("g_",min(dat$trt))]] <- 1 - trtMod$SL.predict

    } else if(!is.null(glm.trt) & is.null(SL.trt)) {
      # set up model formula and data for the treatment regression
      trt_form <- paste("thisY", "~", glm.trt, sep = " ")
      trt_data_in <- as.data.frame(cbind(adjustVars, thisY))

      # fit GLM if Super Learner not requested
      if(!("glm" %in% class(glm.trt)) & !("speedglm" %in% class(glm.trt))) {
        # fit the treatment model
        trtMod <- fast_glm(reg_form = stats::as.formula(trt_form),
                           data = trt_data_in,
                           family = stats::binomial())
      } else {
        trtMod <- glm.trt
      }
      suppressWarnings(
        pred <- predict(trtMod, newdata = trt_data_in, type = "response")
      )
      dat[[paste0("g_", max(dat$trt))]] <- pred
      dat[[paste0("g_", min(dat$trt))]] <- 1 - pred
    }
  } else {
    # fit Super Learner
    if(!is.null(SL.trt)) {
      if(class(SL.trt) != "multi_SuperLearner") {
        trtMod <- multiSuperLearner(Y = dat$trt, X = adjustVars,
                                     SL.library = SL.trt,
                                     verbose = verbose)
      } else {
        trtMod <- SL.trt
      }
      for(a in unique(dat$trt)){
        dat[[paste0("g_",a)]] <- trtMod$SL.predict[[paste0("y = ",a)]]
      }
    } else if(!is.null(glm.trt) & is.null(SL.trt)) {
      ct <- 0
      trtMod <- vector(mode = "list", length = length(unique(dat$trt))-1)
      pred <- vector(mode = "list", length = length(unique(dat$trt))-1)
      sortuniqtrt <- sort(unique(dat$trt))
      for(a in sortuniqtrt[-length(sortuniqtrt)]){
        ct <- ct + 1
        # set up model formula and data for the treatment regression
        trt_form <- paste("thisY", "~", glm.trt, sep = " ")
        trt_data_in <- data.frame(adjustVars, thisY = as.numeric(dat$trt == a))
        # fit GLM if Super Learner not requested
        if(!("glm" %in% class(glm.trt[[1]])) & !("speedglm" %in% class(glm.trt[[1]]))) {
          # fit the treatment model
          trtMod[[ct]] <- fast_glm(reg_form = stats::as.formula(trt_form),
                                   data = trt_data_in[dat$trt >= a, ],
                                   family = stats::binomial())
        } else {
          trtMod[[ct]] <- glm.trt[[ct]]
        }    
        if(ct == 1){
          suppressWarnings(
            pred[[ct]] <- predict(trtMod[[ct]], newdata = trt_data_in, type = "response")
          )
        }else{
          pred[[ct]] <- predict(trtMod[[ct]], newdata = trt_data_in, type = "response") * 
            (1 - Reduce("+",lapply(pred[1:(ct-1)], function(x){ x })))
        }
      }
      for(i in 1:(length(sortuniqtrt)-1)){
        dat[[paste0("g_", sortuniqtrt[i])]] <- pred[[i]]
      }
      dat[[paste0("g_", sortuniqtrt[length(sortuniqtrt)])]] <- 1 - Reduce("+", pred)
    }

  }

  # truncate propensities
  for(a in unique(dat$trt)){
    eval(parse(text = paste0("dat$g_", a, "[dat$g_", a,
                           "< gtol]<- gtol")))  
  }

  # make a column of observed a
  ind <- rep(NA, n)
  for(j in 1:ntrt){
    ind[dat$trt == uniqtrt[j]] <- which(colnames(dat) == paste0("g_",uniqtrt[j]))
  }
  dat$g_obsz <- dat[cbind(seq_along(ind),ind)]

  out <- list()
  out$dat <- dat
  out$trtMod <- NULL
  if(returnModels) out$trtMod <- trtMod
  return(out)
}
