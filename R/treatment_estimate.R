#' SuperLearning for Multinomial Outcomes
#' 
#' @param Y 
#' @param X 
#' @param 

multi_SuperLearner <- function(Y, X, SL.library, 
    id = NULL, verbose = FALSE, control = list(), 
    cvControl = list()){
  # Work flow
  # Order outcomes by prevalence
  # Starting with most prevalent outcome, a_1, fit SuperLearner for
  #     A = a_1 ~ W 
  # Move to next most prevalent outcome, a_2, fit superLearner for
  #     A = a_2 ~ W | A != a_1
  # Map to A = a_2, multiplying probabilities 

  # examples
  n <- 2500
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rbinom(n, 1, 1/2))
  Y <- rbinom(n, 1, plogis(X$x1 + X$x2 * X$x3)) + rbinom(n, 1, plogis(0.4*X$x2^2 + X$x3))
  SL.library <- c("SL.glm","SL.gam","SL.mean","SL.earth","SL.randomForest")
  # Order outcomes by prevalence
  Y_uniq <- sort(unique(Y))
  Y_ord <- Y_uniq[rank(-table(Y))]
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
                                  SL.library = SL.library,
                                  id = id, verbose = verbose, control = control,
                                  cvControl = cvControl)
    if(ct == 1){
      learner_pred_list[[ct]] <- sl_list[[ct]]$Z      
    }else{
      learner_pred_list[[ct]] <- sl_list[[ct]]$Z * 
        Reduce("*",lapply(learner_pred_list[1:(ct-1)], function(x){ 1 - x }))
    }
  }
  # fill in last values
  learner_pred_list[[n_outcome]] <- 1 - Reduce("+", learner_pred_list[1:(n_outcome-1)])

  # now need to compute CV risk
  cvRisk <- rep(NA, n_algo)
  # stupid way to do this
  ind <- rep(n_outcome, n)
  for(j in 1:(n_outcome-1)){
    ind[Y == Y_ord[j]] <- j
  }

  # this will hold the fit predicted probability of each algorithm
  # for the outcome that was observed
  obs_outcome_pred <- matrix(NA, nrow = n, ncol = n_algo)
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
  pred <- matrix(NA, nrow = n, ncol = n_outcome)
  for(i in 1:(n_outcome-1)){
    pred[,i] <- plogis(SuperLearner::trimLogit(sl_list[[i]]$library.predict) %*% 
                       matrix(coef))
  }
  pred[,n_outcome] <- 1 - rowSums(pred[,1:(n_outcome-1)])
  colnames(pred) <- Y_ord
  predY <- apply(pred, 1, function(x){
    Y_ord[which.max(x)]
  })

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

  if(length(unique(dat$trt)) == 1) {
    eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
  } else {
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
  }

  # truncate propensities
  eval(parse(text = paste0("dat$g_", min(dat$trt), "[dat$g_", min(dat$trt),
                           "< gtol]<- gtol")))
  eval(parse(text = paste0("dat$g_", max(dat$trt), "[dat$g_", max(dat$trt),
                           "< gtol]<- gtol")))
  out <- list()
  out$dat <- dat
  out$trtMod <- NULL
  if(returnModels) out$trtMod <- trtMod
  return(out)
}
