#' Extract Influence Curve for Estimated Hazard Functions
#'
#' This function computes the hazard-based efficient influence curve at the
#' final estimate of the fluctuated cause-specific hazard functions and
#' evaluates it on the observed data. The influence-function is computed on the
#' long-format data but is subsequently summed over all timepoints for each
#' observation and the function returns a new short form data set with columns
#' added corresponding to the sum over all timepoints of the estimated
#' efficient influence function evaluated at that observation.
#'
#' @param dataList A list of \code{data.frame} objects. See \code{?makeDataList}
#'        for more information.
#' @param dat A \code{data.frame} in short form. See \code{?makeDataList} for
#'        more information.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param ofInterestJ Numeric vector indicating \code{ftypeOfInterest} that was
#'        passed to \code{hazard_tmle}.
#' @param nJ The number of unique failure types.
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean_tmle}.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param msm.formula TO DO: Add documentation here
#' @param msm.family
#' @param msm.est
#' @param ... Other arguments. Not currently used.
#'
#' @return An object of class \code{data.frame} with columns \code{D.jX.zZ}
#'         added for each value of X in \code{ofInterestJ} and each value of Z
#'         in \code{uniqtrt}. These are the sum over all timepoints of the
#'         estimated efficient influence function evaluated at that observation.
#'

getHazardInfluenceCurve <- function(dataList, dat, allJ, ofInterestJ, nJ,
                                    uniqtrt, t0, verbose, msm.formula, 
                                    msm.family, msm.est = NULL,
                                    msmWeightList, ...) {
  n <- length(dat[,1])
  if(is.null(msm.formula)){
    for(z in uniqtrt) {
      for(j in ofInterestJ) {
        dat[[paste0("margF",j,".z",z,".t0")]] <- 
          mean(dataList[[1]][[paste0("F",j,".z",z,".t0")]][dataList[[1]]$t==min(dataList[[1]]$t)])

        dat[[paste0("F",j,".z",z,".t0")]] <- 
          dataList[[1]][[paste0("F",j,".z",z,".t0")]][dataList[[1]]$t==min(dataList[[1]]$t)]
        thisD <- NULL
        for(jTild in allJ) {
          H <- paste0("H",j,".j",ifelse(jTild==j,"Self","NotSelf"),".z",z)
          thisD <- cbind(thisD, dataList[[1]][[H]] * (dataList[[1]][[paste0("N",jTild)]] - dataList[[1]][[paste0("Q",jTild,"Haz")]]))
        }
        dat[[paste0("D.j",j,".z",z)]] <- unlist(by(rowSums(thisD), dataList[[1]]$id, FUN=sum)) +
            dat[[paste0("F",j,".z",z,".t0")]] - dat[[paste0("margF",j,".z",z,".t0")]]
      }
    }
  }else{
    # compute the msm influence curve...
    outcomeName <- paste0("F",ofInterestJ,".z.t0")
    min.t <- min(dataList[[2]]$t)
    outcomeList <- lapply(dataList[2:length(dataList)], function(l){
      l[l$t == min.t, paste0("F",ofInterestJ,".z.t0")]
    })
    modelMatrixList <- lapply(dataList[2:length(dataList)], function(dl){
      model.matrix(as.formula(paste0(outcomeName, "~", msm.formula)), data = dl[dl$t == min.t, , drop = FALSE])
    })
    modelMatrixObs <- model.matrix(as.formula(paste0(outcomeName, "~", msm.formula)), 
                                   data = dataList[[1]][dataList[[1]]$t == min.t, , drop = FALSE])
    if(msm.family == "binomial"){
      msm.family <- binomial()
    }else if(msm.family == "gaussian"){
      msm.family <- gaussian()
    }
    msm.p <- length(msm.est)
    if(msm.family$family == "gaussian"){
      fittedValueList <- lapply(modelMatrixList, function(mm){
        mm %*% matrix(msm.est)
      })
      derivList <- lapply(modelMatrixList, function(mm){
        rep(1, n)
      })
      fittedValues <- modelMatrixObs %*% matrix(msm.est)
      deriv <- rep(1, n)
    }else if(msm.family$family == "binomial"){
      fittedValueList <- lapply(modelMatrixList, function(mm){
        stats::plogis(mm %*% matrix(msm.est))
      })
      derivList <- lapply(fittedValueList, function(fitted_value){
        fitted_value * (1 - fitted_value)
      })
      fittedValues <- plogis(modelMatrixObs %*% matrix(msm.est))
      deriv <- fittedValues * (1 - fittedValues)
    }
    # CHECK THAT THIS CONTAINS THE CORRECT ENTRY?
    msmWeightList_t1 <- lapply(msmWeightList[3:length(msmWeightList)], function(x){
      x[dataList[[2]]$t == min.t]
    })
    cQ <- Reduce("+", mapply(mw = msmWeightList_t1, 
       mm = modelMatrixList, dl = derivList, FUN = function(mw, mm, dl){
          Reduce("+", lapply(apply(cbind(mw, dl, mm), 1, function(x){ 
            list(x[1] * x[2] * tcrossprod(matrix(x[3:(msm.p+2)])))
          }),"[[",1))
       }, SIMPLIFY = FALSE)) / n
    cQ_inv <- tryCatch(solve(cQ), error = function(){ MASS::ginv(cQ) })
    sumD1_allZ <- Reduce("+", mapply(mw = msmWeightList_t1, fv = fittedValueList, 
       mm = modelMatrixList, ol = outcomeList, function(mm, mw, fv, ol){
          oneZ <- apply(cbind(mw, fv, ol, mm), 1, function(x){
            # weight * (outcome - fitted) * model matrix
            x[1] * (x[3] - x[2]) * x[4:(msm.p+3)]
          })
          return(oneZ)
       }, SIMPLIFY = FALSE))
    # D1 is a matrix with t0 * n columns and msm.p rows 
    D1 <- cQ_inv %*% sumD1_allZ

    # other pieces
    # This seems to be the score equation that is being solved by the targeting.
    # However, it's not giving sane standard error estimates, so perhaps it's the
    # wrong score equation.
    # IDEAS: Check how H is being constructed?
    #        Check cQ 
    #        ** Check who is getting too big IC contributions and what piece of the IC
    #        is driving those too-big values ** 
    thisD <- matrix(0, ncol = msm.p, nrow = length(dataList[[1]][,1]))
    ct <- 0
    for(jTild in allJ) {
      ct <- ct + 1
      H <- paste0("H",ofInterestJ,".j",ifelse(jTild == ofInterestJ, "Self", "NotSelf"), ".", 1:msm.p)
      thisD <- thisD + as.matrix(dataList[[1]][,H]) * as.numeric(dataList[[1]][,paste0("N",jTild)] - dataList[[1]][, paste0("Q",jTild,"Haz")])
    }
    # sum up over each participants times
    tmp <- Reduce("rbind", by(thisD, dataList[[1]]$id, FUN=colSums))

    dat[,paste0("D.j",1:msm.p)] <- t(D1 + tcrossprod(cQ_inv, tmp))
    # Check whether this is correct...
  }
  return(dat)
}
