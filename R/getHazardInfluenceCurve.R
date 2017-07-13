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
#' @param ... Other arguments. Not currently used.
#'
#' @return An object of class \code{data.frame} with columns \code{D.jX.zZ}
#'         added for each value of X in \code{ofInterestJ} and each value of Z
#'         in \code{uniqtrt}. These are the sum over all timepoints of the
#'         estimated efficient influence function evaluated at that observation.
#'

getHazardInfluenceCurve <- function(dataList, dat, allJ, ofInterestJ, nJ,
                                    uniqtrt, t0, verbose, ...) {
  for(z in uniqtrt) {
    for(j in ofInterestJ) {
      dat[[paste0("margF",j,".z",z,".t0")]] <- 
        mean(dataList[[1]][[paste0("F",j,".z",z,".t0")]][dataList[[1]]$t==min(dataList[[1]]$t)])

      dat[[paste0("F",j,".z",z,".t0")]] <- 
        dataList[[1]][[paste0("F",j,".z",z,".t0")]][dataList[[1]]$t==min(dataList[[1]]$t)]
      thisD <- NULL
      for(jTild in allJ) {
        H <- paste0("H",j,".j",ifelse(jTild==j,"Self","NotSelf"),".z",z)
        thisD <- cbind(thisD, dataList[[1]][[H]]/(1-dataList[[1]][[paste0("hazNot",j)]])*
                      (dataList[[1]][[paste0("N",jTild)]] - dataList[[1]][[paste0("Q",jTild,"Haz")]]))
      }
      dat[[paste0("D.j",j,".z",z)]] <- unlist(by(rowSums(thisD), dataList[[1]]$id, FUN=sum)) +
          dat[[paste0("F",j,".z",z,".t0")]] - dat[[paste0("margF",j,".z",z,".t0")]]
    }
  }
  return(dat)
}
