#' Extract Influence Curve for Estimated Hazard Functions
#'
#' @description This function computes the hazard-based efficient influence
#'  curve at the final estimate of the fluctuated cause-specific hazard
#'  functions and evaluates it on the observed data. The influence function is
#'  computed on the long-format data but is subsequently summed over all
#'  timepoints for each observation and the function returns a new short-form
#'  data set with columns added corresponding to the sum over all timepoints of
#'  the estimated efficient influence function evaluated at that observation.
#'
#' @param dataList A list of \code{data.frame} objects. See the documentation
#'  of \code{\link{makeDataList}} for more information.
#' @param dat A \code{data.frame} in short form. See the documentation of
#'  \code{\link{makeDataList}} for more information.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param ofInterestJ Numeric vector indicating \code{ftypeOfInterest} that was
#'  passed to \code{\link{hazard_tmle}}.
#' @param nJ The number of unique failure types.
#' @param uniqtrt The values of \code{trtOfInterest} passed to
#'  \code{\link{mean_tmle}}.
#' @param t0 The timepoint at which \code{\link{survtmle}} was called to
#'  evaluate.
#' @param verbose A \code{logical} indicating whether the function should print
#'  messages to indicate progress.
#' @param att A \code{boolean} indicating whether to compute the ATT estimate,
#'  instead of treatment specific survival curves. This option only works with 
#'  two levels of \code{trt} that are labeled with 0 and 1.
#' @param ... Other arguments. Not currently used.
#' @return An object of class \code{data.frame} with columns \code{D.jX.zZ}
#'  added for each value of X in \code{ofInterestJ} and each value of Z in
#'  \code{uniqtrt}. These are the sum over all timepoints of the estimated
#'  efficient influence function evaluated at that observation.
getHazardInfluenceCurve <- function(dataList, dat, allJ, ofInterestJ, nJ,
                                    uniqtrt, t0, verbose, att, ...) {
  na_ftime_idx <- which(is.na(dat$ftime))
  this_list_idx <- 1
  for (z in uniqtrt) {
    this_list_idx <- this_list_idx + 1
    for (j in ofInterestJ) {
      dat[[paste0("F", j, ".z", z, ".t0")]] <- rep(NA, dim(dat)[1])
      dat[[paste0("F", j, ".z", z, ".t0")]][!is.na(dat$ftime)] <-
        dataList[[1]][[paste0("F", j, ".z", z, ".t0")]][dataList[[1]]$t == min(dataList[[1]]$t)]

      if(any(is.na(dat$ftime))){
        dat[[paste0("F", j, ".z", z, ".t0")]][na_ftime_idx] <- 
          dataList[[this_list_idx]][[paste0("F", j, ".z", z, ".t0")]][
              dataList[[this_list_idx]]$id %in% na_ftime_idx & ( dataList[[this_list_idx]]$t == min(dataList[[this_list_idx]]$t) )
            ]
      }
      if(!att){
        dat[[paste0("margF", j, ".z", z, ".t0")]] <- mean(dat[[paste0("F", j, ".z", z, ".t0")]])
      }else{
        dat[[paste0("margF", j, ".z", z, ".t0")]] <- mean(dat[[paste0("F", j, ".z", z, ".t0")]][dat$trt == 1])
      }

      # non-NA people
      thisD <- NULL
      for (jTild in allJ) {
        H <- paste0("H", j, ".j", ifelse(jTild == j, "Self", "NotSelf"), ".z", z)
        thisD <- cbind(thisD, dataList[[1]][[H]] / (1 - dataList[[1]][[paste0("hazNot", j)]]) *
          (dataList[[1]][[paste0("N", jTild)]] - dataList[[1]][[paste0("Q", jTild, "Haz")]]))
      }
      if(length(na_ftime_idx) > 0){
        tmp_D <- rep(0, dim(dat)[1])
        tmp_D[-na_ftime_idx] <- c(by(rowSums(thisD), dataList[[1]]$id, FUN = sum))
      }else{
        tmp_D <- c(by(rowSums(thisD), dataList[[1]]$id, FUN = sum))
      }
      
      if(!att){
        dat[[paste0("D.j", j, ".z", z)]] <- tmp_D + 
          dat[[paste0("F", j, ".z", z, ".t0")]] - dat[[paste0("margF", j, ".z", z, ".t0")]]
      }else{
        tmp_D <- tmp_D / mean(dat$trt == 1)
        dat[[paste0("D.j", j, ".z", z)]] <- tmp_D + 
          as.numeric(dat$trt == 1) / mean(dat$trt == 1) * 
            (dat[[paste0("F", j, ".z", z, ".t0")]] - dat[[paste0("margF", j, ".z", z, ".t0")]])
      }
    }
  }
  return(dat)
}
