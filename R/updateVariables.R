#' Update TMLEs for Hazard to Cumulative Incidence
#'
#' A helper function that maps hazard estimates into estimates of cumulative
#' incidence and updates the "clever covariates" used by the targeted minimum
#' loss-based estimation fluctuation step.
#'
#' @param dataList A list of \code{data.frame} objects.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param ofInterestJ Numeric vector indicating \code{ftypeOfInterest} that was
#'        passed to \code{hazard_tmle}.
#' @param nJ The number of unique failure types.
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean_tmle}.
#' @param ntrt The number of \code{trt} values of interest.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param ... Other arguments. Not currently used.
#'
#' @return The function returns a list that is exactly the same as the input
#'         \code{dataList}, but with updated columns corresponding with
#'         estimated cumulative incidence at each time and estimated "clever
#'         covariates" at each time.
#'

updateVariables <- function(dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt, t0,
                            verbose, ...) {
  dataList[2:(ntrt + 1)] <- lapply(dataList[2:(ntrt + 1)], function(x, allJ) {
    # total hazard
    Q_dot <- rowSums(cbind(rep(0, nrow(x)), x[, paste0("Q", allJ, "Haz")]))
    # survival at t
    x$S.t <- unlist(by((1 - Q_dot), x$id, FUN = cumprod))
    x$S.t[x$S.t == 0] <- .Machine$double.neg.eps
    # survival at t-1
    S.tminus1 <- c(1, x$S.t[1:(length(x$S.t) - 1)])
    S.tminus1[x$t == 1] <- 1

    for (j in ofInterestJ) {
      # calculate CIF at time t
      x[[paste0("F", j, ".t")]] <- unlist(by(x[, paste0("Q", j, "Haz")] * S.tminus1, x$id, FUN = cumsum))
    }
    x
  }, allJ = allJ)

  # calculate CIF at time t0
  for (j in ofInterestJ) {
    # TO DO: change this to static memory allocation
    Fj.t0.allZ <- vector(mode = "list", length = ntrt)
    for (i in 1:ntrt) {
      t0.mod <- dataList[[i + 1]]$ftime[1]
      Fj.t0.allZ[[i]] <- dataList[[i + 1]][[paste0("F", j, ".t")]][dataList[[i + 1]]$t == t0.mod]
    }

    dataList <- lapply(dataList, function(x, j, uniqtrt, Fj.t0.allZ) {
      for (i in seq_along(uniqtrt)) {
        ind <- tapply(X = x$id, INDEX = x$id, FUN = NULL)
        x[[paste0("F", j, ".z", uniqtrt[i], ".t0")]] <- Fj.t0.allZ[[i]][ind]
      }
      x
    }, j = j, uniqtrt = uniqtrt, Fj.t0.allZ = Fj.t0.allZ)
  }

  # merge into dataList[[1]]
  # indicators of S.t and Fj.t for dataList[[1]]
  colInd <- which(colnames(dataList[[1]]) %in% c("S.t", paste0(
    "F",
    ofInterestJ,
    ".t"
  )))
  # the first time it's called these columns won't exist
  if (length(colInd) == 0) {
    dataList[[1]] <- merge(
      dataList[[1]],
      Reduce(
        rbind,
        dataList[2:(ntrt + 1)]
      )[, c(
        "id", "t", "trt",
        "S.t",
        paste0(
          "F", ofInterestJ,
          ".t"
        )
      )],
      by = c("id", "t", "trt")
    )
  } else {
    # the next times it's called those columns will exist but we want them replaced
    # with the values from dataList[[>1]]
    dataList[[1]] <- merge(
      dataList[[1]][, -colInd],
      Reduce(
        rbind,
        dataList[2:(ntrt + 1)]
      )[, c(
        "id", "t", "trt",
        "S.t",
        paste0(
          "F",
          ofInterestJ,
          ".t"
        )
      )],
      by = c("id", "t", "trt")
    )
  }

  dataList <- lapply(dataList, function(x, allJ) {
    for (j in allJ) {
      if (length(allJ) > 1) {
        x[[paste0("hazNot", j)]] <- rowSums(cbind(rep(0, nrow(x)), x[, paste0("Q", allJ[allJ != j], "Haz")]))
        x[[paste0("hazNot", j)]][x[[paste0("hazNot", j)]] == 1] <- 1 - .Machine$double.neg.eps
      } else {
        x[[paste0("hazNot", j)]] <- 0
      }
    }
    x
  }, allJ = allJ)


  # set up clever covariates needed for fluctuation
  dataList <- lapply(dataList, function(x, ofInterestJ, uniqtrt) {
    for (z in uniqtrt) {
      for (j in ofInterestJ) {
        x[[paste0("H", j, ".jSelf.z", z)]] <-
          (x$ftime >= x$t & x$trt == z) / (x[[paste0("g_", z)]] * x$G_dC) *
            (1 - x[[paste0("hazNot", j)]]) * ((x$t < t0) * (1 - (x[[paste0("F", j, ".z", z, ".t0")]] -
              x[[paste0("F", j, ".t")]]) / c(x$S.t)) + as.numeric(x$t == t0))
        x[[paste0("H", j, ".jNotSelf.z", z)]] <-
          -(x$ftime >= x$t & x$trt == z) / (x[[paste0("g_", z)]] * x$G_dC) *
            (1 - x[[paste0("hazNot", j)]]) * ((x$t < t0) * (x[[paste0("F", j, ".z", z, ".t0")]] -
              x[[paste0("F", j, ".t")]]) / c(x$S.t))
      }
    }
    x
  }, ofInterestJ = ofInterestJ, uniqtrt = uniqtrt)
  #  } else {
  #    dataList <- lapply(dataList, function(x, ofInterestJ, uniqtrt) {
  #    # placebo match
  #      x$H1.jSelf.z0 <- 1 / x$F1.z0.t0 * (x$ftime >= x$t & x$trt == 0) /
  #        (x$g_0 * x$G_dC) * (1 - x$hazNot1) *
  #        ((x$t < t0) * (1 - (x$F1.z0.t0 - x$F1.t)/c(x$S.t)) + (x$t==t0))
  #      x$H1.jNotSelf.z0 <- 1 / x$F1.z0.t0 * - (x$ftime >= x$t & x$trt == 0) /
  #        (x$g_0 * x$G_dC) *(1 - x$hazNot1) * ((x$t < t0) *
  #                                             (x$F1.z0.t0 - x$F1.t) / c(x$S.t))
  #    # vaccine match
  #      x$H1.jSelf.z1 <- -1 / x$F1.z1.t0 * (x$ftime >= x$t & x$trt == 1) /
  #        (x$g_1 * x$G_dC) * (1 - x$hazNot1) *
  #        ((x$t < t0) * (1 - (x$F1.z1.t0 - x$F1.t) / c(x$S.t)) + (x$t == t0))
  #      x$H1.jNotSelf.z1 <- -1 / x$F1.z1.t0 * -(x$ftime >= x$t & x$trt == 1) /
  #        (x$g_1*x$G_dC) * (1 - x$hazNot1) *
  #        ((x$t < t0) * (x$F1.z1.t0 - x$F1.t) / c(x$S.t))
  #    # placebo mismatch
  #      x$H2.jSelf.z0 <- -1 / x$F2.z0.t0 * (x$ftime >= x$t & x$trt == 0) /
  #        (x$g_0 * x$G_dC) * (1 - x$hazNot2) *
  #        ((x$t < t0) * (1 - (x$F2.z0.t0 - x$F2.t) / c(x$S.t)) + (x$t == t0))
  #      x$H2.jNotSelf.z0 <- -1 / x$F2.z0.t0 * -(x$ftime >= x$t & x$trt == 0) /
  #        (x$g_0 * x$G_dC) *(1 - x$hazNot2) *
  #        ((x$t < t0)*(x$F2.z0.t0 - x$F2.t) / c(x$S.t))
  #    # vaccine mismatch
  #      x$H2.jSelf.z1 <- 1 / x$F2.z1.t0 * (x$ftime >= x$t & x$trt == 1) /
  #        (x$g_1 * x$G_dC) * (1 - x$hazNot2) *
  #        ((x$t < t0) * (1 - (x$F2.z1.t0 - x$F2.t)/c(x$S.t)) + (x$t == t0))
  #      x$H2.jNotSelf.z1 <- 1 / x$F2.z1.t0 * -(x$ftime >= x$t & x$trt == 1) /
  #        (x$g_1 * x$G_dC) * (1 - x$hazNot2) *
  #        ((x$t < t0) * (x$F2.z1.t0 - x$F2.t)/c(x$S.t))
  #      x
  #    })
  #  }
  dataList
}
