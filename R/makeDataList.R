#' Convert Short Form Data to List of Wide Form Data
#'
#' The function takes a \code{data.frame} of short format right-censored failure
#' times and reshapes the long format into the wide format needed for calls to
#' both \code{mean_tmle} and \code{hazard_tmle}. The list returned will have a
#' number of entries equal to \code{length(trtOfInterest) + 1}. The first will
#' have number of rows for each observation equal to the \code{ftime}
#' corresponding to that observation. The subsequent entries will have \code{t0}
#' rows for each observation and will set \code{trt} column equal to each value
#' of \code{trtOfInterest} in turn.
#'
#' @param dat The short form \code{data.frame}
#' @param J The unique values of \code{ftype} passed to \code{survtmle}.
#' @param ntrt The number of \code{trt} values of interest.
#' @param uniqtrt The unique values of \code{trtOfInterest} passed to
#'        \code{mean_tmle}.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#' @param bounds Minimum and maximum values to be placed on the \code{ftype}.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom plyr join
#'
#' @return A list of \code{data.frame} objects as described above.
#'

makeDataList <- function(dat, J, ntrt, uniqtrt, t0, bounds = NULL, ...) {
  n <- nrow(dat)
  dataList <- vector(mode = "list", length = ntrt + 1)
  rankftime <- match(dat$ftime, sort(unique(dat$ftime)))

  # first element used for estimation
  dataList[[1]] <- dat[rep(1:nrow(dat), rankftime), ]
  for (j in J) {
    dataList[[1]][[paste0("N", j)]] <- 0
    dataList[[1]][[paste0("N", j)]][cumsum(rankftime)] <- as.numeric(dat$ftype == j)
  }
  dataList[[1]]$C <- 0
  dataList[[1]]$C[cumsum(rankftime)] <- as.numeric(dat$ftype == 0)

  n.row.ii <- nrow(dataList[[1]])
  uniqftime <- unique(dat$ftime)
  orduniqftime <- uniqftime[order(uniqftime)]
  row.names(dataList[[1]])[row.names(dataList[[1]]) %in%
    paste(row.names(dat))] <- paste0(row.names(dat), ".0")
  dataList[[1]]$t <- orduniqftime[as.numeric(paste(unlist(strsplit(
    row.names(dataList[[1]]), ".",
    fixed = TRUE
  ))[seq(2, n.row.ii * 2, 2)])) + 1]

  if (!is.null(bounds)) {
    boundFormat <- data.frame(t = bounds$t)
    for (j in J) {
      if (paste("l", j, sep = "") %in% colnames(bounds)) {
        boundFormat[[paste0("l", j)]] <- bounds[, paste0("l", j)]
      } else {
        boundFormat[[paste0("l", j)]] <- 0
      }
      if (paste("u", j, sep = "") %in% names(bounds)) {
        boundFormat[[paste0("u", j)]] <- bounds[, paste0("u", j)]
      } else {
        boundFormat[[paste0("u", j)]] <- 1
      }
    }
    suppressMessages(
      dataList[[1]] <- plyr::join(
        x = dataList[[1]], y = boundFormat,
        type = "left"
      )
    )
    # if any bounds are missing, add in 0 and 1
    for (j in J) {
      tmp <- is.na(dataList[[1]][, paste0("l", j)])
      dataList[[1]][tmp, paste0("l", j)] <- 0
      tmp <- is.na(dataList[[1]][, paste0("u", j)])
      dataList[[1]][tmp, paste0("u", j)] <- 1
    }
  } else {
    for (j in J) {
      dataList[[1]][[paste0("l", j)]] <- 0
      dataList[[1]][[paste0("u", j)]] <- 1
    }
  }

  # subsequent elements used for prediction
  for (i in seq_len(ntrt)) {
    dataList[[i + 1]] <- dat[sort(rep(1:nrow(dat), t0)), ]
    dataList[[i + 1]]$t <- rep(1:t0, n)
    for (j in J) {
      typejEvents <- dat$id[which(dat$ftype == j)]
      dataList[[i + 1]][[paste0("N", j)]] <- 0
      dataList[[i + 1]][[paste0("N", j)]][dataList[[i + 1]]$id %in% typejEvents &
        dataList[[i + 1]]$t >= dataList[[i + 1]]$ftime] <- 1
    }
    censEvents <- dat$id[which(dat$ftype == 0)]
    dataList[[i + 1]]$C <- 0
    dataList[[i + 1]]$C[dataList[[i + 1]]$id %in% censEvents &
      dataList[[i + 1]]$t >= dataList[[i + 1]]$ftime] <- 1
    dataList[[i + 1]]$trt <- uniqtrt[i]
    # so all indicators pass
    dataList[[i + 1]]$ftime <- t0

    if (!is.null(bounds)) {
      suppressMessages(
        dataList[[i + 1]] <- plyr::join(
          x = dataList[[i + 1]], y = boundFormat,
          type = "left"
        )
      )
      # if any bounds are missing, add in 0 and 1
      for (j in J) {
        tmp <- is.na(dataList[[i + 1]][, paste0("l", j)])
        dataList[[i + 1]][tmp, paste0("l", j)] <- 0
        tmp <- is.na(dataList[[i + 1]][, paste0("u", j)])
        dataList[[i + 1]][tmp, paste0("u", j)] <- 1
      }
    } else {
      for (j in J) {
        dataList[[i + 1]][[paste0("l", j)]] <- .Machine$double.eps
        dataList[[i + 1]][[paste0("u", j)]] <- 1 - .Machine$double.eps
      }
    }
  }
  names(dataList) <- c("obs", uniqtrt)
  return(dataList)
}
