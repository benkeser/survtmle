#' Convert Long Form Data to List of Wide Form Data
#'
#' @description The function takes a \code{data.frame} and \code{list}
#'  consisting of short and long format right-censored failure times. The
#'  function reshapes the long format into the wide format needed for calls to
#'  \code{\link{mean_tmle}}. The list returned by the function will have number
#'  of entries equal to \code{length(trtOfInterest) + 1}. The first will
#'  contain the observed \code{trt} columns and will set \code{C.t} (the
#'  censoring counting process) equal to the observed value of censoring.
#'  Subsequent entries will set \code{trt} equal to each level of
#'  \code{trtOfInterest} and set \code{C.t} to zero for everyone.
#'
#' @param dat The short form \code{data.frame}
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param uniqtrt The values of \code{trtOfInterest} passed to
#'  \code{\link{mean_tmle}}.
#' @param adjustVars A data.frame of adjustment variables that will be used in
#'  estimating the conditional treatment, censoring, and failure (hazard or
#'  conditional mean) probabilities.
#' @param dataList A list of long format \code{data.frame} objects. See the
#'  documentation of \code{\link{makeDataList}} for more details on formatting.
#' @param t0 The timepoint at which \code{\link{survtmle}} was called to
#'  evaluate.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats reshape
#'
#' @return A list of \code{data.frame} objects as described above.
makeWideDataList <- function(dat,
                             allJ,
                             uniqtrt,
                             adjustVars,
                             dataList,
                             t0, ...) {
  wideDataList <- vector(mode = "list", length = length(dataList))
  wideDataList[[1]] <- data.frame(
    wts = dat$wts, trt = dat$trt, dat[, names(adjustVars)],
    stats::reshape(
      dataList[[2]][, !(names(dataList[[2]]) %in%
        c(
          "trt", names(adjustVars), "wts",
          "ftime", "ftype"
        ))],
      direction = "wide",
      timevar = "t", idvar = "id"
    )
  )

  # set Nj0=0 for all j -- makes things easier to run in a loop later
  eval(parse(text = paste0(
    paste0(
      "wideDataList[[1]]$N", allJ, ".0",
      collapse = "<-"
    ),
    "<- wideDataList[[1]]$C.0 <- 0"
  )))

  wideDataList[2:length(dataList)] <- lapply(
    dataList[2:length(dataList)],
    function(x) {
      out <- data.frame(
        wts = dat$wts, dat[, names(adjustVars)],
        stats::reshape(
          x[, !(names(x) %in%
            c(
              "trt", names(adjustVars), "wts", "ftime",
              "ftype"
            ))],
          direction = "wide", timevar = "t", idvar = "id"
        ),
        row.names = NULL
      )
      out[, paste0("C.", seq_len(t0))] <- 0

      # set Nj0=0 for all j -- makes things easier to run in a loop later
      eval(parse(text = paste0(
        paste0("out$N", allJ, ".0", collapse = "<-"),
        "<- out$C.0 <- 0"
      )))
      return(out)
    })
  names(wideDataList) <- c("obs", uniqtrt)

  # append counterfactual values of treatment
  for (z in uniqtrt) wideDataList[[paste0(z)]]$trt <- z

  wideDataList <- lapply(wideDataList, function(x) {
    # make clever covariates
    for (z in uniqtrt) {
      for (t in seq_len(t0)) {
        x[[paste0("H", z, ".", t)]] <-
          (x$trt == z & x[[paste0("C.", t - 1)]] == 0) /
            (x[[paste0("G_dC.", t)]] * x[[paste0("g_", z, ".", t)]])
      }
      x[[paste0("H", z, ".0")]] <- (x$trt == z) / x[[paste0("g_", z, ".", t)]]
    }
    return(x)
  })
  return(wideDataList)
}
