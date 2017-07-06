#' Plot results of cumulative incidence estimates from survtmle
#'
#' Plotting for both raw and smoothened estimates, the latter by isotonic
#' regression of the raw point estimates, of cumulative incidence.
#'
#' @param x object of class \code{tp.survtmle} as produced by a sequence of
#' appropriate calls to \code{survtmle} and \code{timepoints}
#' @param type character describing whether to provide a plot of raw or
#'        smoothened point estimates, with the latter being computed by a call
#'        to \code{stats::isoreg}
#' @param ... additional arguments passed \code{plot} as necessary
#'
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab ggtitle
#' @importFrom tidyr gather
#' @importFrom stats isoreg
#'
#' @return object of class \code{ggplot} containing a lines plot of the raw or
#'         smoothened point estimates of cumulative incidence across a series of
#'         time points of interest.
#'
#' @export
#'
#' @method plot tp.survtmle
#'
plot.tp.survtmle <- function(x, ..., type = c("smooth", "raw")) {

  # check that input for type is appropriate
  type <- match.arg(type)

  # capture call for later use
  call <- match.call(expand.dots = TRUE)

  # extract point estimates from tp.survtmle input object and make tidy
  est <- lapply(x, function(x) {x$est})
  est <- Reduce(cbind, est)

  if(type == "raw") {
    est_in <- as.data.frame(cbind(t(est), seq_len(length(x))))
    colnames(est_in) <- c("0/1", "1/1", "t")
    est_in <- tidyr::gather(est_in, t)
    colnames(est_in) <- c("t", "group", "value")
    plot_in <- est_in
  } else if (type == "smooth") {
    iso <- apply(est, 1, function(y) {
      tmp <- stats::isoreg(y = y, x = seq_len(length(x)))
      tmp$yf
    })
    iso_est <- as.data.frame(cbind(iso, seq_len(length(x))))
    colnames(iso_est) <- c("0/1", "1/1", "t")
    iso_est_in <- tidyr::gather(iso_est, t)
    colnames(iso_est_in) <- c("t", "group", "value")
    plot_in <- iso_est_in
  }
  p <- ggplot2::ggplot(data = plot_in,
                       ggplot2::aes(x = t, y =~ value, colour =~ group)
                      )
  p <- p + ggplot2::geom_line()
  p <- p + ggplot2::xlab("time point")
  p <- p + ggplot2::ylab("cumulative incidence")
  p <- p + ggplot2::ggtitle("cumulative incidence by time")
  p <- p + ggplot2::theme_minimal()
  return(p)
}
