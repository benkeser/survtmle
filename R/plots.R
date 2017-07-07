utils::globalVariables(c("value", "group"))

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
#' @examples
#' library(survtmle)
#' set.seed(341796)
#' n <- 100
#' t_0 <- 10
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, 0.5)
#' T <- rgeom(n,plogis(-4 + W$W1 * W$W2 - A)) + 1
#' C <- rgeom(n, plogis(-6 + W$W1)) + 1
#' ftime <- pmin(T, C)
#' ftype <- as.numeric(ftime == T)
#' suppressWarnings(
#'   fit <- survtmle(ftime = ftime, ftype = ftype,
#'                   adjustVars = W, glm.ftime = "I(W1*W2) + trt + t",
#'                   trt = A, glm.ctime = "W1 + t", method = "hazard",
#'                   verbose = TRUE,  t0 = t_0, maxIter = 2)
#' )
#' tpfit <- timepoints(fit, times = seq_len(t_0))
#' plot(tpfit, type = "smooth")
#'
plot.tp.survtmle <- function(x, ..., type = c("smooth", "raw")) {

  # check that input for type is appropriate
  type <- match.arg(type)

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
                       ggplot2::aes(x = t, y = value, colour = group)
                      )
  p <- p + ggplot2::geom_line()
  p <- p + ggplot2::xlab("Time")
  p <- p + ggplot2::ylab("Cumulative Incidence Estimate")
  p <- p + ggplot2::ggtitle(paste("Cumulative Incidence Amongst Groups",
                                  ifelse(type == "smooth",
                                         "\n (smoothed by isotonic regression)",
                                         "\n (raw estimates)")
                                  )
                           )
  p <- p + ggplot2::theme_bw()
  return(p)
}
