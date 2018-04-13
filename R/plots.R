utils::globalVariables(c(".", "value", "group"))

#' Plot Results of Cumulative Incidence Estimates
#'
#' Step function plots for both raw and smoothed (monotonic) estimates, the
#' latter by isotonic regression of the raw estimates, of cumulative incidence.
#'
#' @param x object of class \code{tp.survtmle} as produced by a sequence of
#'  appropriate calls to \code{survtmle} and \code{timepoints}
#' @param type \code{character} describing whether to provide a plot of raw
#'  ("raw") or monotonic ("iso") estimates in the resultant step function plot,
#'  with the latter being computed by a call to \code{stats::isoreg}
#' @param pal A \code{ggplot2} palette object from the \code{ggsci} package. The
#'  default of \code{scale_color_lancet} is generally appropriate for medical
#'  and epidemiologic applications, though there are situations in which one
#'  might opt to change this. Note that this can also be overridden in the
#'  resultant plot object using standard \code{ggplot2} semantics.
#' @param ... additional arguments passed \code{plot} as necessary
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_step xlab ylab ggtitle
#' @importFrom ggsci scale_color_lancet
#' @importFrom stringr str_length str_sub
#' @importFrom tidyr gather
#' @importFrom dplyr "%>%"
#' @importFrom stats isoreg
#'
#' @return object of class \code{ggplot} containing a step function plot of the
#'  raw or smoothened point estimates of cumulative incidence across a series of
#'  timepoints of interest.
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
#' plot(tpfit)
#
plot.tp.survtmle <- function(x,
                             ...,
                             type = c("iso", "raw"),
                             pal = ggsci::scale_color_lancet()) {

  # check that input for type is appropriate
  type <- match.arg(type)

  # extract point estimates from tp.survtmle input object and make tidy
  est <- lapply(x, function(x) {
    x$est
  })
  est <- Reduce(cbind, est)

  # extract time points of interest by actual values rather than order
  times <- objects(x)
  times_labels <- stringr::str_sub(times, 2, stringr::str_length(times))
  times_labels <- as.numeric(unclass(times_labels))
  times_labels <- times_labels[order(times_labels)] # re-order

  if (type == "raw") {
    raw_est_in <- as.data.frame(cbind(t(est), times_labels))
    colnames(raw_est_in) <- c(gsub(
      " ", "/",
      colnames(raw_est_in)[seq_len(ncol(raw_est_in)
      - 1)]
    ), "t")
    raw_est_in <- tidyr::gather(data = raw_est_in, key = t)
    raw_est_in <- as.data.frame(cbind(
      rep(
        times_labels,
        length(unique(raw_est_in$t))
      ),
      raw_est_in
    ))
    colnames(raw_est_in) <- c("t", "group", "value")
    raw_est_in[, "group"] <- as.factor(raw_est_in[, "group"])
    plot_in <- raw_est_in
  } else if (type == "iso") {
    iso <- apply(est, 1, function(y) {
      tmp <- stats::isoreg(y = y, x = seq_len(length(x)))
      tmp$yf
    })
    iso_est_in <- as.data.frame(cbind(as.matrix(iso), times_labels))
    colnames(iso_est_in) <- c(gsub(
      " ", "/",
      colnames(iso_est_in)[seq_len(ncol(iso_est_in)
      - 1)]
    ), "t")
    iso_est_in <- tidyr::gather(data = iso_est_in, key = t)
    iso_est_in <- as.data.frame(cbind(
      rep(
        times_labels,
        length(unique(iso_est_in$t))
      ),
      iso_est_in
    ))
    colnames(iso_est_in) <- c("t", "group", "value")
    iso_est_in[, "group"] <- as.factor(iso_est_in[, "group"])
    plot_in <- iso_est_in
  }
  # generate output plot
  p <- ggplot2::ggplot(
    data = plot_in,
    ggplot2::aes(x = t, y = value, colour = group)
  )
  if (length(unique(plot_in$t)) > 1) {
    p <- p + ggplot2::geom_step() + ggplot2::geom_point()
  } else {
    p <- p + ggplot2::geom_point()
  }
  p <- p + ggplot2::xlab("Time") +
    ggplot2::ylab("Cumulative Incdicence Estimate") +
    ggplot2::ggtitle(paste(
      "Cumulative Incidence Amongst Groups",
      ifelse(type == "iso",
        "\n (smoothed by isotonic regression)",
        "\n (raw estimates)"
      )
    )) +
    ggplot2::theme_bw() +
    pal
  # ...and print
  return(p)
}
