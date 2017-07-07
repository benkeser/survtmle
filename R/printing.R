#' print.survtmle
#'
#' The print method for an object of class \code{survtmle}
#'
#' @param x An object of class \code{survtmle}
#' @param ... Other options (not currently used)
#'
#' @return Prints the estimates of cumulative incidence and the diagonal
#'         of the estimated covariance matrix.
#'
#' @export
#'
#' @method print survtmle
#'

print.survtmle <- function(x, ...) {
  if(length(x$est) <= 4) {
    print(x[c("est", "var")])
  } else {
    out <- list(est = x$est, var = diag(x$var))
    print(out)
  }
}

#' print.tp.survtmle
#'
#' The print method for a timepoints object of class \code{tp.survtmle}
#'
#' @param x An object of class \code{tp.survtmle}.
#' @param ... Other options (not currently used).
#'
#' @return Prints the estimates of cumulative incidence and the diagonal
#'         of the estimated covariance matrix.
#'
#' @export
#'
#' @method print tp.survtmle
#'

print.tp.survtmle <- function(x, ...) {
  # extract basic information from input `tp.survtmle` object
  len_groups <- as.numeric(unique(lapply(lapply(x, FUN = `[[`, "est"),
                                         FUN = length)))
  names_groups <- unique(lapply(lapply(x, FUN = `[[`, "est"),
                                FUN = rownames))[[1]]

  # structure point estimates
  est_only <- t(matrix(unlist(lapply(x, FUN = `[[`, "est")), ncol = len_groups,
                       byrow = TRUE))
  est_only <- as.data.frame(est_only)
  rownames(est_only) <- names_groups
  colnames(est_only) <- paste0("t", seq_len(ncol(est_only)))

  # structure variance estimates
  vars_only <- lapply(lapply(x, FUN = `[[`, "var"), FUN = diag)
  vars_only <- t(matrix(unlist(vars_only), ncol = len_groups, byrow = TRUE))
  vars_only <- as.data.frame(vars_only)
  rownames(vars_only) <- names_groups
  colnames(vars_only) <- paste0("t", seq_len(ncol(vars_only)))

  # create output object as list to match `print.cuminc` for class from `cmprsk`
  out <- list(est = est_only, var = vars_only)

  print(out)
}
