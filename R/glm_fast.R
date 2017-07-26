#' Wrapper for faster Generalized Linear Models
#'
#' A convenience utility to fit regression models more quickly in the main
#' internal functions for estimation, which usually require logistic regression.
#' Use of \code{speedglm} appears to provide roughly an order of magnitude
#' improvement in speed when compared to \code{glm} in custom benchmarks.
#'
#' @param reg_form Character indicating the regression formula to be computed.
#' @param data Object of class \code{data.frame} containing the data.
#' @param family Object from package \code{stats} indicating error distribution.
#' @param ... Additional arguments passed to \code{glm} or \code{speedglm}.
#'
#' @importFrom speedglm speedglm
#' @importFrom stats glm gaussian binomial
#'
#' @return Object of class \code{glm} or \code{speedglm}.
#'

fast_glm <- function(reg_form, data, family, ...) {
  # a quick type check for safety
  stopifnot(class(reg_form) == "formula")

  # catch the function input in case desired later
  #call <- match.call(expand.dots = TRUE)

  # fit speedglm or glm as appropriate
  out <- tryCatch(
    {
      speedglm::speedglm(formula = reg_form,
                         data = data,
                         family = family,
                         fitted = TRUE,
                         ...)
    },
    error = function(cond) {
      message("'speedglm' ran into an error...reverting to use of 'glm'.")
      # Choose a return value in case of error
      mod <- stats::glm(formula = reg_form,
                        data = data,
                        family = family,
                        ...)
      return(mod)
    }
  )
  return(out)
}
