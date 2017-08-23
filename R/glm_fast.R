#' Wrapper for faster Generalized Linear Models
#'
#' A convenience utility to fit regression models more quickly in the main
#' internal functions for estimation, which usually require logistic regression.
#' Use of \code{speedglm} appears to provide roughly an order of magnitude
#' improvement in speed when compared to \code{glm} in custom benchmarks.
#'
#' @param reg_form Object of class \code{formula} indicating the regression to
#'        be fit.
#' @param data Object of class \code{data.frame} containing the data.
#' @param family Object of class \code{family} from package \code{stats}
#'        indicating the error distribution. Appropriate options are limited to
#'        \code{gaussian} and \code{binomial}.
#' @param ... Additional arguments passed to \code{glm} or \code{speedglm}.
#'
#' @importFrom speedglm speedglm
#' @importFrom stats glm gaussian binomial
#' @importFrom stringr str_split
#'
#' @return Object of class \code{glm} or \code{speedglm}.
#'

fast_glm <- function(reg_form, data, family, ...) {
  # a quick type check for safety
  stopifnot(
    class(reg_form) == "formula" &
    class(data) %in% c("matrix", "data.frame") &
    class(family) == "family"
  )

  # catch the calling function
  calling_fun <- as.character(stringr::str_split(deparse(sys.call(-1)),
                                                 "\\(")[[1]][1])

  # fit speedglm or glm as appropriate
  out <- tryCatch(
    {
      # Obviously, a sparse design matrix is not used when fitting an intercept
      # model. In such cases, 'sparse=TRUE' is an inappropriate choice, though
      # this is only expected when estimating the treatment mechanism.
      speedglm::speedglm(formula = reg_form,
                         data = data,
                         family = family,
                         method = "Cholesky",
                         sparse = ifelse(calling_fun == "estimateTreatment",
                                         FALSE, TRUE),
                         trace = FALSE,
                         ...)
    },
    error = function(cond) {
      message(paste0("'speedglm' ran into an error in ", calling_fun,
                    ".", "'glm' will be used instead."))
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
