#' Faster Yet Robust Generalized Linear Models
#'
#' words, words, words, ...
#'
#' @param reg_form ...
#' @param data ...
#' @param err_fam ...
#' @param start ...
#'
#' @importFrom speedglm speedglm
#' @importFrom stats glm as.formula
#'
#' @return ...
#'

fast_glm <- function(reg_form, data, err_fam, start = NULL) {
  out <- tryCatch(
    {
      speedglm::speedglm(stats::as.formula(reg_form),
                         data = data,
                         family = err_fam,
                         start = start)
    },
    error = function(cond) {
      message("'speedglm' encountered an error, reverting to using 'glm'")
      message("Here's the original error message (from 'speedglm'):")
      message(cond)
      # Choose a return value in case of error
      mod <- stats::glm(stats::as.formula(reg_form),
                        data = data,
                        family = err_fam,
                        start = start)
      return(mod)
    },
    finally = {
    }
  )
  return(out)
}
