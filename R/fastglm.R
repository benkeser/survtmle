#' Faster Yet Robust Generalized Linear Models
#'
#' words, words, words, ...
#'
#' @param reg_form ...
#' @param data ...
#' @param family ...
#' @param start ...
#' @param flavor Character: do you want to live dangerously?
#'
#' @importFrom speedglm speedglm
#' @importFrom stats glm as.formula
#'
#' @return ...
#'

fast_glm <- function(reg_form, data, family, flavor = "fast", start = NULL) {
  if(flavor == "fast") {
    out <- tryCatch(
      {
        speedglm::speedglm(stats::as.formula(reg_form),
                           data = data,
                           family = family,
                           start = start)
      },
      error = function(cond) {
        message("'speedglm' encountered an error, reverting to using 'glm'")
        message("Here's the original error message (from 'speedglm'):")
        message(cond)
        # Choose a return value in case of error
        mod <- stats::glm(stats::as.formula(reg_form),
                          data = data,
                          family = family,
                          start = start)
        return(mod)
      }
      )
    return(out)
  } else {
    out <- stats::glm(stats::as.formula(reg_form),
                      data = data,
                      family = family,
                      start = start)
  }
  return(out)
}
