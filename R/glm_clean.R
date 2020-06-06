#' Clean up outputs from GLM
#'
#' @description Removes superfluous output from the call to
#'  \code{\link[stats]{glm}} that is not needed to perform later predictions.
#'  It is applied as a space saving technique.
#'
#' @param cm An object of class \code{\link[stats]{glm}} or class
#'  \code{\link[speedglm]{speedglm}}.
#'
#' @return An object of class \code{\link[stats]{glm}} or
#'  \code{\link[speedglm]{speedglm}}, but with unnecessary output removed.
cleanglm <- function(cm) {
  cm$y <- NULL
  cm$model <- NULL
  cm$residuals <- NULL
  cm$fitted.values <- NULL
  cm$effects <- NULL
  cm$qr$qr <- NULL
  cm$linear.predictors <- NULL
  cm$weights <- NULL
  cm$prior.weights <- NULL
  cm$data <- NULL
  cm$family$variance <- NULL
  cm$family$dev.resids <- NULL
  cm$family$aic <- NULL
  cm$family$validmu <- NULL
  cm$family$simulate <- NULL
  attr(cm$terms, ".Environment") <- NULL
  attr(cm$formula, ".Environment") <- NULL
  return(cm)
}
