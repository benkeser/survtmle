% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_clean.R
\name{cleanglm}
\alias{cleanglm}
\title{Clean up outputs from GLM}
\usage{
cleanglm(cm)
}
\arguments{
\item{cm}{An object of class \code{\link[stats]{glm}} or class
\code{\link[speedglm]{speedglm}}.}
}
\value{
An object of class \code{\link[stats]{glm}} or
 \code{\link[speedglm]{speedglm}}, but with unnecessary output removed.
}
\description{
Removes superfluous output from the call to
 \code{\link[stats]{glm}} that is not needed to perform later predictions.
 It is applied as a space saving technique.
}
