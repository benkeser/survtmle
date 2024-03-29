% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftime_missing_estimate.R
\name{estimateFtimeMissing}
\alias{estimateFtimeMissing}
\title{Estimate Missing Failure Time mechanism}
\usage{
estimateFtimeMissing(
  dat,
  adjustVars,
  glm.ftimeMissing = NULL,
  SL.ftimeMissing = NULL,
  cvControl,
  returnModels = FALSE,
  verbose = FALSE,
  gtol = 0.001,
  trtOfInterest,
  ...
)
}
\arguments{
\item{dat}{An object of class \code{data.frame}. Must have named column
\code{trt}.}

\item{adjustVars}{An object of class \code{data.frame} that will be used
either as the \code{data} argument in a call to \code{\link[stats]{glm}} or
as \code{X} in a call to \code{\link[SuperLearner]{SuperLearner}}.}

\item{glm.ftimeMissing}{A character specification of the right-hand side of the
equation passed to the \code{\link[stats]{formula}} option of a call to
\code{\link[stats]{glm}} for the missingness failure times model. Ignored if
\code{SL.ftimeMissing} is not equal to \code{NULL}. Use \code{"trt"} to specify
the treatment in this formula (see examples). The formula can additionally
include any variables found in \code{names(adjustVars)}.}

\item{SL.ftimeMissing}{A character vector or list specification to be passed to the
\code{SL.library} in the call to \code{\link[SuperLearner]{SuperLearner}}
for the missingness model for failure time.
It is expected that the wrappers used in the library will play nicely with
the input variables, which will be called \code{"trt"} and 
\code{names(adjustVars)}.}

\item{cvControl}{A \code{list} providing control options to be fed directly
into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
the contents of \code{SuperLearner.CV.control} exactly. For details,
consult the documentation of the \pkg{SuperLearner} package. This is passed
in from \code{\link{mean_tmle}} or \code{\link{hazard_tmle}} via
\code{\link{survtmle}}.}

\item{returnModels}{A \code{logical} indicating whether fitted model objects
should be returned.}

\item{verbose}{A \code{logical} passed to the \code{verbose} option of the
call to \code{\link[SuperLearner]{SuperLearner}}.}

\item{gtol}{The truncation level of predicted missingness probabilities to handle
positivity violations.}

\item{trtOfInterest}{An input specifying which levels of \code{trt} are of
interest. The default value computes estimates for all of the values in
\code{unique(trt)}. Can alternatively be set to a vector of values found in
\code{trt}.}

\item{...}{Other arguments. Not currently used}
}
\value{
dat The input \code{data.frame} object with two added columns
 corresponding with the conditional probability (given \code{adjustVars}) of
 \code{trt==max(trt)} and \code{trt==min(trt)}.

ftimeMissingMod If \code{returnModels = TRUE}, the fitted \code{glm} or
 \code{SuperLearner} object. Otherwise, \code{NULL}
}
\description{
This function computes the conditional probability of having
 non-missing \code{ftime} for each level of treatment using \code{\link[stats]{glm}}
 or \code{\link[SuperLearner]{SuperLearner}}.
}
