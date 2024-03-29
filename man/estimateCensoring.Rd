% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/censoring_estimate.R
\name{estimateCensoring}
\alias{estimateCensoring}
\title{Estimate Censoring Mechanisms}
\usage{
estimateCensoring(
  dataList,
  adjustVars,
  t0,
  SL.ctime = NULL,
  glm.ctime = NULL,
  glm.family,
  cvControl,
  returnModels = FALSE,
  verbose = TRUE,
  gtol = 0.001,
  stratify = FALSE,
  trtOfInterest,
  ...
)
}
\arguments{
\item{dataList}{A list of \code{data.frame} objects as described in the
documentation of \code{\link{makeDataList}}.}

\item{adjustVars}{Object of class \code{data.frame} that contains the
variables to adjust for in the regression.}

\item{t0}{The timepoint at which \code{survtmle} was called to evaluate.
Needed only because the naming convention for the regression if
\code{t == t0} is different than if \code{t != t0}.}

\item{SL.ctime}{A character vector or list specification to be passed to the
\code{SL.library} argument of \code{\link[SuperLearner]{SuperLearner}} for
the outcome regression (either cause-specific hazards or conditional mean).
See the documentation of \code{\link[SuperLearner]{SuperLearner}} for more
information on how to specify valid \code{SuperLearner} libraries. It is
expected that the wrappers used in the library will play nicely with the
input variables, which will be called \code{"trt"} and
\code{names(adjustVars)}.}

\item{glm.ctime}{A character specification of the right-hand side of the
equation passed to the \code{\link[stats]{formula}} option of a call to
\code{\link[stats]{glm}} for the outcome regression (either cause-specific
hazards or conditional mean). Ignored if \code{SL.ctime != NULL}. Use
\code{"trt"} to specify the treatment in this formula (see examples). The
The formula can additionally include any variables found in
\code{names(adjustVars)}.}

\item{glm.family}{The type of regression to be performed if fitting GLMs in
the estimation and fluctuation procedures. The default is "binomial" for
logistic regression. Only change this from the default if there are
justifications that are well understood. This is inherited from the calling
function (either \code{\link{mean_tmle}} or \code{\link{hazard_tmle}}).}

\item{cvControl}{A \code{list} providing control options to be fed directly
into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
the contents of \code{SuperLearner.CV.control} exactly. For details,
consult the documentation of the \pkg{SuperLearner} package. This is passed
in from \code{\link{mean_tmle}} or \code{\link{hazard_tmle}} via
\code{\link{survtmle}}.}

\item{returnModels}{A \code{logical} indicating whether to return the
\code{glm} or \code{SuperLearner} objects used to estimate the nuisance
parameters. Must be set to \code{TRUE} to make downstream calls to
\code{\link{timepoints}} for obtaining estimates at times other than
\code{t0}. See documentation of \code{\link{timepoints}} for more
information.}

\item{verbose}{A \code{logical} indicating whether the function should print
messages to indicate progress.}

\item{gtol}{The truncation level of predicted censoring survival to handle
positivity violations.}

\item{stratify}{If \code{TRUE}, then the censoring hazard model is estimated using only
the observations with \code{trt == trtOfInterest}. Only works if 
\code{length(trtOfInterest) == 1}. If \code{stratify = TRUE} then \code{glm.ftime}
cannot include \code{trt} in the model formula and any learners in \code{SL.ftime}
should not assume a variable named \code{trt} will be included in the candidate 
super learner estimators.}

\item{trtOfInterest}{An input specifying which levels of \code{trt} are of
interest. The default value computes estimates for all values in
\code{unique(trt)}. Can alternatively be set to a vector of values found in
\code{trt}. Ignored unless \code{stratify == TRUE}.}

\item{...}{Other arguments. Not currently used.}
}
\value{
The function returns a list that is exactly the same as the input
 \code{dataList}, but with a column named \code{G_dC} added to it, which is
 the estimated conditional survival distribution for the censoring variable
 evaluated at the each of the rows of each \code{data.frame} in
 \code{dataList}.
}
\description{
Computes an estimate of the hazard for censoring using either
 \code{\link[stats]{glm}} or \code{\link[SuperLearner]{SuperLearner}} based
 on log-likelihood loss. The function then computes the censoring survival
 distribution based on these estimates. The structure of the function is
 specific to how it is called within \code{\link{survtmle}}. In particular,
 \code{dataList} must have a very specific structure for this function to
 run properly. The list should consist of \code{data.frame} objects. The
 first will have the number of rows for each observation equal to the
 \code{ftime} corresponding to that observation. Subsequent entries will
 have \code{t0} rows for each observation and will set \code{trt} column
 equal to each value of \code{trtOfInterest} in turn. One of these columns
 must be named \code{C} that is a counting process for the right-censoring
 variable. The function will fit a regression with \code{C} as the outcome
 and functions of \code{trt} and \code{names(adjustVars)} as specified by
 \code{glm.ctime} or \code{SL.ctime} as predictors.
}
