#' Estimation for the Method of Cause-Specific Hazards
#'
#' @description This function computes an estimate of the cause-specific hazard
#'  functions over all times using either \code{\link[stats]{glm}} or
#'  \code{\link[SuperLearner]{SuperLearner}}. The structure of the function is
#'  specific to how it is called within \code{\link{hazard_tmle}}. In
#'  particular, \code{dataList} must have a very specific structure for this
#'  function to run properly. The list should consist of \code{data.frame}
#'  objects. The first will have the number of rows for each observation equal
#'  to the \code{ftime} corresponding to that observation. Subsequent entries
#'  will have \code{t0} rows for each observation and will set the \code{trt}
#'  column equal to each value of \code{trtOfInterest} in turn. The function
#'  uses the first entry in \code{dataList} to iteratively fit hazard
#'  regression models for each cause of failure. Thus, this \code{data.frame}
#'  needs to have a column called \code{Nj} for each value of j in \code{J}.
#'  The first fit estimates the hazard of \code{min(J)}, while subsequent fits
#'  estimate the pseudo-hazard of all other values of j, where pseudo-hazard is
#'  used to mean the probability of a failure due to type j at a particular
#'  timepoint given no failure of any type at any previous timepoint AND no
#'  failure due to type \code{k < j} at a particular timepoint. The hazard
#'  estimates of causes j' can then be used to map this pseudo-hazard back into
#'  the hazard at a particular time. This is nothing more than the re-framing
#'  of a conditional multinomial probability into a series of conditional
#'  binomial probabilities. This structure ensures that no strata have
#'  estimated hazards that sum to more than one over all possible causes of
#'  failure at a particular timepoint.
#'
#' @param dataList A \code{list} of \code{data.frame} objects.
#' @param J Numeric vector indicating the labels of all causes of failure.
#' @param adjustVars Object of class \code{data.frame} containing the variables
#'  to adjust for in the regression.
#' @param SL.ftime A character vector or list specification to be passed to the
#'  \code{SL.library} argument of \code{\link[SuperLearner]{SuperLearner}} for
#'  the outcome regression (either cause-specific hazards or conditional mean).
#'  See the documentation of \code{\link[SuperLearner]{SuperLearner}} for more
#'  information on how to specify valid \code{SuperLearner} libraries. It is
#'  expected that the wrappers used in the library will play nicely with the
#'  input variables, which will be called \code{"trt"} and
#'  \code{names(adjustVars)}.
#' @param glm.ftime A character specification of the right-hand side of the
#'  equation passed to the \code{\link[stats]{formula}} option of a call to
#'  \code{\link[stats]{glm}} for the outcome regression (either using
#'  cause-specific hazards or conditional mean). Ignored if \code{SL.ftime !=
#'  NULL}. Use \code{"trt"} to specify the treatment in this formula (see
#'  examples). The formula can additionally include any variables found in
#'  \code{names(adjustVars)}.
#' @param glm.family The type of regression to be performed if fitting GLMs in
#'  the estimation and fluctuation procedures. The default is "binomial" for
#'  logistic regression. Only change this from the default if there are
#'  justifications that are well understood. This is inherited from the calling
#'  function (either \code{\link{mean_tmle}} or \code{\link{hazard_tmle}}).
#' @param cvControl A \code{list} providing control options to be fed directly
#'  into calls to \code{\link[SuperLearner]{SuperLearner}}. This should match
#'  the contents of \code{SuperLearner.CV.control} exactly. For details,
#'  consult the documentation of the \pkg{SuperLearner} package. This is passed
#'  in from \code{\link{mean_tmle}} or \code{\link{hazard_tmle}} via
#'  \code{\link{survtmle}}.
#' @param returnModels A \code{logical} indicating whether to return the
#'  \code{glm} or \code{SuperLearner} objects used to estimate the nuisance
#'  parameters. Must be set to \code{TRUE} to make downstream calls to
#'  \code{\link{timepoints}} for obtaining estimates at times other than
#'  \code{t0}. See documentation of \code{\link{timepoints}} for more
#'  information.
#' @param bounds A list of bounds... TODO: Add more description here.
#' @param verbose A \code{logical} indicating whether the function should print
#'  messages to indicate progress.
#' @param trtOfInterest An input specifying which levels of \code{trt} are of
#'  interest. The default value computes estimates for all values in
#'  \code{unique(trt)}. Can alternatively be set to a vector of values found in
#'  \code{trt}. Ignored unless \code{stratify == TRUE}.
#' @param stratify If \code{TRUE}, then the hazard model is estimated using only
#'  the observations with \code{trt == trtOfInterest}. Only works if 
#'  \code{length(trtOfInterest) == 1}. If \code{stratify = TRUE} then \code{glm.ftime}
#'  cannot include \code{trt} in the model formula and any learners in \code{SL.ftime}
#'  should not assume a variable named \code{trt} will be included in the candidate 
#'  super learner estimators. 
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats as.formula predict model.matrix optim glm binomial
#' @importFrom speedglm speedglm
#' @importFrom SuperLearner SuperLearner
#'
#' @return The function returns a list that is exactly the same as the input
#'  \code{dataList}, but with additional columns corresponding to the hazard
#'  pseudo-hazard, and the total hazard summed over all causes \code{k < j}.
#'
#' @export
estimateHazards <- function(dataList,
                            J,
                            adjustVars,
                            SL.ftime = NULL,
                            glm.ftime = NULL,
                            glm.family,
                            cvControl,
                            returnModels,
                            bounds,
                            verbose,
                            trtOfInterest,
                            stratify,
                            ...) {
  
  if(length(trtOfInterest) > 1 & stratify){
    stratify <- FALSE
    warning("stratify option only supported if there is only a single trtOfInterest.")
  }
  if(stratify){
    stratify_include <- dataList[[1]]$trt == trtOfInterest[1]
  }else{
    stratify_include <- rep(TRUE, length(dataList[[1]][,1]))
  }
  ftimeMod <- vector(mode = "list", length = length(J))
  names(ftimeMod) <- paste0("J", J)


  ## determine whether to use linear or logistic regression in GLM fit
  if (!is.null(glm.family)) {
    glm_family <- parse(text = paste0("stats::", glm.family, "()"))
  }

  if (is.null(SL.ftime)) {
    if (is.null(bounds)) {
      for (j in J) {
        # formula
        Qj_form <- sprintf("%s ~ %s", paste("N", j, sep = ""), glm.ftime)

        # sum all events less than current j to see who to include in regression
        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for (i in J[J < j]) {
          NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
        }

        # fit GLM
        if (!("glm" %in% class(glm.ftime[[1]])) &
          !("speedglm" %in% class(glm.ftime[[1]]))) {
          Qj_mod <- fast_glm(
            reg_form = stats::as.formula(Qj_form),
            data = dataList[[1]][NlessthanJ == 0 & stratify_include, , drop = FALSE],
            family = eval(glm_family)
          )
          if (unique(class(Qj_mod) %in% c("glm", "lm"))) {
            Qj_mod <- cleanglm(Qj_mod)
          }
        } else {
          Qj_mod <- glm.ftime[[paste0("J", j)]]
        }
        ftimeMod[[paste0("J", j)]] <- NULL
        if (returnModels) ftimeMod[[paste0("J", j)]] <- Qj_mod

        # get predictions back
        dataList <- lapply(dataList, function(x, j) {
          suppressWarnings(
            x[[paste0("Q", j, "PseudoHaz")]] <- predict(
              Qj_mod,
              newdata = x,
              type = "response"
            )
          )
          if (j != min(J)) {
            x[[paste0("hazLessThan", j)]] <- rowSums(cbind(
              rep(0, nrow(x)),
              x[, paste0("Q", J[J < j], "Haz")]
            ))
            x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]] *
              (1 - x[[paste0("hazLessThan", j)]])
          } else {
            x[[paste0("hazLessThan", j)]] <- 0
            x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]]
          }
          x
        }, j = j)
      }
    } else {
      for (j in J) {
        Qj_form <- sprintf("%s ~ %s", paste("N", j, sep = ""), glm.ftime)
        X <- stats::model.matrix(
          stats::as.formula(Qj_form),
          data = dataList[[1]]
        )

        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for (i in J[J < j]) {
          NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
        }

        dataList <- lapply(dataList, function(x, j) {
          if (j != min(J)) {
            x[[paste0("hazLessThan", j)]] <- rowSums(cbind(
              rep(0, nrow(x)),
              x[, paste0("Q", J[J < j], "Haz")]
            ))
          } else {
            x[[paste0("hazLessThan", j)]] <- 0
          }
          x
        }, j = j)
        Ytilde <- (dataList[[1]][[paste0("N", j)]] -
          dataList[[1]][[paste0("l", j)]]) /
          (pmin(dataList[[1]][[paste0("u", j)]], 1 -
            dataList[[1]][[paste0("hazLessThan", j)]]) -
            dataList[[1]][[paste0("l", j)]])

        if (class("glm.ftime") != "list") {
          include <- (NlessthanJ == 0 & stratify_include)
          Qj_mod <- stats::optim(
            par = rep(0, ncol(X)), fn = LogLikelihood,
            Y = Ytilde[include], X = X[include, drop = FALSE], 
            method = "BFGS", gr = grad,
            control = list(reltol = 1e-7, maxit = 50000)
          )
        } else {
          Qj_mod <- glm.ftime[[paste0("J", j)]]
        }
        if (Qj_mod$convergence != 0) {
          stop("convergence failure")
        } else {
          beta <- Qj_mod$par
          eval(parse(text = paste0("ftimeMod$J", j, " <- Qj_mod")))
          dataList <- lapply(dataList, function(x, j) {
            newX <- stats::model.matrix(stats::as.formula(Qj_form), data = x)
            x[[paste0("Q", j, "PseudoHaz")]] <- plogis(newX %*% beta)
            x[[paste0("Q", j, "Haz")]] <-
              (pmin(x[[paste0("u", j)]], 1 - x[[paste0("hazLessThan", j)]]) -
                x[[paste0("l", j)]]) * x[[paste0("Q", j, "PseudoHaz")]] +
              x[[paste0("l", j)]]
            x
          }, j = j)
        }
      }
    }
  } else if (is.null(glm.ftime)) {
    for (j in J) {
      # add all events less than current j to see who to include in regression
      NlessthanJ <- rep(0, nrow(dataList[[1]]))
      for (i in J[J < j]) {
        NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
      }

      if (class(SL.ftime[[1]]) != "SuperLearner") {
        Qj_mod <- SuperLearner::SuperLearner(
          Y = dataList[[1]][[paste0("N", j)]][NlessthanJ == 0 & stratify_include],
          X = dataList[[1]][
            NlessthanJ == 0 & stratify_include,
            c("t", "trt", names(adjustVars)),
            drop = FALSE
          ],
          id = dataList[[1]]$id[NlessthanJ == 0 & stratify_include],
          family = stats::binomial(),
          SL.library = SL.ftime,
          cvControl = cvControl,
          verbose = verbose
        )
      } else {
        Qj_mod <- SL.ftime[[paste0("J", j)]]
      }
      ftimeMod[[paste0("J", j)]] <- Qj_mod

      # get predictions back
      dataList <- lapply(dataList, function(x, j) {
        suppressWarnings(
          x[[paste0("Q", j, "PseudoHaz")]] <- predict(
            Qj_mod,
            onlySL = TRUE,
            newdata = x[, c("t", "trt", names(adjustVars))]
          )[[1]]
        )
        if (j != min(J)) {
          x[[paste0("hazLessThan", j)]] <- rowSums(cbind(
            rep(0, nrow(x)),
            x[, paste0("Q", J[J < j], "Haz")]
          ))
          x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]] *
            (1 - x[[paste0("hazLessThan", j)]])
        } else {
          x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]]
          x[[paste0("hazLessThan", j)]] <- 0
        }
        x
      }, j = j)
    }
  }
  out <- list()
  out$dataList <- dataList
  out$ftimeMod <- ftimeMod
  return(out)
}
