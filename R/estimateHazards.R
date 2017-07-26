#' Estimation for the Method of Cause-Specific Hazards
#'
#' This function computes an estimate of the cause-specific hazard functions
#' over all times using either \code{glm} or \code{SuperLearner}. The structure
#' of the function is specific to how it is called within \code{hazard_tmle}. In
#' particular, \code{dataList} must have a very specific structure for this
#' function to run properly. The list should consist of \code{data.frame}
#' objects. The first will have the number of rows for each observation equal to
#' the \code{ftime} corresponding to that observation. The subsequent entries
#' will have \code{t0} rows for each observation and will set \code{trt} column
#' equal to each value of \code{trtOfInterest} in turn. The function uses the
#' first entry in \code{dataList} to iteratively fit hazard regression models
#' for each cause of failure. Thus, this \code{data.frame} needs to have a
#' column called \code{Nj} for each value of j in \code{J}. The first fit
#' estimates the hazard of \code{min(J)}, while subsequent fits estimate the
#' pseudo-hazard of all other values of j, where pseudo-hazard is used to mean
#' the probability of a failure due to type j at a particular timepoint given
#' no failure of any type at any previous timepoint AND no failure due to type
#' \code{k < j} at a particular timepoint. The hazard estimates of causes j'
#' can then be used to map this pseudo-hazard back into the hazard at a
#' particular time. This is nothing more than the re-framing of a conditional
#' multinomial probability into a series of conditional binomial probabilities.
#' This structure ensures that no strata have estimated hazards that sum to more
#' than one over all possible causes of failure at a particular timepoint.
#'
#' @param dataList A list of \code{data.frame} objects.
#' @param J Numeric vector indicating the labels of all causes of failure.
#' @param adjustVars Object of class \code{data.frame} that contains the
#'        variables to adjust for in the regression.
#' @param SL.ftime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        outcome regression (either cause-specific hazards or conditional mean).
#'        See \code{?SuperLearner} for more information on how to specify valid
#'        \code{SuperLearner} libraries. It is expected that the wrappers used
#'        in the library will play nicely with the input variables, which will
#'        be called \code{"trt"} and \code{names(adjustVars)}.
#' @param glm.ftime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the outcome regression (either cause-specific hazards or
#'        conditional mean). Ignored if \code{SL.ftime != NULL}. Use \code{"trt"}
#'        to specify the treatment in this formula (see examples). The formula
#'        can additionally include any variables found in
#'        \code{names(adjustVars)}.
#' @param returnModels A boolean indicating whether to return the
#'        \code{SuperLearner} or \code{glm} objects used to estimate the
#'        nuisance parameters. Must be set to \code{TRUE} if the user plans to
#'        use calls to \code{timepoints} to obtain estimates at times other than
#'        \code{t0}. See \code{?timepoints} for more information.
#' @param bounds A list of bounds... TODO: Add more description here.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats as.formula predict model.matrix optim glm binomial
#' @importFrom SuperLearner SuperLearner
#'
#'
#' @return The function returns a list that is exactly the same as the input
#'         \code{dataList}, but with additional columns corresponding to the
#'         hazard, pseudo-hazard, and the total hazard for summed over all
#'         causes \code{k < j}.
#'

estimateHazards <- function(dataList, J, adjustVars,
                            SL.ftime = NULL, glm.ftime = NULL,
                            returnModels, bounds, verbose, ...) {

  ftimeMod <- vector(mode = "list", length = length(J))
  names(ftimeMod) <- paste0("J", J)

  if(is.null(SL.ftime)) {
    if(is.null(bounds)) {
      for(j in J) {
        # formula
        Qj.form <- sprintf("%s ~ %s", paste("N", j, sep = ""), glm.ftime)

        # add up all events less than current j to see who to include in regression
        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for(i in J[J < j]) {
          NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N",i)]]
        }

        # fit GLM
        if(all(class(glm.ftime[[1]]) != "glm")) {
          Qj.mod <- glm(as.formula(Qj.form),
                        data = dataList[[1]][NlessthanJ == 0, ],
                        family = "binomial")
          Qj.mod <- cleanglm(Qj.mod)
        } else {
          Qj.mod <- glm.ftime[[paste0("J",j)]]
        }
        ftimeMod[[paste0("J",j)]] <- if(returnModels){ Qj.mod } else { NULL }

        # get predictions back
        dataList <- lapply(dataList, function(x, j) {
          suppressWarnings(
            x[[paste0("Q",j,"PseudoHaz")]] <- predict(Qj.mod, type = 'response', newdata = x)
          )
          if(j != min(J)) {
            x[[paste0("hazLessThan",j)]] <- rowSums(cbind(rep(0, nrow(x)), 
                                                          x[, paste0('Q', J[J < j], 'Haz')]))
            x[[paste0("Q",j,"Haz")]] <- x[[paste0("Q",j,"PseudoHaz")]] * (1-x[[paste0("hazLessThan",j)]])

          } else {
            x[[paste0("hazLessThan",j)]] <- 0
            x[[paste0("Q",j,"Haz")]] <- x[[paste0("Q",j,"PseudoHaz")]]
          }
          x
        }, j = j)
      }
    } else {
      for(j in J) {
        Qj.form <- sprintf("%s ~ %s", paste("N", j, sep = ""), glm.ftime)
        X <- stats::model.matrix(stats::as.formula(Qj.form),
                                 data = dataList[[1]])

        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for(i in J[J < j]) {
          NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N",i)]] 
        }

        dataList <- lapply(dataList, function(x, j) {
          if(j != min(J)) {
            x[[paste0("hazLessThan",j)]] <- rowSums(cbind(rep(0, nrow(x)), 
                                                          x[, paste0('Q', J[J < j], 'Haz')]))
          } else {
            x[[paste0("hazLessThan",j)]] <- 0
          }
          x
        }, j = j)
        Ytilde <- (dataList[[1]][[paste0("N",j)]] - dataList[[1]][[paste0("l",j)]])/
          (pmin(dataList[[1]][[paste0("u",j)]],1-dataList[[1]][[paste0("hazLessThan",j)]]) -
            dataList[[1]][[paste0("l",j)]])

        if(class("glm.ftime") != "list") {
          Qj.mod <- stats::optim(par = rep(0, ncol(X)), fn = LogLikelihood,
                                 Y = Ytilde, X = X, method = "BFGS", gr = grad,
                                 control = list(reltol = 1e-7, maxit = 50000))
        } else {
          Qj.mod <- glm.ftime[[paste0("J",j)]]
        }
        if(Qj.mod$convergence != 0) {
          stop("convergence failure")
        } else {
          beta <- Qj.mod$par
          eval(parse(text = paste0("ftimeMod$J", j, " <- Qj.mod")))
          dataList <- lapply(dataList, function(x, j) {
            newX <- stats::model.matrix(stats::as.formula(Qj.form), data = x)
            x[[paste0("Q",j,"PseudoHaz")]] <- plogis(newX %*% beta)
            x[[paste0("Q",j,"Haz")]] <- (pmin(x[[paste0("u",j)]], 1 - x[[paste0("hazLessThan",j)]]) - 
                                           x[[paste0("l",j)]])*x[[paste0("Q",j,"PseudoHaz")]] + x[[paste0("l",j)]]
            x
          },j = j)
        }
      }
    }
  } else if(is.null(glm.ftime)) {
    for (j in J) {
      # add all events less than current j to see who to include in regression
      NlessthanJ <- rep(0, nrow(dataList[[1]]))
      for (i in J[J < j]) {
        NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N",i)]]

      }

      if(class(SL.ftime[[1]]) != "SuperLearner") {
        Qj.mod <- SuperLearner(Y = dataList[[1]][[paste0("N",j)]][NlessthanJ == 0],
                               X = dataList[[1]][NlessthanJ == 0,
                                                            c('t', 'trt', names(adjustVars))],
                                        id = dataList[[1]]$id[NlessthanJ == 0],
                                        family = stats::binomial(),
                                        SL.library = SL.ftime,
                                        verbose = verbose)
      } else {
        Qj.mod <- SL.ftime[[paste0("J",j)]]
      }
      ftimeMod[[paste0("J",j)]] <- Qj.mod

      # get predictions back
      dataList <- lapply(dataList, function(x, j){
        suppressWarnings(
        x[[paste0("Q",j,"PseudoHaz")]] <- predict(Qj.mod, onlySL = TRUE,
          newdata = x[,c('t', 'trt', names(adjustVars))])[[1]]
        )
        if(j != min(J)) {
          x[[paste0("hazLessThan",j)]] <- rowSums(cbind(rep(0, nrow(x)), 
                                                        x[, paste0('Q', J[J < j], 'Haz')]))
          x[[paste0("Q",j,"Haz")]] <- x[[paste0("Q",j,"PseudoHaz")]] * 
            (1 - x[[paste0("hazLessThan",j)]])
        } else {
          x[[paste0("Q",j,"Haz")]] <- x[[paste0("Q",j,"PseudoHaz")]]
          x[[paste0("hazLessThan",j)]] <- 0
        }
        x
      }, j = j)
    }
  }
  out <- list(dataList = dataList,
              ftimeMod = if(returnModels)
                ftimeMod
              else
                NULL)
}
