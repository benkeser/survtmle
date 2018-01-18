### Could write a wrapper that takes a super learner object for
### each ftype and computes log sub-density cv risk to determine 
### weights?

# method_CC.haznloglik.R <- function () 
# {
#     computePred = function(predY, coef, control, ...) {
#         if (sum(coef != 0) == 0) {
#             stop("All metalearner coefficients are zero, cannot compute prediction.")
#         }
#         plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*% 
#             matrix(coef[coef != 0]))
#     }
#     computeCoef = function(Z, Y, libraryNames, obsWeights, control, 
#         verbose, id, ...) {
#         dupCols <- which(duplicated(Z, MARGIN = 2))
#         anyDupCols <- length(dupCols) > 0
#         modZ <- Z
#         if (anyDupCols) {
#             warning(paste0(paste0(libraryNames[dupCols], collapse = ", "), 
#                 " are duplicates of previous learners.", " Removing from super learner."))
#             modZ <- modZ[, -dupCols]
#         }
#         modlogitZ <- trimLogit(modZ, control$trimLogit)
#         logitZ <- trimLogit(Z, control$trimLogit)
#         # for each id, compute likelihood at observed failure times
#         modZ_id <- cbind(id, modZ)
#         lik_at_obs_failure <- by(modZ_id, id, function(z){
#             len_haz <- length(z[,1])
#             apply(z, 2, function(haz){
#                 S.tminus1 <- c(1, cumprod(1 - haz[1:(len_haz-1)]))
#                 dens.t <- sum(haz * S.tminus1)
#             })
#         })

#         cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights * 
#             ifelse(Y, plogis(x, log.p = TRUE), plogis(x, log.p = TRUE, 
#                 lower.tail = FALSE))))
#         names(cvRisk) <- libraryNames
#         obj_and_grad <- function(y, x, w = NULL) {
#             y <- y
#             x <- x
#             function(beta) {
#                 xB <- x %*% cbind(beta)
#                 loglik <- y * plogis(xB, log.p = TRUE) + (1 - 
#                   y) * plogis(xB, log.p = TRUE, lower.tail = FALSE)
#                 if (!is.null(w)) 
#                   loglik <- loglik * w
#                 obj <- -2 * sum(loglik)
#                 p <- plogis(xB)
#                 grad <- if (is.null(w)) 
#                   2 * crossprod(x, cbind(p - y))
#                 else 2 * crossprod(x, w * cbind(p - y))
#                 list(objective = obj, gradient = grad)
#             }
#         }
#         lower_bounds = rep(0, ncol(modZ))
#         upper_bounds = rep(1, ncol(modZ))
#         if (anyNA(cvRisk)) {
#             upper_bounds[is.na(cvRisk)] = 0
#         }
#         r <- nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)), 
#             eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds, 
#             ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) - 
#                 1), eval_jac_g_eq = function(beta) rep(1, length(beta)), 
#             opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
#         if (r$status < 1 || r$status > 4) {
#             warning(r$message)
#         }
#         coef <- r$solution
#         if (anyNA(coef)) {
#             warning("Some algorithms have weights of NA, setting to 0.")
#             coef[is.na(coef)] <- 0
#         }
#         if (anyDupCols) {
#             ind <- c(seq_along(coef), dupCols - 0.5)
#             coef <- c(coef, rep(0, length(dupCols)))
#             coef <- coef[order(ind)]
#         }
#         coef[coef < 1e-04] <- 0
#         coef <- coef/sum(coef)
#         out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
#         return(out)
#     }
#     list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
# }