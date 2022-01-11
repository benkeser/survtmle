# run from an R session inside survtmle repo
here::i_am("sandbox/att_test_sim.R")
devtools::load_all(here::here())
future.apply::plan(multicore, workers = 25)

make_data <- function(n){
	W1 <- rbinom(n, 1, 0.5)
	W2 <- rbinom(n, 1, 0.5)
	A <- rbinom(n, 1, plogis(1/2 * W1 - 1/2 * W2))
	T <- rgeom(n, plogis(-2 + -1/2 * A + 1/2 * W1 + 1/2 * W2)) + 1
	C <- rgeom(n, plogis(-3 + 1/2 * W1 + 1/2 * W2)) + 1
	T_tilde <- pmin(T, C)
	failure_ind <- as.numeric(T_tilde == T)
	return(list(
    adjustVars = data.frame(W1 = W1, W2 = W2),
    trt = A,
    ftime = T_tilde,
    ftype = failure_ind
  ))
}

get_truth <- function(n = 1e6, t0 = 4){
	W1 <- rbinom(n, 1, 0.5)
	W2 <- rbinom(n, 1, 0.5)
	A <- rbinom(n, 1, plogis(1/2 * W1 - 1/2 * W2))

	T1 <- rgeom(n, plogis(-2 + -1/2 * 1 + 1/2 * W1 + 1/2 * W2)) + 1
	T0 <- rgeom(n, plogis(-2 + -1/2 * 0 + 1/2 * W1 + 1/2 * W2)) + 1
	
	att <- mean(T1[A == 1] <= t0) - mean(T0[A == 1] <= t0)
	return(att)
}

do_one <- function(n, t0 = 4){
	data_list <- make_data(n)
	fit <- survtmle(
    adjustVars = data_list$adjustVars,
    trt = data_list$trt,
    ftime = data_list$ftime,
    ftype = data_list$ftype,
    glm.trt = "W1 + W2",
    glm.ctime = "W1 + W2",
    glm.ftime = "W1 + W2 + trt",
    att = TRUE, 
    t0 = t0
  )
  att_est <- diff(fit$est)

  g <- t(c(1, -1))
  se <- sqrt( g %*% fit$var %*% t(g) )
  att_ci <- c(att_est) + c(-1.96, 1.96) * c(se)

  return(c(att_est, att_ci))
}

n_sim <- 1000
sample_sizes <- c(500, 1000, 2000, 4000)
result_list <- vector(mode = 'list', length = length(sample_sizes))
idx <- 0
for(n in sample_sizes){
	idx <- idx + 1
	result_list[[idx]] <- future_replicate(n_sim, do_one(n))
}

summarize_result <- function(result, truth){
	bias <- mean(result[1,]) - truth
	coverage <- mean(result[2,] < truth & result[3,] > truth)
	return(c(bias = bias, coverage = coverage))
}

truth <- get_truth()
result_summary <- lapply(result_list, summarize_result, truth = truth)
result_summary_format <- data.frame(Reduce(rbind, result_summary))
result_summary_format$n <- sample_sizes
knitr::kable(result_summary_format, row.names = FALSE)

# |      bias| coverage|    n|
# |---------:|--------:|----:|
# | 0.0048499|    0.973|  500|
# | 0.0030291|    0.977| 1000|
# | 0.0029230|    0.969| 2000|
# | 0.0020602|    0.976| 4000|


make_another_data <- function(n){
	W1 <- rbinom(n, 1, 0.75)
	W2 <- rbinom(n, 1, 0.25)
	A <- rbinom(n, 1, plogis(1 * W1 - 2 * W2))
	T <- rgeom(n, plogis(-2 + -1 * A + 1/4 * W1 - 3/5 * W2)) + 1
	C <- rgeom(n, plogis(-3 - 1/2 * W1 - 3/2 * W2)) + 1
	T_tilde <- pmin(T, C)
	failure_ind <- as.numeric(T_tilde == T)
	return(list(
    adjustVars = data.frame(W1 = W1, W2 = W2),
    trt = A,
    ftime = T_tilde,
    ftype = failure_ind
  ))
}

get_another_truth <- function(n = 1e6, t0 = 5){
	W1 <- rbinom(n, 1, 0.75)
	W2 <- rbinom(n, 1, 0.25)
	A <- rbinom(n, 1, plogis(1 * W1 - 2 * W2))
	T1 <- rgeom(n, plogis(-2 + -1 * 1 + 1/4 * W1 - 3/5 * W2)) + 1
	T0 <- rgeom(n, plogis(-2 + -1 * 0 + 1/4 * W1 - 3/5 * W2)) + 1
	att <- mean(T1[A == 1] <= t0) - mean(T0[A == 1] <= t0)
	return(att)
}


do_another_one <- function(n, t0 = 5){
	data_list <- make_another_data(n)
	fit <- survtmle(
    adjustVars = data_list$adjustVars,
    trt = data_list$trt,
    ftime = data_list$ftime,
    ftype = data_list$ftype,
    glm.trt = "W1 + W2",
    glm.ctime = "W1 + W2",
    glm.ftime = "W1 + W2 + trt",
    att = TRUE, 
    t0 = t0
  )
  att_est <- diff(fit$est)

  g <- t(c(1, -1))
  se <- sqrt( g %*% fit$var %*% t(g) )
  att_ci <- c(att_est) + c(-1.96, 1.96) * c(se)

  return(c(att_est, att_ci))
}

another_result_list <- vector(mode = 'list', length = length(sample_sizes))
idx <- 0
for(n in sample_sizes){
	idx <- idx + 1
	another_result_list[[idx]] <- future_replicate(n_sim, do_another_one(n))
}

another_truth <- get_another_truth()
another_result_summary <- lapply(another_result_list, summarize_result, truth = another_truth)
another_result_summary_format <- data.frame(Reduce(rbind, another_result_summary))
another_result_summary_format$n <- sample_sizes
knitr::kable(another_result_summary_format, row.names = FALSE)

# |      bias| coverage|    n|
# |---------:|--------:|----:|
# | 0.0073568|    0.964|  500|
# | 0.0060761|    0.958| 1000|
# | 0.0032362|    0.966| 2000|
# | 0.0010769|    0.972| 4000|
