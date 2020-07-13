# library(plyr)
library(ggplot2)
library(knitr)

Rcpp::sourceCpp("3_cluster_py.cpp")


logEPPF_PY <- function(theta, alpha, frequencies) {
  if (any(alpha < 0, theta <= -alpha)) {
    return(-Inf)
  }


  # Sample size
  n <- sum(frequencies)
  # Number of frequencies
  K <- length(frequencies)
  # Loglikelihood
  loglik <- sum(log(theta + 1:(K - 1) * alpha)) - lgamma(theta + n) + lgamma(theta + 1) + sum(lgamma(frequencies - alpha)) - K * lgamma(1 - alpha)

  loglik
}

logEPPF_DP <- function(theta, frequencies) {

  # Sample size
  n <- sum(frequencies)
  # Number of frequencies
  K <- length(frequencies)

  # Loglikelihood
  loglik <- K * log(theta) - lgamma(theta + n) + lgamma(theta)

  loglik
}

max_EPPF_PY <- function(frequencies) {
  start <- c(1, 0.5) # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_PY(theta = param[1], alpha = param[2], frequencies = frequencies),
    lower = c(-Inf, -Inf), upper = c(Inf, 1 - 1e-10)
  )
  return(out)
}


max_EPPF_DP <- function(frequencies) {
  start <- 1 # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_DP(theta = param, frequencies = frequencies),
    lower = 1e-10, upper = Inf
  )
  return(out)
}


tau1_dp <- function(m1, n, theta, N) {
  m1 * (n + theta - 1) / (N + theta - 1)
}

tau1_py <- function(m1, n, theta, alpha, N) {
  # exp(log(m1) + lgamma(theta - 1 + alpha + N) - lgamma(theta + n - 1 + alpha) - lgamma(theta + N) + lgamma(theta + n))
  exp(log(m1) + lgamma(theta + alpha + N) - lgamma(theta + n - 1 + alpha) - lgamma(theta + N + 1) + lgamma(theta + n)
    + log(theta + N) - log(theta + N + alpha - 1))
}


tau1_py_sim <- function(frequencies, theta, alpha, N, R = 100, verbose = TRUE) {

  # Main quantities
  n <- sum(frequencies)
  m1 <- sum(frequencies == 1)
  K_n <- length(frequencies)

  K_sim <- plyr::aaply(matrix(0, R, 2), 1, function(x) cluster_py_C(N - n, 1 - alpha, theta + n), .progress = "text")
  # K_sim <- replicate(R,cluster_py_C(N - n, 1 - alpha, theta+n))

  sim <- rhyper(R, (theta + n) / (1 - alpha) - 1, K_sim, m1)
  sim

  cat("Estimated mean: ", round(mean(sim), 2), ". Monte Carlo se: ", round(sd(sim) / sqrt(R), 2), ".\n", sep = "")
  sim
}


expected_cl_py <- function(n, alpha, theta) {
  n <- as.integer(n)
  if (alpha == 0) {
    out <- theta * sum(1 / (theta - 1 + 1:n))
  } else {
    out <- 1 / alpha * exp(lgamma(theta + alpha + n) - lgamma(theta + alpha) - lgamma(theta + n) + lgamma(theta + 1)) - theta / alpha
  }

  return(out)
}


expected_m_dp <- function(m, n, theta) {
  out <- log(theta) + lchoose(n, m) + lgamma(m) + lgamma(theta + n - m) - lgamma(theta + n)
  exp(out)
}
expected_m_dp <- Vectorize(expected_m_dp, vectorize.args = "m")


expected_m_py <- function(m, n, alpha, theta) {
  out <- log(theta) + lchoose(n, m) + lgamma(m - alpha) - lgamma(1 - alpha) - lgamma(theta + n) + lgamma(theta) + lgamma(theta + alpha + n - m) - lgamma(theta + alpha)
  exp(out)
}
expected_m_py <- Vectorize(expected_m_py, vectorize.args = "m")

frequency_check_PY <- function(frequencies) {
  n <- sum(frequencies)
  fit_PY <- max_EPPF_PY(frequencies)
  sigma <- fit_PY$par[2]


  M_l <- as.numeric(table(factor(frequencies, levels = 1:n)))
  P_l <- M_l / sum(M_l)

  idx <- which(P_l > 0)

  data_plot <- data.frame(Size = idx, P_l = P_l[idx], Theoretical = exp(log(sigma) + lgamma(idx - sigma) - lgamma(1 - sigma) - lfactorial(idx)))
  p <- ggplot(data = data_plot, aes(x = Size, y = P_l)) +
    geom_point() +
    geom_line(aes(y = Theoretical), color = "blue", linetype = "dashed") +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    xlab(expression(M[ln])) +
    ylab(expression(M[ln] / K[n]))
  p
}