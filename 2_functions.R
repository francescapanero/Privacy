# library(plyr)
library(ggplot2)
library(knitr)
library(EnvStats)
library(MASS)
library(xtable)
library(gridExtra)
library(dplyr)
library(reticulate) # Library that interface R to python
use_condaenv() # Use python anaconda

random <- import("numpy.random") # Import python libraries

Rcpp::sourceCpp("3_cluster_py.cpp")

# Alternative implementation of rhyper, but we allow m to be a real number
qhyper2 <- function(p, m, n, k) {
  x <- 0:k
  lprobs <- lchoose(m, x) + lchoose(n, k - x) - lchoose(m + n, k) # choose(m, x) choose(n, k-x) / choose(m+n, k)
  lprobs <- lprobs
  probs <- exp(lprobs)
  cdf <- cumsum(probs)
  findInterval(p, sort(cdf))
}

# Alternative implementation of rhyper, but we allow m to be a real number
rhyper2 <- function(nn, m, n, k) {
  x <- 0:k
  lprobs <- lchoose(m, x) + lchoose(n, k - x) - lchoose(m + n, k) # choose(m, x) choose(n, k-x) / choose(m+n, k)
  lprobs <- lprobs - max(lprobs)
  probs <- exp(lprobs)
  sample(x, nn, replace = TRUE, prob = probs)
}

logEPPF_PY <- function(theta, alpha, frequencies) {
  if (any(alpha < 0, theta <= -alpha + 1e-04)) {
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
    lower = c(-Inf, 1e-16), upper = c(Inf, 1 - 1e-10)
  )
  return(out)
}

max_EPPF_PY_theta <- function(frequencies, alpha) {
  start <- 1 # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(theta) -logEPPF_PY(theta, alpha = alpha, frequencies = frequencies),
    lower = -Inf, upper = Inf
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

  sim <- sapply(K_sim, function(x) rhyper2(R, (theta + n) / (1 - alpha) - 1, x, m1))
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
  fit_DP <- max_EPPF_DP(frequencies)
  fit_PY <- max_EPPF_PY(frequencies)
  theta <- fit_PY$par[1]
  alpha <- fit_PY$par[2]
  theta_DP <- fit_DP$par


  M_l <- as.numeric(table(factor(frequencies, levels = 1:n)))
  # P_l <- M_l / sum(M_l)

  idx <- 1:(which.min(M_l) - 1) # which(P_l > 0)

  data_plot <- data.frame(Size = idx, M_l = M_l[idx], 
                          Theoretical_PY = expected_m_py(idx, n = n, alpha = alpha, theta = theta),
                          Theoretical_DP = expected_m_dp(idx, n = n, theta = theta_DP))
  
  p <- ggplot(data = data_plot, aes(x = Size, y = M_l)) +
    geom_point(size=0.8) +
    geom_line(aes(y = Theoretical_PY), color = "black") +
    geom_line(aes(y = Theoretical_DP), color = "black", linetype = "dotted") +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    xlab("r") +
    ylab(expression(m[r]))
  p
}


tau1_bs <- function(freq, N) {
  k <- length(freq)
  n <- sum(freq)
  bar_f <- n / k
  K_N <- N / bar_f # estimate number of categories in population
  alpha <- mle_negbin(freq) # number of successes neg bin

  # out <- fitdistr(freq, "negative binomial")
  # alpha <- out$estimate[1]
  beta <- 1 / (K_N * alpha) # probability neg bin
  tau_1_skin <- k * ((1 + n * beta) / (1 + N * beta))^(1 + alpha)
  tau_1_bet <- n * (1 + beta * N)^(-alpha - 1)
  return(c(tau_1_bet, tau_1_skin))
}


mle_negbin <- function(x) {
  tol <- 0.001
  iter <- 10^5
  a <- 0.001
  b <- 20

  k <- length(x) # number of classes
  n <- sum(x)
  fa <- sum(psigamma(x + a)) - k * psigamma(a) + k * log(1 - n / (n + k * a))
  fb <- sum(psigamma(x + b)) - k * psigamma(b) + k * log(1 - n / (n + k * b))

  while (fa * fb > 0) {
    b <- b + 100
    fb <- sum(psigamma(x + b)) - k * psigamma(b) + k * log(1 - n / (n + k * b))
  }
  i <- 1
  while (i <= iter & abs(a - b) > tol) {
    c <- (a + b) / 2
    fc <- sum(psigamma(x + c)) - k * psigamma(c) + k * log(1 - n / (n + k * c))
    if (fc * fb < 0) {
      a <- c
      fa <- fc
    }
    else {
      b <- c
      fb <- fc
    }
    i <- i + 1
  }
  if (abs(a - b) > tol) {
    print("Convergence failed")
  }
  return(c)
}

g_i <- function(i, lambda, pois_param, m) {
  if (m[i + 1] != 0) {
    return((-1)^i * exp(log(i + 1) + i * log(lambda) + log(ppois(i - 1, pois_param, lower.tail = FALSE)) + log(m[i + 1])))
  } else {
    return(0)
  }
}

tau1_np_pois <- function(N, n, freq) {
  lambda <- (N - n) / n
  m <- rep(0, n)
  a <- data.frame(freq) %>%
    group_by(freq) %>%
    summarise(count = n())
  for (j in 1:n) {
    if (j %in% a$freq) m[j] <- a$count[match(j, a$freq)]
  }
  ind_max <- tail(which(m != 0), 1)
  pois_param <- log(n / (2 * lambda - 1)) / (4 * lambda)
  ind <- c(0:(ind_max - 1))
  b <- sapply(ind, g_i, pois_param = pois_param, lambda = lambda, m = m)
  return(sum(b))
}

f_i <- function(alpha, i, cumsum) {
  return(alpha / (i - alpha) * rev(cumsum_m)[i + 1])
}

f_alpha <- function(alpha, ind_max, cumsum_m, K_n) {
  ind <- c(1:(ind_max - 1))
  a <- sapply(ind, f_i, alpha = alpha, cumsum = cumsum_m)
  a <- sum(a) - K_n
  return(a)
}

PY_alpha <- function(freq, ind_max, cumsum_m, K_n) {
  a <- uniroot(f_alpha, c(0, 1), ind_max = ind_max, cumsum_m = cumsum_m, K_n = K_n)
  return(a$root)
}

dataset_creation_zipf <- function(n, zipf_param, N) {
  
  # Old implementation using truncated zipf laws
  # probs <- 1 / (1:H)^(zipf_param)
  # points_full <- factor(sample(1:H, N, prob = probs, replace = TRUE), levels = 1:H)
  
  # Using python zipf sampler
  points_full <- factor(random$zipf(zipf_param, N))
  points_obs <- sample(points_full, n, replace = FALSE)
  
  # Observed frequencies
  freq_full <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed > 0]
  
  list(
    frequencies = freq_observed, n = sum(freq_observed), N = sum(freq_full),
    K_n = sum(freq_observed > 0), K_N = sum(freq_full > 0),
    m1 = sum(freq_observed == 1), percentage = n / N, zipf_param = zipf_param, true_tau1 = true_tau1
  )
}

dataset_creation_geom <- function(n, N, p) {
  points_full <- factor(rgeom(N, p = p))
  points_obs <- sample(points_full, n, replace = FALSE)
  
  # Observed frequencies
  freq_full <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed > 0]
  
  list(
    frequencies = freq_observed, n = sum(freq_observed), N = sum(freq_full),
    K_n = sum(freq_observed > 0), K_N = sum(freq_full > 0),
    m1 = sum(freq_observed == 1), percentage = n / N, true_tau1 = true_tau1
  )
}

dataset_creation_probs <- function(n, N, probs) {
  
  H <- length(probs)
  points_full <- factor(sample(1:H, N, prob = probs, replace = TRUE))
  points_obs <- sample(points_full, n, replace = FALSE)
  
  # Observed frequencies
  freq_full <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed > 0]
  
  list(
    frequencies = freq_observed, n = sum(freq_observed), N = sum(freq_full),
    K_n = sum(freq_observed > 0), K_N = sum(freq_full > 0),
    m1 = sum(freq_observed == 1), percentage = n / N, true_tau1 = true_tau1
  )
}
