library(tidyverse)
library(ggpubr)

rm(list = ls())

load("data/IPMUS.RData")
source("2_functions.R")

# -------------------------------
# 5% dataset
# -------------------------------

dataset <- data_5perc

# MODEL CHECKING ---------------

out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

# Comparison between M_l and the asymptotic formula
# frequency_check_PY(dataset$frequencies)

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# Table creation
tab <- rbind(
  PY = expected_m_py(1:10, dataset$n, out_PY$par[2], out_PY$par[1]),
  DP = expected_m_dp(1:10, dataset$n, out_DP$par[1]),
  Data = M_l[1:10]
)

colnames(tab) <- 1:10
kable(tab, digits = 0)
xtable(tab, digits=0)


p_check <- frequency_check_PY(dataset$frequencies)
p_check <- p_check + ggtitle(paste("Dataset percentage: 5%")) + theme(plot.title = element_text(hjust = 0.5, size = 10))
p_check

# ggsave(p_check, height = 5, width = 4 * 2.5, file = "check_5.eps", device = "eps")


# ---------------------------
# PY estimation
# ---------------------------

alpha <- 0.01 # Credible intervals percentage

out_PY <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)

PY_lower <- quantile(PY_sim, alpha / 2)
PY_upper <- quantile(PY_sim, 1 - alpha / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper2(alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper2(1 - alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Bethlehem and Skinner estimators
estim <- tau1_bs(dataset$frequencies, dataset$N)
tau1_bet <- estim[1]
tau1_skin <- estim[2]
tau1_cam <- tau1_np_pois(dataset$N, dataset$n, dataset$frequencies)

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = ""),
  tau1_bet = tau1_bet, tau1_skin = tau1_skin, tau1_cam = tau1_cam
))


# Summary rounded
result <- data.frame(
  N = dataset$N,
  percentage = round(dataset$percentage * 100, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = round(tau1_PY, 0), CI_PY = paste("[", round(PY_lower, 0), ", ", round(PY_upper, 0), "]", sep = ""),
  tau1_dp = round(tau1_DP, 0), CI_DP = paste("[", round(DP_lower, 0), ", ", round(DP_upper, 0), "]", sep = ""),
  tau1_bet = round(tau1_bet, 0), tau1_skin = round(tau1_skin, 0), tau1_cam = round(tau1_cam, 0)
)

# change col names
result %>%
  kable(., col.names = c("N", "% sample", "# classes", "# sample uniques", "tau_1 true", "PY", "CI PY", "DP", "CI DP", "Bethelehem", "Skinner", "Camerlenghi"))

# Latex table
result %>%
  kable(., "latex",
    escape = F, booktabs = T, linesep = "", align = "c",
    col.names = c("N", "$\\%$ sample", "$\\#$ classes", "# sample uniques", "$\\tau_1$ true", "PY", "CI PY", "DP", "CI DP", "Bethelehem", "Skinner", "Camerlenghi")
  )

# Plot estimates and confidence intervals
type <- c("PY", "DP", "B", "S", "C")
estimates <- c(tau1_PY, tau1_DP, tau1_bet, tau1_skin, tau1_cam)
df <- data.frame(
  type = factor(type, levels = type[order(estimates)]), estim = estimates,
  lower_CI = c(PY_lower, DP_lower, NA, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA, NA)
)

p <- ggplot(df, aes(type, estim)) +
  scale_x_discrete(limits = c("B", "DP", "PY", "S", "C"))
p <- p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() +
  theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle(paste("5%")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
p
ggsave(p, file = "img/5perc.eps", device = "eps")

# -------------------------------
# 10% dataset
# -------------------------------

dataset <- data_10perc

# MODEL CHECKING ---------------

out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

# Comparison between M_l and the asymptotic formula
# frequency_check_PY(dataset$frequencies)

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# Table creation
tab <- rbind(
  PY = expected_m_py(1:10, dataset$n, out_PY$par[2], out_PY$par[1]),
  DP = expected_m_dp(1:10, dataset$n, out_DP$par[1]),
  Data = M_l[1:10]
)

colnames(tab) <- 1:10
kable(tab, digits = 0)
xtable(tab, digits=0)

p1_check <- frequency_check_PY(dataset$frequencies)
p1_check <- p1_check + ggtitle(paste("Dataset percentage: 10%")) + theme(plot.title = element_text(hjust = 0.5, size = 10))
p1_check
# ggsave(p1_check, height = 5, width = 4 * 2.5, file = "img/check_10.eps", device = "eps")

# Figure 3 in the paper
p_check_tog <- list()
p_check_tog[[1]] <- p_check
p_check_tog[[2]] <- p1_check
a <- do.call(grid.arrange, c(p_check_tog, ncol = 2))

ggsave(a, height = 4, width = 10, file = "img/check.eps", device = "eps")

# ---------------------------
# PY estimation
# ---------------------------

alpha <- 0.01 # Credible intervals percentage

out_PY <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)

PY_lower <- quantile(PY_sim, alpha / 2)
PY_upper <- quantile(PY_sim, 1 - alpha / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper2(alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper2(1 - alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Bethlehem and Skinner estimators
estim <- tau1_bs(dataset$frequencies, dataset$N)
tau1_bet <- estim[1]
tau1_skin <- estim[2]
tau1_cam <- tau1_np_pois(dataset$N, dataset$n, dataset$frequencies)

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = "")
))

# Summary rounded
result <- data.frame(
  N = dataset$N,
  percentage = round(dataset$percentage * 100, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = round(tau1_PY, 0), CI_PY = paste("[", round(PY_lower, 0), ", ", round(PY_upper, 0), "]", sep = ""),
  tau1_dp = round(tau1_DP, 0), CI_DP = paste("[", round(DP_lower, 0), ", ", round(DP_upper, 0), "]", sep = ""),
  tau1_bet = round(tau1_bet, 0), tau1_skin = round(tau1_skin, 0), tau1_cam = round(tau1_cam, 0)
)

# change col names
result %>%
  kable(., col.names = c("N", "% sample", "# classes", "# sample uniques", "tau_1 true", "PY", "CI PY", "DP", "CI DP", "Bethelehem", "Skinner", "Camerlenghi"))

# Latex table
result %>%
  kable(., "latex",
    escape = F, booktabs = T, linesep = "", align = "c",
    col.names = c("N", "$\\%$ sample", "$\\#$ classes", "# sample uniques", "$\\tau_1$ true", "PY", "CI PY", "DP", "CI DP", "Bethelehem", "Skinner", "Camerlenghi")
  )

# Plot estimates and confidence intervals
type <- c("PY", "DP", "B", "S", "C")
estimates <- c(tau1_PY, tau1_DP, tau1_bet, tau1_skin, tau1_cam)
df <- data.frame(
  type = factor(type, levels = type[order(estimates)]), estim = estimates,
  lower_CI = c(PY_lower, DP_lower, NA, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA, NA)
)

p1 <- ggplot(df, aes(type, estim)) +
  scale_x_discrete(limits = c("B", "DP", "PY", "S", "C"))
p1 <- p1 + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() +
  theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle(paste("10%")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
p1

# Figure 2 in the paper
p_together <- list()
p_together[[1]] <- p
p_together[[2]] <- p1
b <- do.call(grid.arrange, c(p_together, ncol = 2))
ggsave(b, height = 5, width = 4 * 2.5, file = "ipmus", device = "eps")
