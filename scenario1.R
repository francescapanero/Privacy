rm(list=ls())
library(knitr)
library(poweRlaw)

source('functions.R')

# -------------------------------------------
# Scenario 1 - Low sigma
# -------------------------------------------
   
H     <- 3*10^6
N     <- 10^6
n     <- 10^6/5
sigma <- 1.5

dataset <- dataset_creation(n = n, alpha = sigma, H = H, N = N)
# ADDED FRA
tablepl <- table(dataset$frequencies)
tablepl <- data.frame("frequency"=as.integer(names(tablepl)), "count"=as.numeric(tablepl))
ggplot(data=tablepl, aes(frequency, log(count))) + geom_col() + ggtitle('pl')
ggplot(data=tablepl, aes(frequency, count)) + geom_col() + ggtitle('pl')
datapl <- data.frame('x' = tablepl$frequency, 'y' =1/(tablepl$frequency)^(sigma)*n)
data <- bind_cols(tablepl, datapl)
ggplot(data = data) + geom_point(aes(x,y), color='red') + geom_col(aes(frequency, count))  + ggtitle('pl')
m_pl = displ$new(dataset$frequencies)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl, main="Power law 1.5")
lines(m_pl, col=2, main="Power law 1.5")


kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 sigma = dataset$sigma))

# CA application
dataset  <- read.table("Dati/campione_cal_1.txt", header=FALSE)
dataset <- list( frequencies = dataset$V1[dataset$V1 > 0])
dataset$n <- sum(dataset$frequencies)
dataset$percentage <- 0.1
dataset$N <- 1150934
dataset$K_n <- length(dataset$frequencies)
dataset$K_N <- NA
dataset$H <- 3600000
dataset$true_tau1 <- 15
dataset$sigma <- NA
dataset$m1 <- sum(dataset$frequencies == 1)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 sigma = dataset$sigma))

# -----------------------
# 1 ] DP estimation
# -----------------------

max_EPPF_DP(dataset$frequencies)

# Theta estimation via MLE
theta_hat <- max_EPPF_DP(dataset$frequencies)$par

# Estimation of tau1_dp
dataset$tau1_dp <- tau1_dp(dataset$m1, dataset$n, theta_hat,dataset$N)

dataset$true_tau1
dataset$tau1_dp

# -----------------------
# 2 ] PY estimation
# -----------------------

param_hat <- max_EPPF_PY(dataset$frequencies)
theta_hat <- param_hat$par[1]
alpha_hat <- param_hat$par[2]

dataset$tau1_py <- tau1_py(dataset$m1, dataset$n, theta_hat, alpha_hat, dataset$N)

dataset$true_tau1
dataset$tau1_py 


#rm(list=setdiff(ls(), "dataset"))
kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n = dataset$K_n, K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 = dataset$true_tau1, tau_dp = dataset$tau1_dp, 
                 tau_py = dataset$tau1_py, #tau_dm = dataset$tau1_dm, 
                 sigma = dataset$sigma))

save.image("scenario1.RData")
