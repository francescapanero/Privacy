rm(list = ls())

source("2_functions.R")

# Scenario 1 ----------------------

N <- 10^6 # The "L" is crucial, otherwise it is not recognized as an integer
n <- 10^5
m1 <- 2000

# Definizione degli intervalli di punti
alpha <- seq(0.05, 0.95, length = 200)
theta <- exp(seq(5, 15, length = 200))
# Ottenimento della griglia tramite prodotto cartesiano
parvalues <- expand.grid(alpha, theta) # Calcolo valori di verosimiglianza
llikvalues <- apply(parvalues, 1, 
                    function(x) tau1_py(m1 = m1, n = n, alpha = x[1], theta = x[2], N = N))

# Ri-organizzazione dei valori della log-verosimiglianza in forma matriciale
llikvalues <- matrix(llikvalues, nrow = length(alpha), ncol = length(theta), byrow = F)

par(mfrow=c(1,2))
contour(alpha, log(theta), llikvalues,
        xlab = expression(alpha), ylab = expression(log(theta)),
        levels = seq(from = 300, to = 1900, by = 300))
 plot(log(theta), tau1_dp(m1 = m1, n = n, theta = theta, N = N), type = "l", 
     xlab = expression(log(theta)), ylab = expression(tau[1]))

N <- 1000 # The "L" is crucial, otherwise it is not recognized as an integer
n <- 100
m1 <- 10

# Definizione degli intervalli di punti
alpha <- seq(0.1, 0.95, length = 200)
theta <- exp(seq(0, 8.5, length = 200))
# Ottenimento della griglia tramite prodotto cartesiano
parvalues <- expand.grid(alpha, theta) # Calcolo valori di verosimiglianza
llikvalues <- apply(parvalues, 1, 
                    function(x) tau1_py(m1 = m1, n = n, alpha = x[1], theta = x[2], N = N))

# Ri-organizzazione dei valori della log-verosimiglianza in forma matriciale
llikvalues <- matrix(llikvalues, nrow = length(alpha), ncol = length(theta), byrow = F)

contour(alpha, log(theta), llikvalues,
        xlab = expression(alpha), ylab = expression(log(theta)),
        levels = seq(from = 0, to = 19, by = 1))
plot(log(theta), tau1_dp(m1 = m1, n = n, theta = theta, N = N), type = "l", 
     xlab = expression(log(theta)), ylab = expression(tau[1]))
