# Additional code

# # ------------------------
# # California dataset 
# # ------------------------

# cal3 <- read.table("data/campione_cal_3.txt", header=FALSE)
# freqcal3 <- table(cal3)
# tablefreqcal3 <- data.frame("frequency"=as.integer(names(freqcal3)), "count"=as.numeric(freqcal3))
# ggplot(data=tablefreqcal3, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Cal3 log freq count')
# m_pl = displ$new(cal3$V1)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl) 
# lines(m_pl, col=2)





# ------------
# Trying different cross-classifications (you'll see a lot of commented attempts...)
# ------------

# -----------
# FINER: CITY, everything else same
# -----------

# index of the uniques in the population. Do not uncomment the following lines!
# You can directly download the file of the frequencies uncommenting the last 2 lines.
# id_uniq <- datasetsub %>%
#   add_count(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#            HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
#            INCTOT, FTOTINC, VETSTAT) %>% filter(n==1) %>% select(id)
# id_uniq <- id_uniq$id
# write.table(id_uniq, "data/id_uniq.txt", sep="\t", col.names="id")
# id_uniq <- read.table("data/id_uniq.txt", header=TRUE)
# id_uniq <- id_uniq$id

# Frequencies counts and plots.
# do not uncomment these following lines! 
# You can directly download the file of the frequencies uncommenting the last 2 lines.
# table_city <- datasetsub %>%
#   group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#             HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
#             INCTOT, FTOTINC, VETSTAT) %>% count()
# freq_city <- table_city$n
# write.table(freq_city, "data/freq_city1.txt", sep="\t", col.names="count")
# freq_city <- read.table("data/freq_city1.txt", header=TRUE)
# freq_city <- freq_city$count

# tablefreq_city <- table(freq_city)
# tablefreq_city <- data.frame("frequency"=as.factor(names(tablefreq_city)), "count"=as.numeric(tablefreq_city))
# ggplot(data=tablefreq_city, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('City log frequencies counts (61 classes)')
# m_pl = displ$new(freq_city)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl, main="US 2018 city (61 classes)")
# lines(m_pl, col=2, main="US 2018 city (61 classes)")

# ------------
# "COARSER": region and not state or city,
# three attempts
# ------------
# tablereg <- datasetsub %>% group_by(REGION, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#                                  HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
#                                  INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_reg <- tablereg$count
# tablefreq_reg <- table(freq_reg)
# tablefreq_reg <- data.frame("frequency"=as.factor(names(tablefreq_reg)), "count"=as.numeric(tablefreq_reg))
# ggplot(data=tablefreq_reg, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Region log frequencies counts')
# m_pl = displ$new(freq_reg)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl, main="US 2018 region")
# lines(m_pl, col=2, main="US 2018 region")
# 
# tablereg6 <- datasetsub %>% group_by(REGION, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE) %>% count()
# freq_reg6 <- tablereg6$n
# m_pl6 = displ$new(freq_reg6)
# est6 = estimate_xmin(m_pl6)
# m_pl6$setXmin(est6)
# plot(m_pl6, main="US 2018 reg")
# lines(m_pl6, col=2, main="US 2018 reg")
# id_uniq6 <- datasetsub %>%
#   add_count(REGION, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE) %>% filter(n==1) %>% select(id)
# id_uniq6 <- id_uniq6$id


# ---------------------------
# other three RANDOM attemps
# ---------------------------

# table_nospace <- datasetsub %>% group_by(STATEICP, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, 
#                                         HCOVANY, EDUC, EMPSTAT, OCC,
#                                         INCTOT, FTOTINC) %>% summarize(count=n())
# freq_nospace <- table_nospace$count
# m_pl = displ$new(freq_nospace)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl, main="US 2018 no space variables")
# lines(m_pl, col=2, main="US 2018 no space variables")
# id_uniq_nospace <- datasetsub %>%
#   add_count(HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, 
#             HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
#             INCTOT, FTOTINC, VETSTAT) %>% filter(n==1) %>% select(id)
# id_uniq_nospace <- id_uniq_nospace$id
# 
# 
# table_nospace <- datasetsub %>% group_by(STATEICP, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, 
#                                         HCOVANY, EDUC, EMPSTAT, OCC,
#                                         INCTOT, FTOTINC) %>% summarize(count=n())
# freq_nospace <- table_nospace$count
# m_pl = displ$new(freq_nospace)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl, main="US 2018 no space variables")
# lines(m_pl, col=2, main="US 2018 no space variables")
# id_uniq_nospace <- datasetsub %>%
#   add_count(STATEICP, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, 
#             HCOVANY, EDUC, EMPSTAT, OCC,
#             INCTOT, FTOTINC) %>% filter(n==1) %>% select(id)
# id_uniq_nospace <- id_uniq_nospace$id
# 
# 
# tablestate <- datasetsub %>% group_by(EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                      INCTOT, FTOTINC, VETSTAT) %>% count()
# freq_state <- sort(tablestate$n)
# 
# m_pl = displ$new(freq_state)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl, main="US 2018")
# id_state <- datasetsub %>% add_count(EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                     INCTOT, FTOTINC, VETSTAT) %>% filter(n==1) %>% select(id)
# id_state <- id_state$id






# # old code when I didn't have the wonderful function loop_estimates
# N = dim(datasetsub)[1]
# percentage0 = 0.1
# n0 = as.integer(percentage0*N)
# ind0 <- sample(1:N, n0, replace = FALSE)
# sample_n0 <- datasetsub[ind0,]
# id_uniqsamp1_state <- sample_n0 %>% add_count(REGION, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST,
#                                               HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
#                                               INCTOT, FTOTINC, VETSTAT) %>% 
#                                   filter(n==1) %>% select(id)
# id_uniqsamp1_state <- id_uniqsamp1_state$id
# table_samp_state <- sample_n0 %>% group_by(REGION, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST,
#                                            HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
#                                            INCTOT, FTOTINC, VETSTAT) %>% count()
# freq_samp_state <- table_samp_state$n
# # m_pl = displ$new(freq_samp_state)
# # est = estimate_xmin(m_pl)
# # m_pl$setXmin(est)
# # plot(m_pl, main="US 2018 state SAMPLE")
# # compute tau_1 true
# id_alluniq_state <- intersect(id_uniqsamp1_state, id_state)
# tau1_true_state <- length(id_alluniq_state)
# # comparison uniques in sample AND pop, uniques in sample, uniques in population
# tau1_true_state # pop and sample uniques
# length(id_state)/N # pop uniques
# length(id_uniqsamp1_state)/n0 # sample uniques
# # number of uniques in sample
# m1_state <- sum(freq_samp_state==1)
# # # if you wanna plot the frequency counts
# # tablefreq_state <- table(freq_samp_state)
# # tablefreq_state <- data.frame("frequency"=as.factor(names(tablefreq_state)), "count"=as.numeric(tablefreq_state))
# # ggplot(data=tablefreq_state, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('City log frequencies counts sample 10%')
# # Stime
# out_PY <- max_EPPF_PY(freq_samp_state)
# tau1_PY <- tau1_py(m1_state, n0, out_PY$par[1], out_PY$par[2], N)
# PY1_sim <- tau1_py_sim(freq_samp_state, out_PY$par[1], out_PY$par[2], N)
# PY1_lower <- quantile(PY1_sim, 0.005)
# PY1_upper <- quantile(PY1_sim, 0.995)
# out_DP <- max_EPPF_DP(freq_samp_state)
# tau1_DP <- tau1_dp(m1_state, n0, out_DP$par, N)
# DP1_sim  <- tau1_py_sim(freq_samp_state, out_DP$par, 0, N)
# DP1_lower <- quantile(DP1_sim, 0.005)
# DP1_upper <- quantile(DP1_sim, 0.995)
# # tau1_NP_pois <- tau1_np_pois(n1, N, freq_n1)
# # tau1_NP_bin <- tau1_np_bin(n1, N, freq_n1)
# # Binomial approximation
# tau1_py_binom1 <- m1_state * (n0 / N)^(1 - out_PY$par[2])
# lower_py_binom1 <- qbinom(0.005, m1_state, (n0 / N)^(1 - out_PY$par[2]))
# upper_py_binom1 <- qbinom(0.995, m1_state, (n0 / N)^(1 - out_PY$par[2]))
# kable(data.frame(n = n0, N = N, percentage = percentage0, 
#                  m1 = m1_state, 
#                  true_tau1 = tau1_true_state, 
#                  K_n = dim(table_samp_state)[1],
#                  tau1_py = tau1_PY, CI_PY = paste("[", PY1_lower,", ", PY1_upper,"]",sep=""), 
#                  tau1_py_binapprox = tau1_py_binom1,
#                  CI_PY_binapprox = paste("[", lower_py_binom1,", ", upper_py_binom1,"]",sep=""),
#                  tau1_dp = tau1_DP, CI_DP = paste("[", DP1_lower,", ", DP1_upper,"]",sep="")))


# --------------
# Simulations
# --------------

# H     <- 3*10^6
# N     <- 10^6
# n     <- 10^6/5
# alpha <- 1.5
# data <- dataset_creation(n, alpha, H, n*20)
# 
# freq <- data$frequencies
# n <- data$n
# N <- data$N
# K_n = data$K_n
# m1 <- data$m1
# tau1_true <- data$true_tau1
# 
# out_PY <- max_EPPF_PY(freq)
# tau1_PY <- tau1_py(m1, n, out_PY$par[1], out_PY$par[2], N)
# PY1_sim <- tau1_py_sim(freq, out_PY$par[1], out_PY$par[2], N)
# PY1_lower <- quantile(PY1_sim, 0.025)
# PY1_upper <- quantile(PY1_sim, 0.975)
# out_DP <- max_EPPF_DP(freq)
# tau1_DP <- tau1_dp(m1, n, out_DP$par, N)
# DP1_sim  <- tau1_py_sim(freq, out_DP$par, 0, N)
# DP1_lower <- quantile(DP1_sim, 0.025)
# DP1_upper <- quantile(DP1_sim, 0.975)
# # tau1_NP_pois <- tau1_np_pois(n1, N, freq_n1)
# # tau1_NP_bin <- tau1_np_bin(n1, N, freq_n1)
# 
# # Binomial approximation
# tau1_py_binom1 <- m1 * (n / N)^(1 - out_PY$par[2])
# lower_py_binom1 <- qbinom(0.005, m1, (n / N)^(1 - out_PY$par[2]))
# upper_py_binom1 <- qbinom(0.995, m1, (n / N)^(1 - out_PY$par[2]))
# 
# m_pl = displ$new(freq)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl)
# 
# kable(data.frame(n = n, N = N, percentage = n / N, 
#                  m1 = m1, 
#                  true_tau1 = tau1_true, 
#                  K_n = K_n,
#                  tau1_py = tau1_PY, CI_PY = paste("[", PY1_lower,", ", PY1_upper,"]",sep=""), 
#                  tau1_py_binapprox = tau1_py_binom1,
#                  CI_PY_binapprox = paste("[", lower_py_binom1,", ", upper_py_binom1,"]",sep=""),
#                  tau1_dp = tau1_DP, CI_DP = paste("[", DP1_lower,", ", DP1_upper,"]",sep="")))


