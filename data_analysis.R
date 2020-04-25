setwd("~/Documents/Privacy_git/Privacy")

library(tidyverse)
library(dplyr)
library(knitr)
library(poweRlaw)

source("functions.R")


# -------------
# California
# -------------

# show california behavior

cal1 <- read.table("Dati/campione_cal_1.txt", header=FALSE)
freqcal1 <- table(cal1)
tablefreqcal1 <- data.frame("frequency"=as.factor(names(freqcal1)), "count"=as.numeric(freqcal1))
ggplot(data=tablefreqcal1, aes(reorder(frequency,-count), log(count))) + geom_col() + ggtitle('Cal1 log freq count')
m_pl = displ$new(cal1$V1)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl) 
lines(m_pl, col=2)

cal2 <- read.table("Dati/campione_cal_2.txt", header=FALSE)
freqcal2 <- table(cal2)
tablefreqcal2 <- data.frame("frequency"=as.integer(names(freqcal2)), "count"=as.numeric(freqcal2))
ggplot(data=tablefreqcal2, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Cal2 log freq count')
m_pl = displ$new(cal2$V1)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl) 
lines(m_pl, col=2)

cal3 <- read.table("Dati/campione_cal_3.txt", header=FALSE)
freqcal3 <- table(cal3)
tablefreqcal3 <- data.frame("frequency"=as.integer(names(freqcal3)), "count"=as.numeric(freqcal3))
ggplot(data=tablefreqcal3, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Cal3 log freq count')
m_pl = displ$new(cal3$V1)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl) 
lines(m_pl, col=2)

cal4 <- read.table("Dati/campione_cal_4.txt", header=FALSE)
freqcal4 <- table(cal4)
tablefreqcal4 <- data.frame("frequency"=as.integer(names(freqcal4)), "count"=as.numeric(freqcal4))
ggplot(data=tablefreqcal4, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Cal4 log freq count')
m_pl = displ$new(cal4$V1)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl) 
lines(m_pl, col=2)

cal5 <- read.table("Dati/campione_cal_5.txt", header=FALSE)
freqcal5 <- table(cal5)
tablefreqcal5 <- data.frame("frequency"=as.integer(names(freqcal5)), "count"=as.numeric(freqcal5))
ggplot(data=tablefreqcal5, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Cal5 log freq count')
m_pl = displ$new(cal5$V1)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl) 
lines(m_pl, col=2)



# -------------
# IPMUS DATA
# -------------

prova2 <- read.csv("Dati/usa_00002.csv", header = T)
head(prova2)
prova2sub <-prova2[prova2$AGE > 20,]
prova2sub <- prova2sub %>% mutate(id = 1:n())

# ------------
# Trying different cross-classifications
# ------------

# # "COARSER": region and not state or city, 
# # race and not raced, 
# # educ and not educd
# # empstat and not empstatd,
# # vetstat and not vetstatd
# length(levels(as.factor(prova2sub$REGION))) # 9
# tablereg <- prova2sub %>% group_by(REGION, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#                                  HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                  INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_reg <- tablereg$count
# tablefreq_reg <- table(freq_reg)
# tablefreq_reg <- data.frame("frequency"=as.factor(names(tablefreq_reg)), "count"=as.numeric(tablefreq_reg))
# ggplot(data=tablefreq_reg, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('Region log frequencies counts')


# # "LESS COARSE": STATE, everything else same
# # TAKE HOME: doesn't change much from region, actually the tails are finer
# length(levels(as.factor(prova2sub$STATE))) # 51
# tablestate <- prova2sub %>% group_by(STATEICP, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#                                  HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                  INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_state <- tablestate$count
# tablefreq_state <- table(freq_state)
# tablefreq_state <- data.frame("frequency"=as.factor(names(tablefreq_state)), "count"=as.numeric(tablefreq_state))
# ggplot(data=tablefreq_state, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('State log frequencies counts')


# FINER: CITY, everything else same
length(levels(as.factor(prova2sub$CITY))) # 103
# index of the uniques in the population
ind_uniq <- prova2sub %>% 
            group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
                HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
                INCTOT, FTOTINC, VETSTAT) %>% mutate(count=n()) %>% 
            filter(count == 1) %>% select(id) 
id_uniq <- ind_uniq$id

# run these only if you want to get the plots of the frequencies counts
tablecity <- prova2sub %>% group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
                                    HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
                                    INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
freq_city <- tablecity$count
tablefreq_city <- table(freq_city)
tablefreq_city <- data.frame("frequency"=as.factor(names(tablefreq_city)), "count"=as.numeric(tablefreq_city))
ggplot(data=tablefreq_city, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('City log frequencies counts (61 classes)')
m_pl = displ$new(freq_city)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl, main="US 2018 city (61 classes)") 
lines(m_pl, col=2, main="US 2018 city (61 classes)")


# # FINER: city and raced
# # TAKE HOME: doesn't change much from CITY and race alone.
# tablecityrace <- prova2sub %>% group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACED,
#                                     HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                     INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_cityrace <- tablecityrace$count
# tablefreq_cityrace <- table(freq_cityrace)
# tablefreq_cityrace <- data.frame("frequency"=as.factor(names(tablefreq_cityrace)), "count"=as.numeric(tablefreq_cityrace))
# ggplot(data=tablefreq_cityrace, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('CityRace log frequencies counts')


# # CITY and EDUCD
# tablecityeduc <- prova2sub %>% group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#                                         HCOVANY, EDUCD, SCHLTYPE, EMPSTAT, OCC, 
#                                         INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_cityeduc <- tablecityeduc$count
# tablefreq_cityeduc <- table(freq_cityeduc)
# tablefreq_cityeduc <- data.frame("frequency"=as.factor(names(tablefreq_cityeduc)), "count"=as.numeric(tablefreq_cityeduc))
# ggplot(data=tablefreq_cityeduc, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('CityEduc log frequencies counts')


# # Without city
# tablenocity <- prova2sub %>% group_by(HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
#                                     HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                     INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_nocity <- tablenocity$count
# tablefreq_nocity <- table(freq_nocity)
# tablefreq_nocity <- data.frame("frequency"=as.factor(names(tablefreq_nocity)), "count"=as.numeric(tablefreq_nocity))
# ggplot(data=tablefreq_nocity, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('No geographic log frequencies counts')


# # Without city and VALUEH and FAMSIZE
# tablenocity1 <- prova2sub %>% group_by(HHINCOME, NCHILD, SEX, AGE, MARST, RACE,
#                                       HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
#                                       INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
# freq_nocity1 <- tablenocity1$count
# tablefreq_nocity1 <- table(freq_nocity1)
# tablefreq_nocity1 <- data.frame("frequency"=as.factor(names(tablefreq_nocity1)), "count"=as.numeric(tablefreq_nocity1))
# ggplot(data=tablefreq_nocity1, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('No geographic log frequencies counts')


# # Without city and VALUEH and FAMSIZE and SCHLTYPE and FTOTINC
# tablenocity2 <- prova2sub %>% group_by(HHINCOME, NCHILD, SEX, AGE, MARST, RACE,
#                                        HCOVANY, EDUC, EMPSTAT, OCC, 
#                                        INCTOT, VETSTAT) %>% summarize(count=n())
# freq_nocity2 <- tablenocity2$count
# tablefreq_nocity2 <- table(freq_nocity2)
# tablefreq_nocity2 <- data.frame("frequency"=as.factor(names(tablefreq_nocity2)), "count"=as.numeric(tablefreq_nocity2))
# ggplot(data=tablefreq_nocity2, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('No geographic log frequencies counts')



# # Without city and VALUEH and FAMSIZE and SCHLTYPE and FTOTINC and OCC
# tablenocity3 <- prova2sub %>% group_by(HHINCOME, NCHILD, SEX, AGE, MARST, RACE,
#                                        HCOVANY, EDUC, EMPSTAT,
#                                        INCTOT, VETSTAT) %>% summarize(count=n())
# freq_nocity3 <- tablenocity3$count
# tablefreq_nocity3 <- table(freq_nocity3)
# tablefreq_nocity3 <- data.frame("frequency"=as.factor(names(tablefreq_nocity3)), "count"=as.numeric(tablefreq_nocity3))
# ggplot(data=tablefreq_nocity3, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('No geographic log frequencies counts')


# # Without city and VALUEH and FAMSIZE and SCHLTYPE and FTOTINC and OCC and INCTOT and HCOVANY
# tablenocity4 <- prova2sub %>% group_by(HHINCOME, NCHILD, SEX, AGE, MARST, RACE,
#                                        EDUC, EMPSTAT, VETSTAT) %>% summarize(count=n())
# freq_nocity4 <- tablenocity4$count
# tablefreq_nocity4 <- table(freq_nocity4)
# tablefreq_nocity4 <- data.frame("frequency"=as.factor(names(tablefreq_nocity4)), "count"=as.numeric(tablefreq_nocity4))
# ggplot(data=tablefreq_nocity4, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('No geographic log frequencies counts (171 classes)')
# m_pl = displ$new(freq_nocity4)
# est = estimate_xmin(m_pl)
# m_pl$setXmin(est)
# plot(m_pl, main="US 2018 no geographical (171 classes)") 
# lines(m_pl, col=2, main="US 2018 no geographical (171 classes)")



# ---------------
# Now let's take freq_city and let's take some samples
# ---------------


# ---------------
# Sample 5%
# ---------------

N = dim(prova2sub)[1]
percentage0 = 0.05
n0 = as.integer(percentage0*N)
ind0 <- sample(1:N, n0, replace = FALSE)
sample_n0 <- prova2sub[ind0,]

table_n0 <- sample_n0 %>% dplyr::group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
                                   HCOVANY, EDUCD, SCHLTYPE, EMPSTAT, OCC, 
                                   INCTOT, FTOTINC, VETSTAT) %>% dplyr::summarize(count=dplyr::n())
ind_uniqsamp0 <- sample_n0 %>%
  dplyr::group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
           HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
           INCTOT, FTOTINC, VETSTAT) %>% dplyr::mutate(count=dplyr::n()) %>% 
  filter(count ==1) %>% select(id) 
id_uniqsamp0 <- ind_uniqsamp0$id

# compute tau_1 true
id_alluniq0 <- intersect(id_uniqsamp0, id_uniq)
tau1_true0 <- length(id_alluniq0)

# comparison uniques in sample AND pop, uniques in sample, uniques in population
tau1_true0
length(id_uniq)
length(id_uniqsamp0)

# m1 = uniques in sample
m1_0 <- length(id_uniqsamp0)

# compute and plot frequencies counts to input in estimation
freq_n0 <- table_n0$count
tablefreq_n0 <- table(freq_n0)
tablefreq_n0 <- data.frame("frequency"=as.factor(names(tablefreq_n0)), "count"=as.numeric(tablefreq_n0))
ggplot(data=tablefreq_n0, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('City log frequencies counts sample 5%')

# Estimates of tau1 and confidence intervals 
out_PY0 <- max_EPPF_PY(freq_n0)
tau1_PY0 <- tau1_py(m1_0, n0, out_PY0$par[1], out_PY0$par[2], N)
PY0_sim <- tau1_py_sim(freq_n0, out_PY0$par[1], out_PY0$par[2], N)
PY0_lower <- quantile(PY0_sim, 0.025)
PY0_upper <- quantile(PY0_sim, 0.975)

out_DP0 <- max_EPPF_DP(freq_n0)
tau1_DP0 <- tau1_dp(m1_0, n0, out_DP0$par, N)
DP0_sim  <- tau1_py_sim(freq_n0, out_DP0$par, 0, N)
DP0_lower <- quantile(DP0_sim, 0.025)
DP0_upper <- quantile(DP0_sim, 0.975)
# tau1_NP_pois0 <- tau1_np_pois(n0, N, freq_n0)
# tau1_NP_bin0 <- tau1_np_bin(n0, N, freq_n0)


kable(data.frame(n = n0, N = N, percentage = n0/N, 
                 K_n = dim(table_n0)[1], m1 = m1_0, 
                 true_tau1 = tau1_true0, 
                 tau1_py = tau1_PY0, CI_PY = paste("[", PY0_lower,", ", PY0_upper,"]",sep=""),
                 tau1_dp = tau1_DP0, CI_PY = paste("[", DP0_lower,", ", DP0_upper,"]",sep="")))


# ---------------
# Sample 10%
# ---------------

N = dim(prova2sub)[1]
percentage1 = 0.10
n1 = as.integer(percentage1*N)
ind1 <- sample(1:N, n1, replace = FALSE)
sample_n1 <- prova2sub[ind1,]

table_n1 <- sample_n1 %>% group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
                                    HCOVANY, EDUCD, SCHLTYPE, EMPSTAT, OCC, 
                                    INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
ind_uniqsamp <- sample_n1 %>%
  group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
                          HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
                          INCTOT, FTOTINC, VETSTAT) %>% mutate(count=n()) %>% 
  filter(count ==1) %>% select(id) 
id_uniqsamp <- ind_uniqsamp$id

# compute tau_1 true
id_alluniq <- intersect(id_uniqsamp, id_uniq)
tau1_true <- length(id_alluniq)

# comparison uniques in sample AND pop, uniques in sample, uniques in population
tau1_true # pop and sample uniques
length(id_uniq) # pop uniques
length(id_uniqsamp) # sample uniques

# number of uniques in sample
m1 <- length(id_uniqsamp)

# frequencies counts
freq_n1 <- table_n1$count

# if you wanna plot the frequency counts
tablefreq_n1 <- table(freq_n1)
tablefreq_n1 <- data.frame("frequency"=as.factor(names(tablefreq_n1)), "count"=as.numeric(tablefreq_n1))
ggplot(data=tablefreq_n1, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('City log frequencies counts sample 10%')

# Stime
out_PY <- max_EPPF_PY(freq_n1)
tau1_PY <- tau1_py(m1, n1, out_PY$par[1], out_PY$par[2], N)
PY1_sim <- tau1_py_sim(freq_n1, out_PY$par[1], out_PY$par[2], N, R=1000)
PY1_lower <- quantile(PY1_sim, 0.025)
PY1_upper <- quantile(PY1_sim, 0.975)
out_DP <- max_EPPF_DP(freq_n1)
tau1_DP <- tau1_dp(m1, n1, out_DP$par, N)
DP1_sim  <- tau1_py_sim(freq_n1, out_DP$par, 0, N)
DP1_lower <- quantile(DP1_sim, 0.025)
DP1_upper <- quantile(DP1_sim, 0.975)
# tau1_NP_pois <- tau1_np_pois(n1, N, freq_n1)
# tau1_NP_bin <- tau1_np_bin(n1, N, freq_n1)


kable(data.frame(n = n1, N = N, percentage = percentage1, 
                 K_n = dim(table_n1)[1], m1 = m1, 
                 true_tau1 = tau1_true, 
                 tau1_py = tau1_PY, CI_PY = paste("[", PY1_lower,", ", PY1_upper,"]",sep=""), 
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP1_lower,", ", DP1_upper,"]",sep="")))


# ---------------
# Sample 20%
# ---------------

n2 = as.integer(0.20*N)
ind2 <- sample(1:N, n2, replace = FALSE)
sample_n2 <- prova2sub[ind2,]
sample_n2 <- bind_cols(sample_n2, 'original_id' = ind2)
table_n2 <- sample_n2 %>% group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
                                   HCOVANY, EDUCD, SCHLTYPE, EMPSTAT, OCC, 
                                   INCTOT, FTOTINC, VETSTAT) %>% summarize(count=n())
freq_n2 <- table_n2$count
tablefreq_n2 <- table(freq_n2)
tablefreq_n2 <- data.frame("frequency"=as.factor(names(tablefreq_n2)), "count"=as.numeric(tablefreq_n2))
ggplot(data=tablefreq_n2, aes(reorder(frequency, -count), log(count))) + geom_col() + ggtitle('City log frequencies counts sample 20%')
m_pl = displ$new(freq_n2)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl, main="sample city 20%") 
lines(m_pl, col=2, main="sample city 20%")

ind2_uniqsamp <- sample_n2 %>%
  group_by(CITY, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST, RACE,
           HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC, 
           INCTOT, FTOTINC, VETSTAT) %>% mutate(count=n()) %>% 
  filter(count==1) %>% select(original_id) 
id2_uniqsamp <- ind2_uniqsamp$original_id

# tau1_true
id2_alluniq <- intersect(id2_uniqsamp, id_uniq)
tau1_true2 <- length(id2_alluniq)

# comparison uniques in sample AND pop, uniques in sample, uniques in population
tau1_true2
length(id_uniq)
length(id2_uniqsamp)

# m1
m1_2 <- length(id2_uniqsamp)

# Estimates
out_PY2 <- max_EPPF_PY(freq_n2)
tau1_PY2 <- tau1_py(m1_2, n2, out_PY2$par[1], out_PY2$par[2], N)
PY2_sim <- tau1_py_sim(freq_n2, out_PY2$par[1], out_PY2$par[2], N)
PY2_lower <- quantile(PY2_sim, 0.025)
PY2_upper <- quantile(PY2_sim, 0.975)

out_DP2 <- max_EPPF_DP(freq_n1)
tau1_DP2 <- tau1_dp(m1_2, n2, out_DP2$par, N)
DP2_sim  <- tau1_py_sim(freq_n2, out_DP2$par, 0, N)
DP2_lower <- quantile(DP2_sim, 0.025)
DP2_upper <- quantile(DP2_sim, 0.975)
# tau1_NP_pois2 <- tau1_np_pois(n2, N, freq_n2)
# tau1_NP_bin2 <- tau1_np_bin(n2, N, freq_n2)


kable(data.frame(n = n2, N = N, percentage = n2/N, 
                 K_n = dim(table_n2)[1], m1 = m1_2, 
                 true_tau1 = tau1_true2, 
                 tau1_py = tau1_PY2, CI_PY = paste("[", PY2_lower,", ", PY2_upper,"]",sep=""),
                 tau1_dp = tau1_DP2, CI_DP = paste("[", DP2_lower,", ", DP2_upper,"]",sep="")))



# -----------------
# Cal 3
# -----------------

cal3 <- read.table("Dati/campione_cal_3.txt", header=FALSE)
freqcal3 <- sort(cal3$V1)
table_cal3 <- table(cal3)
m1_cal3 = sum(cal3==1)
tau1_truecal3 = 469

ncal3 = 115093
Ncal = ncal3*10
K_ncal3 = length(table_cal3)

# Stima parametri
out_PYcal3 <- max_EPPF_PY(cal3)
tau1_PYcal3 <- tau1_py(m1_cal3, ncal3, out_PYcal3$par[1], out_PYcal3$par[2], Ncal)
out_DPcal3 <- max_EPPF_DP(cal3)
tau1_DPcal3 <- tau1_dp(m1_cal3, ncal3, out_DPcal3$par, Ncal)
cal3PY_sim <- tau1_py_sim(freqcal3, out_PYcal3$par[1], out_PYcal3$par[2], N)
cal3PY_lower <- quantile(cal3PY_sim, 0.025)
cal3PY_upper <- quantile(PY2_sim, 0.975)

kable(data.frame(n = ncal3, N = Ncal, percentage = ncal3/Ncal, 
                 K_n = K_ncal3,  K_N = 'NA', m1 = m1_cal3, 
                 true_tau1 = tau1_truecal3, 
                 tau1_py = tau1_PYcal3, tau1_dp = tau1_DPcal3))
