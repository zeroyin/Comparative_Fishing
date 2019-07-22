# ###############################################
# Simulate data for model testing
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(ggplot2)


# ----- Presets -----#

# strata list
n_strat <- 20
strat_list <- seq(1, n_strat)

# station list
n_stn <- 20
# stn_strat <- 1 + rmultinom(1, n_station-n_strat, rep(1/n_strat,n_strat)) # random number of stations per stratum
n_stn_strat <- rep(1, n_stn) # fixed number of stations per stratum
stn_list <- seq(1, n_stn)

# length list
n_len <- 50
len_list <- seq(10, 10 + n_len, length.out = n_len)


# ----- Population Density -----#

# population density same within stratum
## use flat beta disribution
## use spatial distribution inferred from survey

dens <- rbeta(n_strat, 8, 5)

# plot
hist(dens)

# length composition per stratum
## flat uniform distribution
## simulate a smooth curve
## use smoothed curve from survey

dens_mat <- matrix(rep(dens, n_len), n_stn, n_len)

# len_comp <- dbeta(len_list/70,5,4)
# dens_mat <- dens %o% len_comp

# plot
plot(len_list, len_comp, type = "l")


# ----- Gear Catchability -----#

# cachability is more than gear selectivity 
## parametric selectivity models: logistic, lognormal etc
## designed selectivity 

s_A <- dlnorm(len_list, meanlog = log(40), sdlog = 0.4)
q_A <- s_A/sum(s_A)

s_B <-  dlnorm(len_list, meanlog = log(30), sdlog = 0.3)
q_B <- s_B/sum(s_B)

#  true conversion: 
rho <- q_A/q_B

# plot
matplot(len_list, cbind(q_A, q_B), type = "l", lty = c(1,2), col = "black")
legend("topleft", lty = c(1,2), legend = c("A", "B"))

plot(len_list, rho, type = "l")


# ----- Catch at Length -----#

# Poisson sampling variation
area <- 10000
N_A <- N_B <- matrix(NA, n_stn, n_len)
for(i.stn in 1:n_stn){
    N_A[i.stn,] <- rpois(n = n_len, lambda = q_A*dens_mat[i.stn,]*area)
    N_B[i.stn,] <- rpois(n = n_len, lambda = q_B*dens_mat[i.stn,]*area)
}

# plot
matplot(len_list, t(N_A), type = "l", lty = "dashed", col = "black")
lines(len_list, colMeans(N_A), col = "orange")
matplot(len_list, t(N_B), type = "l", col = "black")
lines(len_list, colMeans(N_B), col = "orange")

matplot(len_list, t(N_A/(N_A+N_B)), type = "p", cex = 0.2, col = "black")
lines(len_list, colMeans(N_A/(N_A+N_B)), col = "orange")
lines(len_list, boot::inv.logit(log(rho)), col = "black")
legend("topleft", pch = c(1,NA,NA), lty = c(0,1,1), col = c("black","black","orange"), c("mu by stn", "mean obs mu", "true mu"))    
    
