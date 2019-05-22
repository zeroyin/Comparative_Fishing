# #################################
# Test beta negative binomial
# simulation in R and estimation using TMB
# #################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Workspace/test")

# ### Data simulation

# number of stations
n_i <- 100

# density across locations
r <- 10 # measure of dispersion
d <- 30 # density mean
dens <- rgamma(n = n_i, shape = r, scale = d/r)

# cond prob of catch by B across stations
alpha = 20
beta = 20
p <- rbeta(n_i, alpha, beta)
rho <- p/(1-p)

# catch at location i: poisson
q_A <- rep(1, n_i)
q_B <- q_A * rho
mu_A <- dens * q_A
mu_B <- dens * q_B
n_A <- rpois(n_i, mu_A)
n_B <- rpois(n_i, mu_B)



plot(1:n_i, n_A, ylim = range(c(n_A, n_B)))
points(1:n_i, n_B, col = "red", pch = 2)

plot(n_A,n_B, main = cor(n_A, n_B))


# ### Fit model



