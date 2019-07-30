# #######################################
# Test beta negative binomial estimator
# simulation in R and estimation using TMB
# #######################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test")

# ### Data simulation

# number of stations
n_i <- 100

# density across locations
r <- 10 # measure of dispersion
d <- 30 # density mean
dens <- rgamma(n = n_i, shape = r, scale = d/r)

# conditional prob of catch by B across stations
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


# plot data
plot(1:n_i, n_A, ylim = range(c(n_A, n_B)))
points(1:n_i, n_B, col = "red", pch = 2)

# plot correlation 
plot(n_A,n_B, main = cor(n_A, n_B))


# ### Fit model

# --------------------------------------
# ### fit a negative multinomial model
# simulation model = estimation model
# well recovered

library(TMB)

data = list(
    n_A = n_A,
    n_B = n_B)
parameters = list(
    log_d = 0,
    log_r = 0,
    log_dens = rep(0, n_i),
    log_q_A = 0,
    log_q_B = 0)

version <- "negative_multinomial"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))
obj = MakeADFun(data=data,
                parameters=parameters,
                map=list("log_q_A" = factor(NA)),
                random=c("log_dens"),
                DLL=version,
                silent = F)
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)

c(rep$value, rep$sd)


# plot results
plot(exp(rep$par.random), dens)
abline(0, 1, col = "red")



