# #######################################
# Test negative binomial estimator
# simulation in R and estimation using TMB
# #######################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_basic_estimator/")

# ### Data simulation

# number of stations
n_i <- 100

# density at location i
r <- 10 # measure of dispersion: small r causes convergence issues
d <- 5 # density mean
dens <- rgamma(n = n_i, shape = r, scale = d/r)

# mean catch at location i
rho <- 1.5
q_A <- 1
q_B <- q_A*rho
mu_A <- q_A * dens
mu_B <- q_B * dens

# catch at location i
n_A <- rpois(n_i, mu_A)
n_B <- rpois(n_i, mu_B)

# plot data
plot(1:n_i, n_A, ylim = range(c(n_A, n_B)))
points(1:n_i, n_B, col = "red", pch = 2)

# plot(n_A,n_B, main = cor(n_A, n_B))


# ### Fit model

# --------------------------------------
# ### fit a negative multinomial model
# small r and small n_i causes convergence issues
# consistent bias in estimated rho possibly due to data simulation non-agreement

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


# --------------------------------------
# ### fit a beta negative multinomial model
# only converges with an upper limit on beta distn parameters
# because alpha, beta -> Inf
# try different parameterization?
# algorithm very sensitive to data

library(TMB)

data = list(
  n_A = n_A,
  n_B = n_B)
parameters = list(
  log_alpha = 0,
  log_beta = 0,
  log_rho = rep(0, n_i))

version <- "beta_binomial"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))
obj = MakeADFun(data=data, 
                parameters=parameters,
                random=c("log_rho"),
                DLL=version,
                silent = F)
opt <- nlminb(obj$par,obj$fn,obj$gr,
              control = list(),
              lower = c(log(0.1),log(0.1)),
              upper = c(log(100),log(100)))
rep <- sdreport(obj)

t(rbind(rep$value, rep$sd))





