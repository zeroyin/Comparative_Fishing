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
len_list <- seq(11, 10 + n_len, length.out = n_len)


# ----- Population Density -----#

# population density same within stratum
## use flat beta disribution
## use spatial distribution inferred from survey

dens <- rbeta(n_strat, 4,4)

# plot
hist(dens)

# length composition per stratum
## flat uniform distribution
## simulate a smooth curve
## use smoothed curve from survey

# dens_mat <- matrix(rep(dens, n_len), n_stn, n_len)

len_comp <- dbeta(len_list/70,1.8,2)
dens_mat <- dens %o% len_comp

# plot
plot(len_list, len_comp, type = "l")


# ----- Gear Catchability -----#

# cachability is more than gear selectivity 
## parametric selectivity models: logistic, lognormal etc
## designed selectivity 

s_A <-  dlnorm(len_list, meanlog = log(20), sdlog = 0.4)
    # cos(len_list/8)+1
q_A <- s_A/sum(s_A)

s_B <-  dlnorm(len_list, meanlog = log(40), sdlog = 0.3)
q_B <- s_B

#  true conversion: 
rho <- q_A/q_B

# plot
matplot(len_list, cbind(q_A, q_B), type = "l", lty = c(1,2), col = "black")
legend("topleft", lty = c(1,2), legend = c("A", "B"))


plot(len_list, q_A/(q_A+q_B), type = "l")

plot(len_list, rho, type = "l")


# ----- Catch at Length -----#

# area for scaling catch numbers
area <- 1000
N_A <- N_B <- matrix(NA, n_stn, n_len)

# Poisson sampling variation
for(i.stn in 1:n_stn){
    N_A[i.stn,] <- rpois(n = n_len, lambda = q_A*dens_mat[i.stn,]*area*rlnorm(1,-0.18,0.6))
    N_B[i.stn,] <- rpois(n = n_len, lambda = q_B*dens_mat[i.stn,]*area*rlnorm(1,-0.18,0.6))
}

# combined catch 
N <- N_A + N_B


# plot

jpeg(paste(sep = "-","res/truth_catch.jpg"),
     res = 300, width = 8, height = 10, units = "in")
par(mfrow = c(2,1))
matplot(len_list, t(N_A), type = "l", lty = "dashed", col = "black")
lines(len_list, colMeans(N_A), col = "orange")
matplot(len_list, t(N_B), type = "l", col = "black")
lines(len_list, colMeans(N_B), col = "orange")
dev.off()




# ########################################
# Fit models
# ########################################


# load TMB model
library(TMB)
version <- "betabinom_re_vector"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))


i.model = "BB7"

# model related: nll index and df index
ind_nll<- 1 + switch(
    i.model,
    "ZB2"=c(5,12), "ZB3"=c(0,5,12), 
    "BI0"=10,"BI1"=c(4,10),"BI2"=c(0,10),"BI3"=c(0,4,10),"BI4"=c(0,2,3,10),
    "BB0"=11,"BB1"=c(4,11),"BB2"=c(0,11),"BB3"=c(0,1,11),"BB4"=c(0,4,11),"BB5"=c(0,1,4,11),"BB6"=c(0,2,3,11),"BB7"=c(0,1,2,3,11)
)
ind_df <- switch(
    i.model,
    "ZB2"=4,"ZB3"=5,
    "BI0"=1,"BI1"=2,"BI2"=3,"BI3"=4,"BI4"=7,
    "BB0"=2,"BB1"=3,"BB2"=4,"BB3"=6,"BB4"=5,"BB5"=7,"BB6"=8,"BB7"=10
)


# basis and penalty matrices for cubic spline: default to 10 knots
library(mgcv)
cs <- smooth.construct(
    object = s(len, bs = "cr"),
    data = cbind(len = len_list, catch = colSums(N)) %>% as.data.frame(),
    knots = NULL
)

n_f <- 2
n_r <- cs$df - n_f
eigende <- eigen(cs$S[[1]])

# input for TMB

# vectorize
v_A <- as.vector(t(N_A))
v_B <- as.vector(t(N_B))
v_index <- which(v_A>0 | v_B>0)
data = list(
    A = v_A[v_index],
    B = v_B[v_index],
    indstn = ceiling(v_index/n_len)-1,
    indlen = rep(1:n_len, n_stn)[v_index]-1,
    offset = rep(0, length(v_index)),
    Xf = cs$X %*% eigende$vectors[,1:n_f+n_r],
    Xr = cs$X %*% eigende$vectors[,1:n_r],
    d = eigende$value[1:n_r]
)

parameters = list(
    beta = rep(0, n_f),
    b = rep(0, n_r),
    gamma = rep(0, n_f),
    g = rep(0, n_r),
    delta = matrix(0, n_stn, n_f),
    epsilon = matrix(0, n_stn, n_r),
    log_s_b = log(10),
    log_s_g = log(10),
    log_s_epsilon = log(10),
    chol_delta = c(1,1,2)
)
map = list(
    g = rep(factor(NA),n_r)
)

obj = MakeADFun(data=data,
                parameters=parameters,
                map=map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = T)
opt <- nlminb(obj$par,obj$fn,obj$gr)


# check convergence, maximum gradient and positive definite
# if(exists("opt")){
#     gra <- obj$gr()
#     hes <- eigen(optimHess(par=opt$par, fn=obj$fn, gr=obj$gr))$values
#     if(max(abs(gra)) < 0.1 & min(hes) > -0.1){
        aic <- 2*sum(obj$report()$nll) + 2*ind_df
        rep <- try(sdreport(obj))
        
        # estimate and std of mu and phi and rho
        est.mu <- matrix(NA,n_len,n_stn); est.mu[v_index] <- obj$report()$mu; est.mu <- t(est.mu)
        
        
        
        jpeg(paste(sep = "-","res/est_mu.jpg"),
             res = 300, width = 8, height = 6, units = "in")
        matplot(len_list, t(N_A/N), type = "p", cex = 0.2, col = "black")
        lines(len_list, colMeans(N_A/N, na.rm = T), col = "orange")
        lines(len_list, boot::inv.logit(log(rho)), col = "black")
        lines(len_list, colMeans(est.mu, na.rm = T), ylim = c(0,1), col = "blue")
        legend("bottomright", pch = c(1,NA,NA), lty = c(0,1,1), col = c("black","orange","black"), c("mu by stn", "mean obs mu", "true mu"))    
        dev.off()
        
        
        jpeg(paste(sep = "-","res/est_mu_stn.jpg"),
             res = 300, width = 12, height = 10, units = "in")
        par(mfrow = c(5,4))
        for(i in 1:n_stn){
            plot(len_list, N_A[i,]/N[i,], type = "p", cex = 0.2, col = "black")
            lines(len_list, est.mu[i,], col = "gray")
        }
        par(mfrow = c(1,1))
        dev.off()
        
        
#     }
# }

