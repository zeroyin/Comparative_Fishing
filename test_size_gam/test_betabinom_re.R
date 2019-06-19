# #########################################
# test beta binomial model:
# length GAM, with random effects
# model in TMB
# #########################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)
library(ggplot2)

load("data-NED2013.RData")

# organize data for model

i.species <- 23
b.len <- 1

d <- d.length %>% 
    filter(species == i.species) %>%
    transmute(
        station = factor(station),
        gear = factor(gear), 
        len = (floor(len/b.len)*b.len+b.len/2), # length grouping
        catch = catch) %>%
    group_by(station, gear, len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup() %>%
    complete(len, # remove using full_seq()
             station, 
             gear, 
             fill = list(catch = 0)) # complete zero observations
    

# length at center of each bin, vector for plotting
# lenseq <- full_seq(d$len, b.len)
lenseq <- unique(d$len)


# data for offset: same within a tow
# use zeros as this data is tow-standardized
d.offset <- d %>% 
    distinct(station) %>%
    mutate(offset = 0)


# basis and penalty matrices for cubic spline
# over length bins
# location of knots might influence smooth functions (and convergence)
# default to 10 knots resulting in 2 fixed and 8 random effects
library(mgcv)
cs <- smooth.construct(
    object = s(len, bs = "cr"),
    data = d %>% group_by(len) %>% summarise(catch = sum()),
    knots = NULL
)

# cs <- smooth.construct(
#     object = s(len, bs = "cr", k = 6),
#     data = d %>% group_by(len) %>% summarise(catch = sum()),
#     knots =  list(len = seq(min(lenseq),max(lenseq),length.out = 6)))

n_f <- 2
n_r <- cs$df - n_f
eigende <- eigen(cs$S[[1]])

# input for TMB
nlen = length(lenseq)
nstation = nlevels(d$station)

data = list(
    A = d %>% filter(gear == 9 ) %>% spread(len, catch) %>% select(-station, -gear) %>% as.matrix(),
    B = d %>% filter(gear == 15) %>% spread(len, catch) %>% select(-station, -gear) %>% as.matrix(),
    offset = outer(d.offset$offset,rep(1,length(lenseq))),
    Xf = cs$X %*% eigende$vectors[,1:n_f+n_r],
    Xr = cs$X %*% eigende$vectors[,1:n_r],
    d = eigende$value[1:n_r]
)
parameters = list(
    beta = rep(0, n_f),
    b = rep(0, n_r),
    gamma = rep(0, n_f),
    g = rep(0, n_r),
    delta = matrix(0, nstation, n_f),
    epsilon = matrix(0, nstation, n_r),
    log_s_b = log(10),
    log_s_g = log(10),
    log_s_epsilon = log(10),
    chol_delta = c(1,0,1) # use chol decomp in vector form
)
map <- list(
    # log_s_b = factor(NA),
    # log_s_g = factor(NA),
    # log_s_epsilon = factor(NA),
    # chol_delta = rep(factor(NA), length(parameters$chol_delta))
)

# run TMB model
library(TMB)
version <- "test_betabinom_re"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)

opt


# plot results

# estimate and std of mu and phi and rho
est <- summary(rep, "report")[,"Estimate"]
std <-  summary(rep, "report")[,"Std. Error"]
est.mean_mu <- est[names(est) == "mean_mu"]
std.mean_mu <- std[names(std) == "mean_mu"]
est.mean_phi <- est[names(est) == "mean_phi"]
std.mean_phi <- std[names(std) == "mean_phi"]
est.mean_log_rho <- est[names(est) == "mean_log_rho"]
std.mean_log_rho <- std[names(std) == "mean_log_rho"]
est.mu <- matrix(est[names(est) == "mu"], nrow = nstation, ncol = nlen)
std.mu <- matrix(std[names(std) == "mu"], nrow = nstation, ncol = nlen)
est.rho <- matrix(est[names(est) == "rho"], nrow = nstation, ncol = nlen)
std.rho <- matrix(std[names(std) == "rho"], nrow = nstation, ncol = nlen)



# estimated mu, phi and rho for each station
jpeg(paste(sep = "-","beta-binom-re/NED2013/estimates","species",i.species,"lenbin",b.len,"CI95_zscore.jpg"),
     res = 300, width = 6, height = 10, units = "in")
par(mfrow=c(3,1))
plot(lenseq, est.mean_mu, ylim = c(0,1), type = "l")
for(i in 1:nstation){lines(lenseq, obj$report()$mu[i,], col = "gray")}
lines(lenseq, est.mean_mu + 1.96*std.mean_mu, col = "blue", lty = "dashed")
lines(lenseq, est.mean_mu - 1.96*std.mean_mu, col = "blue", lty = "dashed")
plot(lenseq, est.mean_phi, type = "l")
plot(lenseq, est.mean_log_rho, ylim = c(-4,4), type = "l")
abline(a = 0, b = 0, col = "red")
for(i in 1:nstation){lines(lenseq, obj$report()$eta_mu[i,], col = "gray")}
lines(lenseq, est.mean_log_rho + std.mean_log_rho, col = "blue", lty = "dashed")
lines(lenseq, est.mean_log_rho - std.mean_log_rho, col = "blue", lty = "dashed")
dev.off()


# estimated mu with observations for each station
# CI is for observation
# jpeg(paste0("beta-binom-re/mu_by_station-species_", i.species, ".jpg"),
#     res = 300, width = 12, height = 10, units = "in")
# par(mfrow=c(4,4))
# for(i in 1:nstation){
#     plot(lenseq, data$A[i,]/(data$A[i,]+data$B[i,]), ylim = c(0, 1), ylab = "Prop. of Gear 9", main = paste0("Station ", levels(d$station)[i]))
#     lines(lenseq, est.mu[i,])
#     lines(lenseq, gamlss.dist::qBB(0.05, mu = est.mu[i,], sigma = obj$report()$phi[i,], bd = 1000)/1000, col = "blue")
#     lines(lenseq, gamlss.dist::qBB(0.95, mu = est.mu[i,], sigma = obj$report()$phi[i,], bd = 1000)/1000, col = "blue")
# }
# dev.off()


# estimated mu with observations for each station
# SD is from TMB sd
jpeg(paste(sep = "-","beta-binom-re/NED2013/mu_by_station","species",i.species,"lenbin",b.len,"TMB_SD.jpg"),
     res = 300, width = 24, height = 20, units = "in")
par(mfrow=c(10,8))
for(i in 1:nstation){
    plot(lenseq, data$A[i,]/(data$A[i,]+data$B[i,]), ylim = c(0, 1), ylab = "Prop. of Gear 9", main = paste0("Station ", levels(d$station)[i]))
    lines(lenseq, est.mu[i,])
    lines(lenseq, est.mu[i,] - std.mu[i,], lty = "dashed")
    lines(lenseq, est.mu[i,] + std.mu[i,], lty = "dashed")
}
dev.off()


# estimated rho with observations for each station: 
# calculation of CI: what to do
jpeg(paste(sep = "-","beta-binom-re/NED2013/rho_by_station","species",i.species,"lenbin",b.len,"TMB_SD.jpg"),
     res = 300, width = 24, height = 20, units = "in")
par(mfrow=c(10,8))
for(i in 1:nstation){
    plot(lenseq, data$A[i,]/data$B[i,], ylim = range(0, 5), ylab = "Conversion: 9/15", main = paste0("Station ", levels(d$station)[i]))
    lines(lenseq, est.rho[i,])
    abline(a = 1, b = 0, col = "red")
    lines(lenseq, est.rho[i,] - std.rho[i,], lty = "dashed")
    lines(lenseq, est.rho[i,] + std.rho[i,], lty = "dashed")
}
dev.off()



# ##############################################
# Produce simplified figures for documentation
# ##############################################



