# #########################################
# Validate/compare my model with Tim's model
# using redfish data from Tim
# #########################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)
library(ggplot2)

load("../../Materials/From_Millar-TMB/redfish.RData")

i.species <- "REDFISH"
b.len <- 1

# read data
d <- redfish.res$spp.dat.aug %>%
    transmute(station = id, len = LENGTH, A = ALBATROSS_RECNUMLEN, B = BIGELOW_RECNUMLEN) %>% 
    gather(vessel, catch, -station, -len) %>%
    filter(catch > 0) %>%
    transmute(
        station = factor(station),
        vessel = factor(vessel), 
        len = (floor(len/b.len)*b.len+b.len/2), # length grouping
        catch = catch) %>%
    group_by(station, vessel, len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup() %>%
    complete(len,
             station, 
             vessel, 
             fill = list(catch = 0)) # complete zero observations

# length at center of each bin, vector for plotting
lenseq <- unique(d$len)


# data for offset: same within a tow
d.offset <- redfish.res$spp.dat.aug %>%
    filter(ALBATROSS_RECNUMLEN + BIGELOW_RECNUMLEN > 0) %>%
    transmute(station = factor(id), offset = offst) %>% 
    group_by(station) %>%
    slice(1) %>%
    ungroup() 


# basis and penalty matrices for cubic spline
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
    A = d %>% filter(vessel == "B") %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
    B = d %>% filter(vessel == "A") %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
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
    # beta = rep(factor(NA),2),
    # gamma = rep(factor(NA),2),
    # log_s_b = factor(NA),
    # log_s_g = factor(NA),
    # log_s_epsilon = factor(NA)
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
jpeg(paste0("beta-binom-re/estimates-CI95_zscore-species_", i.species, ".jpg"),
     res = 300, width = 6, height = 10, units = "in")
par(mfrow=c(3,1))
plot(lenseq, est.mean_mu, ylim = c(0,1), type = "l")
for(i in 1:nstation){lines(lenseq, obj$report()$mu[i,], col = "gray")}
lines(lenseq, est.mean_mu + 1.96*std.mean_mu, col = "blue", lty = "dashed")
lines(lenseq, est.mean_mu - 1.96*std.mean_mu, col = "blue", lty = "dashed")
plot(lenseq, est.mean_phi, type = "l")
plot(lenseq, est.mean_log_rho, ylim = c(-5,5), type = "l")
abline(a = 1, b = 0, col = "red")
for(i in 1:nstation){lines(lenseq, obj$report()$eta_mu[i,], col = "gray")}
lines(lenseq, est.mean_log_rho + std.mean_log_rho, col = "blue", lty = "dashed")
lines(lenseq, est.mean_log_rho - std.mean_log_rho, col = "blue", lty = "dashed")
dev.off()


# compare with Tim's code
# need to run Tim's model before this section

n=as.vector(t(data$A+data$B));
pred <- rep(NA,length(n))
pred[n>0]=betabin.7$report()$mu
mupred <- matrix(pred, dim(data$A)[1], byrow = T)


jpeg(paste0("beta-binom-re/mu_by_station-TMB_SD-species_", i.species, ".jpg"),
     res = 300, width = 20, height = 20, units = "in")
par(mfrow=c(10,10))
for(i in 1:nstation){
    plot(lenseq, data$A[i,]/(data$A[i,]+data$B[i,]), pch =21, cex = 0.2, ylim = c(0, 1), ylab = "Prop. of Gear 9", main = paste0("Station ", levels(d$station)[i]))
    points(lenseq, est.mu[i,], col = "blue", cex = 0.2)
    lines(lenseq, est.mu[i,], col = "blue")
    # lines(lenseq, est.mu[i,] - std.mu[i,], lty = "dashed")
    # lines(lenseq, est.mu[i,] + std.mu[i,], lty = "dashed")
    points(lenseq, mupred[i,], col = "red", cex = 0.2)
    lines(lenseq, mupred[i,], col = "red")
}
dev.off()



# tow-aggregated mu
jpeg(paste0("beta-binom-re/mu_Tim-Yihao-species_", i.species, ".jpg"),
     res = 200, width = 10, height = 8, units = "in")
plot(lenseq, colMeans(est.mu, na.rm = T), ylim = c(0,1), type = "l", col = "blue")
lines(lenseq,  colMeans(mupred, na.rm = T), ylim = c(0,1), type = "l", col = "red")
for(i in 1:nstation){
    points(lenseq+0.004*i-0.2, data$A[i,]/(data$A[i,]+data$B[i,]), pch = 19, cex =0.1)
}
dev.off()


# predicted rho
# compare with Tim's results
jpeg(paste0("beta-binom-re/rho_Tim-Yihao-species_", i.species, ".jpg"),
     res = 200, width = 10, height = 8, units = "in")
plot(lenseq, est.mean_log_rho, ylim = c(-5,5), type = "l", col = "blue")
lines(redfish.res$log_rho_table[,1],(redfish.res$log_rho_table[,2]), col = "red")
for(i in 1:nstation){
    points(lenseq+0.004*i-0.2, log(data$A[i,]/data$B[i,])-data$offset[i,], pch = 19, cex =0.1)
} 
points(lenseq, rep(-5, length(lenseq)), cex = 0.05*colSums(data$A==0))
points(lenseq, rep(5, length(lenseq)), cex = 0.05*colSums(data$B==0))

# add mean and median
x=log(data$A/data$B)-data$offset; x[x==Inf]=NaN; x[x==-Inf]=NaN;
points(lenseq, apply(x, 2, function(x) mean(x, na.rm = T)), col = "green", pch =8)
points(lenseq, apply(x, 2, function(x) median(x, na.rm = T)), col = "yellow", pch = 19)

dev.off()


