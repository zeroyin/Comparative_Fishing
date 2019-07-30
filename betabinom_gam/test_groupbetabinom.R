# #########################################
# test beta binomial model with subgrouping
# #########################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/betabinom_gam/")

library(dplyr)
library(tidyr)
library(ggplot2)

load("../read_data/data-NED2005.RData")

# organize data for model

i.species <- 204
b.len <- 1

d <- d.length %>% 
    filter(species == i.species) %>%
    transmute(
        station = station,
        vessel = vessel, 
        len = (floor(len/b.len)*b.len+b.len/2), # length grouping
        catch = catch
    ) %>%
    group_by(station, vessel, len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup() %>%
    complete(
        len, # remove using full_seq()
        station, 
        vessel,  
        fill = list(catch = 0),
    ) %>% # complete zero observations
    inner_join( # add grouping index
        d.length %>%
            mutate(group = strat) %>%
            distinct(station, group),
        by = "station") %>%
    mutate(station = factor(station),
           vessel = factor(vessel),
           group = as.integer(factor(group)))


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

# run TMB model
library(TMB)
version <- "test_groupbetabinom"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))


# input for TMB
nlen = length(lenseq)
nstation = nlevels(d$station)
nind = n_distinct(d$group)
data = list(
    A = d %>% filter(vessel == "NED") %>% spread(len, catch) %>% select(-station, -vessel, -group) %>% as.matrix(),
    B = d %>% filter(vessel == "TEL") %>% spread(len, catch) %>% select(-station, -vessel, -group) %>% as.matrix(),
    offset = outer(d.offset$offset,rep(1,length(lenseq))),
    Xf = cs$X %*% eigende$vectors[,1:n_f+n_r],
    Xr = cs$X %*% eigende$vectors[,1:n_r],
    d = eigende$value[1:n_r],
    ind = distinct(d, station, group)$group,
    n_ind = nind
)
parameters = list(
    beta = rep(0, n_f),
    b = rep(0, n_r),
    gamma = rep(0, n_f),
    g = rep(0, n_r),
    delta = matrix(0, nind, n_f),
    epsilon = matrix(0, nind, n_r),
    log_s_b = log(10),
    log_s_g = log(10),
    log_s_epsilon = log(10),
    chol_delta = c(1,0,1) # use chol decomp in vector form
)
map <- list()
obj = MakeADFun(data=data,
                parameters=parameters,
                map=map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
opt <- nlminb(obj$par,obj$fn,obj$gr)

i.model <- "GB"
# check convergence, maximum gradient and positive definite
if(exists("opt")){
    if(!opt$convergence){
        gra <- obj$gr()
        hes <- eigen(optimHess(par=opt$par, fn=obj$fn, gr=obj$gr))$values
        if(max(abs(gra)) < 0.1 & min(hes) > -0.1){
            aic <- 2*sum(obj$report()$nll) + 2*10
            rep <- try(sdreport(obj))
            res <- list(obj = obj, opt = opt, aic = aic, gra = gra, hes = hes)
            save(res, file = paste0("beta-binom-gr/NED2005/", i.species, "-",i.model,".rda"))
        }
    }
}


# #################################################
# Results: report and visualize

# rep <- sdreport(obj)


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
jpeg(paste(sep = "-","beta-binom-gr/NED2005/estimates","species",i.species,"lenbin",b.len,"CI95_zscore.jpg"),
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
# SD is from TMB sd
jpeg(paste(sep = "-","beta-binom-gr/NED2005/mu_by_station","species",i.species,"lenbin",b.len,"TMB_SD.jpg"),
     res = 300, width = 24, height = 20, units = "in")
par(mfrow=c(10,10))
for(i in 1:nstation){
    plot(lenseq, data$A[i,]/(data$A[i,]+data$B[i,]), ylim = c(0, 1), ylab = "Prop. of Gear 9", main = paste0("Station ", levels(d$station)[i]))
    lines(lenseq, est.mu[i,])
    lines(lenseq, est.mu[i,] - std.mu[i,], lty = "dashed")
    lines(lenseq, est.mu[i,] + std.mu[i,], lty = "dashed")
}
dev.off()


# estimated rho with observations for each station: 
# calculation of CI: what to do
jpeg(paste(sep = "-","beta-binom-gr/NED2005/rho_by_station","species",i.species,"lenbin",b.len,"TMB_SD.jpg"),
     res = 300, width = 24, height = 20, units = "in")
par(mfrow=c(10,10))
for(i in 1:nstation){
    plot(lenseq, data$A[i,]/data$B[i,], ylim = range(0, 5), ylab = "Conversion: 9/15", main = paste0("Station ", levels(d$station)[i]))
    lines(lenseq, est.rho[i,])
    abline(a = 1, b = 0, col = "red")
    lines(lenseq, est.rho[i,] - std.rho[i,], lty = "dashed")
    lines(lenseq, est.rho[i,] + std.rho[i,], lty = "dashed")
}
dev.off()



