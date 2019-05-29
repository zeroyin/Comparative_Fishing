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

load("data-NED2016.RData")

# organize data for model

i.species <- 11

d <- d.length %>% 
    filter(species == i.species) %>%
    transmute(
        station = factor(station),
        gear = factor(gear), 
        len = (floor(len/5)), # length grouping
        catch = c) %>%
    group_by(station, gear, len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup() %>%
    complete(len = full_seq(len, 1), station, gear, fill = list(catch = 0))


# basis abd penalty matrices for cubic spline
# over length bins
library(mgcv)
cs <- smooth.construct(
    object = s(len, bs = "cr"), 
    data = d %>% group_by(len) %>% summarise(catch = sum()), 
    knots = NULL)

# default to 10 knots resulting in 2 fixed and 8 random effects
n_f <- 2
n_r <- 8
eigende <- eigen(cs$S[[1]])

# input for TMB
nlen = nlevels(as.factor(d$len))
nstation = nlevels(d$station)

data = list(
    A = d %>% filter(gear == 9 ) %>% spread(len, catch) %>% select(-station, -gear) %>% as.matrix(),
    B = d %>% filter(gear == 15) %>% spread(len, catch) %>% select(-station, -gear) %>% as.matrix(),
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
    log_s_b = log(1),
    log_s_g = log(1),
    log_s_epsilon = log(1),
    C_delta = diag(1,2)
)
map <- list(
    log_s_b = factor(NA),
    log_s_g = factor(NA),
    log_s_epsilon = factor(NA)
)

# run TMB model
library(TMB)
version <- "test_betaninom_re"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))
obj = MakeADFun(data=data,
                parameters=parameters,
                # map = map,
                DLL=version,
                random = c("delta", "epsilon"),
                silent = F)
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)

opt

# plot results
plot(d$len, d$A/d$N, ylim = c(0, 1), ylab = "Prop. of Gear 9")
lines(d$len, obj$report()$mu)


