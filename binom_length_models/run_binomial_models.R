# ###############################################
# run a suite of binomial and beta binomial models:
# summarize and compare results
# note: NED2005: same gear (9) different vessel
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)
library(ggplot2)

load("data-NED2005.RData")

# ------------------------------------------------
# organize data for model

i.species <- 11
b.len <- 1

d <- d.length %>% 
    filter(species == i.species) %>%
    transmute(
        station = factor(station),
        vessel = factor(vessel), 
        len = (floor(len/b.len)*b.len+b.len/2), # length grouping
        catch = catch) %>%
    group_by(station, vessel, len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup() %>%
    complete(len, # remove using full_seq()
             station, 
             vessel, 
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

# run TMB model
library(TMB)
version <- "test_betabinom_re"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))


# input for TMB
nlen = length(lenseq)
nstation = nlevels(d$station)
data = list(
    A = d %>% filter(vessel == "NED") %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
    B = d %>% filter(vessel == "TEL") %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
    offset = outer(d.offset$offset,rep(1,length(lenseq))),
    Xf = cs$X %*% eigende$vectors[,1:n_f+n_r],
    Xr = cs$X %*% eigende$vectors[,1:n_r],
    d = eigende$value[1:n_r]
)


# ------------------------------------------------
# run a suite of models
# ------------------------------------------------

# 
parameters1 = list(
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
map1 <- list(
    chol_delta = factor(rep(NA,3))
)

obj1 = MakeADFun(data=data,
                 parameters=parameters1,
                 map = map1,
                 DLL=version,
                 random = c("b", "g", "delta", "epsilon"),
                 silent = F)
opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr)



