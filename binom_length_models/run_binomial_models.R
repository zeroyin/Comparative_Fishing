# ###############################################
# run a suite of binomial and beta binomial models:
# summarize and compare results
# note: NED2005: same gear (9) different vessel
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/binom_length_models/")

library(dplyr)
library(tidyr)
library(ggplot2)

load("../read_data/data-NED2005.RData")

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
version <- "test_betabinom_re_v2"
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

# ------------------------------------------------
# run a suite of models
# ------------------------------------------------

# function to run model given mapped "obj"
run_model <- function(model, obj){
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    if(exists("opt")){
        if(!opt$convergence){
            rep <- sdreport(obj)
            res <- list(obj = obj, opt = opt, rep = rep)
            save(res, file = paste0("res-",model,".rda"))
            
            # estimate and std of mu and phi and rho
            est <- summary(rep, "report")[,"Estimate"]
            std <-  summary(rep, "report")[,"Std. Error"]
            est.mean_mu <- est[names(est) == "mean_mu"]
            std.mean_mu <- std[names(std) == "mean_mu"]
            est.mean_phi <- est[names(est) == "mean_phi"]
            std.mean_phi <- std[names(std) == "mean_phi"]
            est.mean_log_rho <- est[names(est) == "mean_log_rho"]
            std.mean_log_rho <- std[names(std) == "mean_log_rho"]
            
            jpeg(paste(sep = "-",model,"estimates","species",i.species,"lenbin",b.len,"CI95_zscore.jpg"),
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
        }
    }
}



# BB2
rm("map","obj","opt","rep")
map = list(
    delta = factor(matrix(NA, nstation, n_f)),
    epsilon = factor(matrix(NA, nstation, n_r)),
    g = factor(rep(NA,n_r)),
    log_s_g = factor(NA),
    log_s_epsilon = factor(NA),
    chol_delta = factor(rep(NA,3))
)
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
run_model("BB2", obj)



# BB3
rm("map","obj","opt","rep")
map = list(
    delta = factor(matrix(NA, nstation, n_f)),
    epsilon = factor(matrix(NA, nstation, n_r)),
    log_s_epsilon = factor(NA),
    chol_delta = factor(rep(NA,3))
)
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
run_model("BB3", obj)




# BB4
rm("map","obj","opt","rep")
map = list(
    epsilon = factor(matrix(NA, nstation, n_r)),
    g = factor(rep(NA,n_r)),
    log_s_g = factor(NA),
    log_s_epsilon = factor(NA),
)
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
run_model("BB4", obj)




# BB5
rm("map","obj","opt","rep")
map = list(
    epsilon = factor(matrix(NA, nstation, n_r)),
    log_s_epsilon = factor(NA)
)
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
run_model("BB5", obj)



# BB6
rm("map","obj","opt","rep")
map <- list(
    g = factor(rep(NA,n_r)),
    log_s_g = factor(NA)
)
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
run_model("BB6", obj)



# BB7
rm("map","obj","opt","rep")
map <- list()
obj = MakeADFun(data=data,
                parameters=parameters,
                map = map,
                DLL=version,
                random = c("b", "g", "delta", "epsilon"),
                silent = F)
run_model("BB7", obj)

