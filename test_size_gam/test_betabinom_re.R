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
    transmute(station, gear, len, c) %>%
    spread(gear, c, fill = 0) %>%
    # mutate(len = (floor(len/10))) %>% # length grouping
    group_by(len) %>% 
    summarise(A = (sum(`9`)), B = (sum(`15`)))%>%
    mutate(N = A+B) %>%
    complete(len = full_seq(len, 1), fill = list(A=0,B=0,N=0)) %>%
    print(n=Inf)

# basic functions

library(splines)
x <- splineDesign(knots = d$len, x = 1:100, outer.ok = T)

plot(x[,5])

# input for TMB
nlam = 3
nlen = dim(d)[1]
data = list(n_A = d$A,
            n_B = d$B,
            X = matrix(1, nlen, nlam),
            Xd = matrix(1, nlen, nlam))
parameters = list(beta = rep(1, nlam),
                  betad = rep(1, nlam))


# run TMB model
library(TMB)
version <- "test_betaninom_re"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))
obj = MakeADFun(data=data,
                parameters=parameters,
                DLL=version,
                silent = F)
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)

opt

# plot results
plot(d$len, d$A/d$N, ylim = c(0, 1), ylab = "Prop. of Gear 9")
lines(d$len, obj$report()$mu)


