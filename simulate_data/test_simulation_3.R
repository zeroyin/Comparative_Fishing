# ###############################################
# Simulate data for model testing
# realistic simulation:
# Catch of A from Maritime survey
# Catch of B = Catch of A * conversion
# 
# Smoothing survey catch by A is neccessary for 
# generating reasonable values for catch by B
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

# -----------------------------------------------
# Catch of A from Maritime survey
# N_a(s,l) <- dens(s,l)*q_A(s,l) 

load("../read_data/MARsurvey.rda")
d.A <- d %>% 
    mutate(station = ID) %>%
    filter(species == 11) %>% 
    group_by(station,len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup %>%
    complete(station, len=full_seq(len,1), fill = list(catch = 0))


# -----------------------------------------------
# True conversion: ratio of two logistic catchability functions

logistic <- function(len,L,k,x0) L/(1+exp(-k*(len-x0)))
# logistic <- function(len,alpha, beta, gamma) gamma * exp(alpha+beta*len)/(1+exp(alpha+beta*len))
len <- full_seq(d.A$len,1)
q_A <- logistic(len, 1, 0.2, 30)
q_B <- logistic(len, 1, 0.6, 40)
matplot(len, cbind(q_A,q_B))

rho <- q_B/q_A
plot(len, rho, type = "l", col = "blue")
abline(h = 1)


# -----------------------------------------------
# Smooth catch by A:
# (1) fit parametric curve by station
# (2) smoothing
# 
# Simulate catch of B using conversion
# N_B(s,l) <- dens(s,l)*q_A(s,l)*rho(s,l)*error

d.B <- d.A %>%
    group_by(station) %>%
    mutate(catch = rollmean(catch, 10, fill = 0)*rho*rlnorm(length(len), 0, 1)) %>%
    ungroup()


bind_rows(d.A %>% mutate(gear="A"), d.B %>% mutate(gear="B")) %>%
    filter(catch > 0) %>%
    ggplot() +
    geom_tile(aes(as.factor(station), len, fill = catch)) +
    scale_fill_continuous(low = "white", high = "red", trans = "log10", na.value = "white", limits = c(1, NA)) +
    theme(axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          panel.background = element_blank(),
          panel.grid= element_blank()) +
    facet_wrap(~gear, ncol = 1)


N_A <- d.A %>% spread(len, catch) %>% select(-station) %>% as.matrix()
N_A <- d.B %>% spread(len, catch) %>% select(-station) %>% as.matrix()














