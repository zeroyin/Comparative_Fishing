# ###############################################
# Simulate data for model testing
# realistic simulation:
# Catch of A from Maritime survey
# Catch of B = Catch of A * conversion
# 
# Smoothing survey catch by A is neccessary for 
# generating reasonable values for catch by B
# 
# Resample within homogeneous stratum
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

load("../read_data/MARsurvey.rda")
load("../read_data/MARgrid_1km.rda")
# dens.base <- d %>% 
#     filter(species == 11) %>%
#     group_by(len) %>% 
#     summarise(catch = mean(catch)) %>%
#     ungroup() %>%
#     complete(len=full_seq(len,1), fill = list(catch = NA)) %>%
#     na.fill("extend") %>%
#     as.data.frame()

# variance for each stratum
d %>% 
    filter(species == 11) %>% 
    select(ID, len, catch) %>%
    complete(ID, len=full_seq(len,1), fill = list(catch = 0)) %>%
    group_by(ID) %>%
    mutate(A = na.fill(rollmean(catch, 5, fill = NA), "extend")) %>%
    ungroup %>%
    inner_join(d %>% filter(species == 11) %>% select(ID, strat), by = "ID") %>%
    group_by(strat, len) %>%
    mutate(A = mean(A)) %>%
    ungroup %>% 
    mutate(diff_catch = (catch - A)/A) %>%
    group_by(strat) %>%
    summarise(catch_sd = sd(diff_catch, na.rm = T)) %>%
    ungroup() %>%
    mutate(strat = as.factor(strat)) %>%
    inner_join(MARgrid, by = "strat") %>%
    ggplot() +
    geom_point(aes(x = lon, y = lat, color = catch_sd)) +
    scale_color_continuous(high = "red", low = "white") +
    coord_quickmap(xlim = c(-69,-57), ylim = c(41, 47)) +
    theme_bw()+
    theme(legend.position="bottom", axis.title = element_blank())



# -----------------------------------------------
# Moving average density*q_A from Maritime survey
# Assuming homogeneity within stratum: summarize by strata

d.dens <- d %>% 
    filter(species == 11) %>% 
    select(ID, len, catch) %>%
    complete(ID, len=full_seq(len,1), fill = list(catch = 0)) %>%
    group_by(ID) %>%
    mutate(A = ceiling(na.fill(rollmean(catch, 5, fill = NA), "extend")*1000)/1000) %>%
    ungroup %>%
    inner_join(d %>% filter(species == 11) %>% select(ID, strat), by = "ID") %>%
    group_by(strat, len) %>%
    summarize(A = mean(A)) %>%
    ungroup
    

# -----------------------------------------------
# True conversion: ratio of two logistic catchability functions

logistic <- function(len,L,k,x0) L/(1+exp(-k*(len-x0)))
# logistic <- function(len,alpha, beta, gamma) gamma * exp(alpha+beta*len)/(1+exp(alpha+beta*len))
len <- full_seq(d.dens$len,1)
q_A <- logistic(len, 1, 0.1, 30)
q_B <- logistic(len, 1, 0.2, 20)


jpeg(filename = "rho.jpg",width = 600,height = 800)
par(mfrow=c(2,1))
matplot(len, cbind(q_A,q_B), type = "l", col = "blue")
legend("bottomright", lty = c(1,2), col="blue", legend = c("q_A", "q_B"))   

rho <- q_A/q_B
plot(len, log(rho), type = "l", col = "blue")
abline(h = 0)
dev.off()

# -----------------------------------------------
# density*q_B = density*q_A*rho
# spatial variation for rho?
d.dens <- d.dens %>%
    group_by(strat) %>%
    mutate(B = A * rho) %>%
    ungroup


d.dens %>%
    gather(type, value, -strat, -len) %>%
    filter(value > 0) %>%
    ggplot() +
    geom_tile(aes(as.factor(strat), len, fill = value)) +
    scale_fill_continuous(low = "white", high = "red", trans = "log10", na.value = "white", limits = c(1, NA)) +
    theme(axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          panel.background = element_blank(),
          panel.grid= element_blank()) +
    facet_wrap(~type, ncol = 1)

    
    
# -----------------------------------------------
# Simulate catch of A from resampling of smoothed density*q_A
# Simulate catch of B similarly:
# 
# N_a(s,l) <- dens(s,l)*q_A(s,l)*error
# N_B(s,l) <- dens(s,l)*q_A(s,l)*rho(s,l)*error


# one pair for each stratum
d.catch <- d.dens %>%
    mutate(station = strat) %>%
    group_by(station) %>%
    mutate(c.A = rpois(length(len), A*rlnorm(length(len),-1^2/2,1)),
           c.B = rpois(length(len), B*rlnorm(length(len),-1^2/2,1))) %>%
    ungroup()


d.catch %>%
    select(len, station, c.A, c.B, A, B) %>%
    gather(type, value, -station, -len) %>%
    filter(value > 0) %>%
    ggplot() +
    geom_tile(aes(as.factor(station), len, fill = value)) +
    scale_fill_continuous(low = "white", high = "red", trans = "log10", na.value = "white", limits = c(1, NA)) +
    theme(axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          panel.background = element_blank(),
          panel.grid= element_blank()) +
    facet_wrap(~type, ncol = 2) 
ggsave("dens-catch.jpg", width = 10, height = 6)


# -----------------------------------------------
# Fit model

dir.create("res/output/", recursive = T)

# load TMB model
library(TMB)
dyn.load(dynlib("binom"))
# load function for model fitting
source("fit_model.R")


N_A <- d.catch %>% select(len, station, c.A) %>% spread(len, c.A) %>% select(-station) %>% as.matrix()
N_B <- d.catch %>% select(len, station, c.B) %>% spread(len, c.B) %>% select(-station) %>% as.matrix()

simu <- list(N_A=N_A,N_B=N_B,len_list=len)

model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4),paste0("ZB", 2:3))
for(i.model in model_vec){
    res <- try(fit_model(i.model, simu))
    if(!is.null(res)){if(!inherits(res, "try-error")){
        save(res, file = paste0("res/output/res-",i.model,".rda"))
    }}
}




# AIC table
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4),paste0("ZB", 2:3))
aic_mat <- rep(NA, length(model_vec))
names(aic_mat) <- model_vec
for(i_model in 1:length(model_vec)){
    res_file <- paste0("res/output/res-",model_vec[i_model],".rda")
    if(file.exists(res_file)){
        load(res_file)
        aic_mat[i_model] <- res$aic
        rm("res", "res_file")
    }
}

t(round(aic_mat, digits = 0)) %>%
    write.csv(file = "res/aic_table.csv")



jpeg(paste0("res/est_mu-",".jpg"),res = 600, width = 10, height = 8, units = "in")
matplot(simu$len_list, t(simu$N_A/(simu$N_A + simu$N_B)), type = "p", cex = 0.2, col = "black", xlab = "len", ylab = NA)
lines(simu$len_list, colMeans(simu$N_A/(simu$N_A + simu$N_B), na.rm = T), col = "orange")
lines(simu$len_list, boot::inv.logit(log(1/rho)), col = "black")
for(i.model in model_vec){
    res_file <- paste0("res/output/res-",i.model,".rda")
    if(file.exists(res_file)){
        load(res_file)
        est.mean_mu <- res$obj$report()$mean_mu
        lines(simu$len_list, est.mean_mu, type = "l",  col = "blue")
        text(x=simu$len_list[1]-1, y=est.mean_mu[1], labels = i.model,
             col = ifelse(i.model == model_vec[which.min(aic_mat)], "red", "blue"))
        rm("res", "res_file", "est.mean_mu")
    }
}
legend("bottomright", pch = c(1,NA,NA), lty = c(0,1,1), 
       col = c("black","orange","black"),
       legend = c("observed mu", "sample mean of mu", "true mu"))    
dev.off()

