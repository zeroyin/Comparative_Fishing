# ###############################################
# Simulate data for model testing
# Run a set of simulations, assess performance
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)


# ###############################################

# function to simulate data
simu_data <- function(data_survey, i.species = 11, nrep_stn = 1){
    
    # -----------------------------------------------
    # Moving average density*q_A from Maritime survey
    # Assuming homogeneity within stratum: summarize by strata
    
    d.dens <- d %>% 
        filter(species == i.species) %>% 
        select(ID, len, catch) %>%
        complete(ID, len=full_seq(len,1), fill = list(catch = 0)) %>%
        group_by(ID) %>%
        mutate(A = ceiling(na.fill(rollmean(catch, 5, fill = NA), "extend")*1000)/1000) %>%
        ungroup %>%
        inner_join(d %>% filter(species == i.species) %>% select(ID, strat), by = "ID") %>%
        group_by(strat, len) %>%
        summarize(A = mean(A)) %>%
        ungroup
    
    
    # -----------------------------------------------
    # True conversion: ratio of two logistic catchability functions
    
    len <- full_seq(d.dens$len,1)
    
    # logistic <- function(len,L,k,x0) L/(1+exp(-k*(len-x0)))
    # # logistic <- function(len,alpha, beta, gamma) gamma * exp(alpha+beta*len)/(1+exp(alpha+beta*len))
    # q_A <- logistic(len, 1, 0.1, 30)
    # q_B <- logistic(len, 1, 0.2, 20)
    # 
    # mu <- q_B/(q_A+q_B)
    # rho <- q_A/q_B
    
    mu <- 1-1/(exp(0.04*len-0.5)+1)
    rho <- exp(boot::logit(mu))
    
    
    # -----------------------------------------------
    # density*q_B = density*q_A*rho
    # spatial variation for rho?
    d.dens <- d.dens %>%
        group_by(strat) %>%
        mutate(B = A * rho) %>%
        ungroup
    
    
    # -----------------------------------------------
    # Simulate catch of A from resampling of smoothed density*q_A
    # Simulate catch of B similarly:
    # 
    # N_a(s,l) <- dens(s,l)*q_A(s,l)*error
    # N_B(s,l) <- dens(s,l)*q_A(s,l)*rho(s,l)*error
    
    
    # one pair for each stratum
    use.seq <- seq(0:(max(len)+3))
    use.sd <- dbeta((use.seq-min(use.seq))/(max(use.seq)-min(use.seq)),0.2,0.2)[len]
    d.catch1 <- d.dens %>%
        mutate(station = paste0(strat,1)) %>%
        group_by(station) %>%
        mutate(c.A = rpois(length(len), A*rlnorm(length(len),-use.sd^2/2,use.sd)),
               c.B = rpois(length(len), B*rlnorm(length(len),-use.sd^2/2,use.sd))) %>%
        ungroup()
    d.catch2 <- d.dens %>%
        mutate(station = paste0(strat,2)) %>%
        group_by(station) %>%
        mutate(c.A = rpois(length(len), A*rlnorm(length(len),-use.sd^2/2,use.sd)),
               c.B = rpois(length(len), B*rlnorm(length(len),-use.sd^2/2,use.sd))) %>%
        ungroup()
    d.catch <- bind_rows(d.catch1, d.catch2)
    
    # remove length gaps
    N_A <- d.catch %>% select(len, station, c.A) %>% spread(len, c.A) %>% select(-station) %>% as.matrix()
    N_B <- d.catch %>% select(len, station, c.B) %>% spread(len, c.B) %>% select(-station) %>% as.matrix()
    
    simu <- list(N_A=N_A,N_B=N_B,len_list=len, mu = mu, rho = rho)
    
    return(simu)
    
}

# function to fit model, given model name
source("fit_model.R")


# ##################################################


# load TMB model
library(TMB)
version <- "binom"
# compile(paste0(version,".cpp"))
dyn.load(dynlib(version))

dir.create("res2/output/",recursive = T)
dir.create("res2/est_mu/",recursive = T)

# Read survey data
load("../read_data/MARsurvey.rda")
load("../read_data/MARgrid_1km.rda")

# run repeated simulations
n_seed <- 1
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
for(i.seed in 1:n_seed){
    set.seed(i.seed)
    simu <- simu_data(data_survey = d, i.species = 11)
    save(simu, file = paste0("res2/output/data-",i.seed,".rda"))
    for(i.model in model_vec){
        res <- try(fit_model(i.model, simu))
        if(!is.null(res)){if(!inherits(res, "try-error")){
            save(res, file = paste0("res2/output/res-",i.seed,"-",i.model,".rda"))
        }}
    }
}


# organize results: AIC table

aic_mat <- matrix(NA, n_seed, length(model_vec))
colnames(aic_mat) <- model_vec
for(i.seed in 1:n_seed){
    for(i.model in model_vec){
        res_file <- paste0("res2/output/res-",i.seed,"-",i.model,".rda")
        if(file.exists(res_file)){
            load(res_file)
            aic_mat[i.seed, i.model] <- res$aic
            rm("res", "res_file")
        }
    }
}

t(round(aic_mat, digits = 0))

round(aic_mat - apply(aic_mat,1,function(x){min(x,na.rm = T)}), digits = 0) %>%
    write.csv(file = "res2/aic_table.csv")


aic.summary <- data.frame(
    nconv = as.integer(apply(aic_mat, 2, function(x){sum(!is.na(x))})),
    nbest = as.integer(colSums(aic_mat==apply(aic_mat,1,function(x){min(x,na.rm = T)}), na.rm = T)),
    row.names = model_vec)

xtable::xtable(t(aic.summary))



# organize results: estimated mu

for(i.seed in 1:n_seed){
    
    load(paste0("res2/output/data-",i.seed, ".rda"))
    
    jpeg(paste0("res2/est_mu/est_mu-",i.seed,".jpg"),res = 600, width = 10, height = 8, units = "in")
    matplot(simu$len_list, t(simu$N_A/(simu$N_A + simu$N_B)), type = "p", cex = 0.2, col = "black", xlab = "len", ylab = NA)
    lines(simu$len_list, colMeans(simu$N_A/(simu$N_A + simu$N_B), na.rm = T), col = "orange")
    lines(simu$len_list, simu$mu, col = "black")
    for(i.model in model_vec){
        res_file <- paste0("res2/output/res-",i.seed,"-",i.model,".rda")
        if(file.exists(res_file)){
            load(res_file)
            est.mean_mu <- res$obj$report()$mean_mu
            col = ifelse(i.model == model_vec[which.min(aic_mat[i.seed,])], "red", "blue")
            lines(simu$len_list, est.mean_mu, type = "l",  col = col)
            text(x=simu$len_list[1]-1, y=est.mean_mu[1], labels = i.model, col = col)
            rm("res", "res_file", "est.mean_mu", "col")
        }
    }
    legend("bottomright", pch = c(1,NA,NA), lty = c(0,1,1), 
           col = c("black","orange","black"),
           legend = c("obs mu/stn", "obs mu/mean", "true mu"))    
    dev.off()
    
    rm("simu")
}




# organize results:
# Metric to assess variance structure

load("res2/output/data-1.rda")
len <- simu$len_list
mu <- simu$mu

# fill rate/observation rate
n_enc_mat <- matrix(NA, n_seed, length(len))
sd_stn_mat <- matrix(NA, n_seed, length(len))
for(i.seed in 1:n_seed){
    dat_file <- paste0("res2/output/data-",i.seed, ".rda")
    load(dat_file)
    n_enc_mat[i.seed,] <- colMeans((simu$N_A+simu$N_B)>0, na.rm = T)
    sd_stn_mat[i.seed,] <- apply(simu$N_B/(simu$N_A+simu$N_B)-simu$mu,2,function(x)sd(x, na.rm = T))
    rm("simu")
}



# organize results: residuals

# AIC
aic_mat <- read.csv("res2/aic_table.csv")[,-1]

summary_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4),"AIC")
wmresid_vec <- rep(NA, length(summary_vec)); names(wmresid_vec) <- summary_vec
for(i.model in summary_vec){
    
    resid_mat <- matrix(NA, n_seed, length(len))
    
    # compute residuals
    for(i.seed in 1:n_seed){
        
        if(i.model == "AIC"){
            # AIC selected best model
            model <- model_vec[which.min(aic_mat[i.seed,])]
        }else{
            model <- i.model
        }
        
        dat_file <- paste0("res2/output/data-",i.seed, ".rda")
        res_file <- paste0("res2/output/res-",i.seed,"-",model,".rda")
        load(dat_file)
        if(file.exists(res_file)) {
            load(res_file)
            resid_mat[i.seed,] <- res$obj$report()$mean_mu - simu$mu
        }
        rm("simu","res", "mean_mu")
    }
    
    # residual plots
    resid_mat %>%
        as.data.frame() %>%
        `colnames<-`(len)  %>%
        gather() %>%
        `colnames<-`(c("len","resid")) %>%
        mutate(len = as.integer(len)) %>%
        ggplot(aes(x = len, y = resid)) +
        geom_point( color = "blue", size = 0.2, alpha = 0.2) +
        geom_hline(yintercept = 0, color = "black") +
        stat_summary(fun.y = median, geom="point", color = "blue") +
        theme_bw()
    ggsave(filename = paste0("res2/resid_median-",i.model,"-100simu.jpg"),width = 10, height = 6)
    
    
    # mean standard deviation weighed by encounter rate 
    wmresid_vec[i.model] <- mean(abs(resid_mat) * n_enc_mat , na.rm = T)
    
    rm("resid_mat")
    
}



c(
    wmresid_vec,
    "enc_rate" = mean(n_enc_mat, na.rm = T), # ecnounter rate
    "stn_sd" = mean(sd_stn_mat, na.rm = T), # station-variation rate
    "flat_mu" = sd(mu)/mean(mu) # flatness of mu
) %>% 
    write.csv(file = "res2/wmresid.csv")




