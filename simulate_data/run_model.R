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
library(TMB)

# ###############################################

# function to simulate data
simu_data <- function(d, i.species, mu_method=2, phi_method = 1,  n_rep_stn){
    
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
    
    len <- full_seq(d.dens$len,1)
    strat <- unique(d.dens$strat)
    
    nlen <- length(len)
    nstrat <- length(strat)
    
    # number of stations
    nstn <- n_rep_stn * nstrat
    
    
    # -----------------------------------------------
    # True conversion: 
    
    if (mu_method == 1){
        # ratio of two logistic catchability functions
        logistic <- function(len,L,k,x0) L/(1+exp(-k*(len-x0)))
        q_A <- logistic(len, 1, 0.2, len[floor(length(len)*1/3)])
        q_B <- logistic(len, 1, 0.1, len[floor(length(len)*1/2)])
        mu <-  matrix(rep(q_B/(q_A+q_B),nstn),nstn,nlen, byrow = T)
        rho <- exp(boot::logit(mu))
        
    }else if (mu_method == 2){
        # exponential or linear: a + b*len; a+b*exp(-c*len)
        mu <- matrix(rep(1 - 1/(exp(0.7653+2.4198*exp(-0.2607*len))+1),nstn),nstn,nlen, byrow = T)
        rho <- exp(boot::logit(mu))
        
    }else if (mu_method == 3){
        # station-level variation of shape
        mu <- rho <- matrix(NA, nstn, nlen)
        dimnames(mu)[[2]] <- dimnames(rho)[[2]] <- len
        for(i.stn in 1:nstn){
            mu[i.stn,] <- 1 - 1/(exp(0.03*rlnorm(1,0-0.3^2/2,0.3)*len-2*rlnorm(1,0-0.3^2/2,0.3))+1)
            rho[i.stn,] <- exp(boot::logit(mu[i.stn,]))
        }
    }
   

    
    # heteroschedecity
    if (phi_method == 1){
        phi <- rep(0.4, nlen)  
        
    }else if (phi_method == 2){
        phi_len <- seq(0:(max(len)+3))
        phi <- dbeta((phi_len-min(phi_len))/(max(phi_len)-min(phi_len)),0.5,0.5)[len]
        
    }
    
    # -----------------------------------------------
    # density*q_B = density*q_A*rho
    # 
    # Simulate catch of A from resampling of smoothed density*q_A
    # Simulate catch of B similarly:
    # 
    # N_a(s,l) <- dens(s,l)*q_A(s,l)*error
    # N_B(s,l) <- dens(s,l)*q_A(s,l)*rho(s,l)*error
    
    
    list.catch <- list()
    for(i.strat in 1:nstrat){
        for(i.rep in 1:n_rep_stn){
            i.stn <- i.rep+(i.strat-1)*n_rep_stn
            list.catch[[i.stn]] <- d.dens[d.dens$strat == strat[i.strat], ] %>%
                mutate(station = paste0(strat[i.strat],i.rep)) %>%
                mutate(B = A / rho[i.stn,]) %>%
                mutate(c.A = rpois(length(len), A*rlnorm(length(len),-phi^2/2,phi)),
                       c.B = rpois(length(len), B*rlnorm(length(len),-phi^2/2,phi)))
        }
    }
    
    d.catch <- bind_rows(list.catch) %>% mutate(len = as.integer(len))
    
    
    # remove length gaps
    N_A <- d.catch %>% select(len, station, c.A) %>% spread(len, c.A) %>% select(-station) %>% as.matrix()
    N_B <- d.catch %>% select(len, station, c.B) %>% spread(len, c.B) %>% select(-station) %>% as.matrix()
    
    simu <- list(N_A=N_A,N_B=N_B,len_list=len, mu = mu, rho = rho)
    
    return(simu)
    
}

# function to fit model, given model name
source("fit_model.R")

# summarize results
write_results <- function(n_seed, save.i.scenario, save.dir.output, save.dir.est_mu){
    
    # run repeated simulations
    model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
    
    # organize results: AIC table
    aic_mat <- matrix(NA, n_seed, length(model_vec))
    colnames(aic_mat) <- model_vec
    for(i.seed in 1:n_seed){
        for(i.model in model_vec){
            res_file <- paste0(save.dir.output, "res-",i.seed,"-",i.model,".rda")
            if(file.exists(res_file)){
                load(res_file)
                aic_mat[i.seed, i.model] <- res$aic
                rm("res", "res_file")
            }
        }
    }
    
    round(aic_mat - apply(aic_mat,1,function(x){min(x,na.rm = T)}), digits = 0) %>%
        write.csv(file = paste0(save.i.scenario, "aic_table.csv"))
    
    data.frame(
        nconv = as.integer(apply(aic_mat, 2, function(x){sum(!is.na(x))})),
        nbest = as.integer(colSums(aic_mat==apply(aic_mat,1,function(x){min(x,na.rm = T)}), na.rm = T)),
        row.names = model_vec) %>%
        write.csv(file = paste0(save.i.scenario, "aic_summary.csv"))
    
    
    
    # organize results: estimated mu
    
    for(i.seed in 1:n_seed){
        
        load(paste0(save.dir.output, "/data-",i.seed, ".rda"))
        
        jpeg(paste0(save.dir.est_mu, "est_mu-",i.seed,".jpg"),res = 600, width = 10, height = 8, units = "in")
        matplot(simu$len_list, t(simu$N_A/(simu$N_A + simu$N_B)), type = "p", cex = 0.05, pch =20, col = "black", xlab = "len", ylab = NA)
        matplot(simu$len_list, t(simu$mu), col = "grey", type = "l", lty = "dashed", add = T)
        lines(simu$len_list, colMeans(simu$N_A/(simu$N_A + simu$N_B), na.rm = T), col = "orange")
        for(i.model in model_vec){
            res_file <- paste0(save.dir.output, "res-",i.seed,"-",i.model,".rda")
            if(file.exists(res_file)){
                load(res_file)
                est.mean_mu <- res$obj$report()$mean_mu
                col = ifelse(i.model == model_vec[which.min(aic_mat[i.seed,])], "red", "blue")
                lines(simu$len_list, est.mean_mu, type = "l",  col = col)
                text(x=simu$len_list[1]-1, y=est.mean_mu[1], labels = i.model, col = col)
                rm("res", "res_file", "est.mean_mu", "col")
            }
        }
        legend("bottom", pch = c(1,NA,NA), lty = c(0,1,1), 
               col = c("black","orange","black"),
               legend = c("obs mu/stn", "obs mu/mean", "true mu"))    
        dev.off()
        
        rm("simu")
    }
    
    
    # organize results:
    # Metric to assess variance structure
    
    load(paste0(save.dir.output, "data-1.rda"))
    len <- simu$len_list
    
    # fill rate/observation rate
    p_enc_mat <- matrix(NA, n_seed, length(len))
    mu_lin <- sd_len <- sd_stn <- p_enc_A <- p_enc_B <- rep(NA, n_seed)
    for(i.seed in 1:n_seed){
        dat_file <- paste0(save.dir.output, "/data-",i.seed, ".rda")
        load(dat_file)
        p_enc_mat[i.seed,] <- colMeans((simu$N_A+simu$N_B)>0)
        p_enc_A[i.seed] <- mean((simu$N_A)>0, na.rm = T)
        p_enc_B[i.seed] <- mean((simu$N_B)>0, na.rm = T)
        sd_stn[i.seed] <- mean(apply(as.data.frame(simu$N_B/(simu$N_A+simu$N_B))-simu$mu,2,function(x) sd(x, na.rm=T)),na.rm=T)
        sd_len[i.seed] <- mean(apply(simu$N_B/(simu$N_A+simu$N_B),1,function(x) sd(x, na.rm=T)),na.rm=T)
        mu_lin[i.seed] <- mean(apply(simu$mu, 1, function(x) deviance(lm(x~simu$len_list))))
        rm("simu")
    }
    
    
    
    # organize results: residuals
    
    # AIC
    aic_mat <- read.csv(paste0(save.i.scenario, "aic_table.csv"))[,-1]
    
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
            
            dat_file <- paste0(save.dir.output, "/data-",i.seed, ".rda")
            res_file <- paste0(save.dir.output, "res-",i.seed,"-",model,".rda")
            load(dat_file)
            if(file.exists(res_file)) {
                load(res_file)
                resid_mat[i.seed,] <- colMeans(as.data.frame(res$obj$report()$mu) - simu$mu)
            }
            rm("simu","res")
        }
        
        
        if(!all(is.na(resid_mat))){
            
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
                ylim(-0.3,0.3) +
                theme_bw()
            ggsave(filename = paste0(save.i.scenario,"resid_median-",i.model,"-100simu.jpg"),
                   width = 10, height = 6)
            
        }
        
        
        # mean standard deviation weighed by encounter rate 
        wmresid_vec[i.model] <- mean(abs(resid_mat)*p_enc_mat/rowMeans(p_enc_mat), na.rm = T)
        
        rm("resid_mat")
        
    }
    
    
    c(
        wmresid_vec,
        "enc_prop" = mean(p_enc_mat, na.rm = T), # ecnounter rate
        "stn_sd" = mean(sd_stn, na.rm = T), # station-variation rate
        "len_sd" = mean(sd_len, na.rm = T), # length conversion-variation rate
        "mu_lin" = mean(mu_lin, na.rm = T) # rss of true mu - lm(mu~len)
    ) %>% 
        write.csv(file = paste0(save.i.scenario, "wmresid.csv"))
    
}

# function to run model
run_model <- function(i.species, mu_method, phi_method, n_seed = 100){
    
    # load TMB model
    version <- "binom"
    # compile(paste0(version,".cpp"))
    dyn.load(dynlib(version))
    
    # Read survey data
    load("../read_data/MARsurvey.rda")
    load("../read_data/MARgrid_1km.rda")
    
    # run repeated simulations
    model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
    
    # id of scenario
    save.i.scenario <- paste0("resampling/res-",i.species,"-",mu_method,"-",phi_method,"/")
    
    save.dir.output <- paste0("C:/Users/yinyi/Documents/workspace/simulation/",save.i.scenario,"/output/")
    save.dir.est_mu <- paste0("C:/Users/yinyi/Documents/workspace/simulation/",save.i.scenario,"/est_mu/")
    dir.create(save.dir.output, recursive = T)
    dir.create(save.dir.est_mu,recursive = T)
    dir.create(save.i.scenario, recursive = T)
    
    for(i.seed in 1:n_seed){
        set.seed(i.seed)
        simu <- simu_data(d = d, i.species, mu_method, phi_method, n_rep_stn = 2)
        save(simu, file = paste0(save.dir.output, "data-",i.seed,".rda"))
        for(i.model in model_vec){
            res <- try(fit_model(i.model, simu))
            if(!is.null(res)){if(!inherits(res, "try-error")){
                save(res, file = paste0(save.dir.output, "res-",i.seed,"-",i.model,".rda"))
            }}
        }
    }
    
    write_results(n_seed, save.i.scenario, save.dir.output, save.dir.est_mu)
    
}


library(foreach)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

comb <- expand.grid(1:3,1:2,c(14,23,201,204,10,11))
foreach(i = 1:nrow(comb), .packages = c("dplyr", "tidyr", "ggplot2", "zoo", "TMB")) %dopar% {
    run_model(i.species = comb[i,3],mu_method = comb[i,1],phi_method = comb[i,2])
}

