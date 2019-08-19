# ###############################################
# Simulate data for model testing
# Run a set of simulations, assess performance
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(ggplot2)


# ###############################################

# function to simulate data
simu_data <- function(area = 1000){
    
    # station list
    n_stn <- 20
    n_stn_strat <- rep(1, n_stn) # fixed number of stations per stratum
    stn_list <- seq(1, n_stn)
    
    # length list
    n_len <- 50
    len_list <- seq(11, 10 + n_len, length.out = n_len)
    
    # length composition per station
    len_comp <- matrix(NA, n_stn,n_len)
    for(i.stn in 1:n_stn){
        len_comp[i.stn, ] <- dnorm(seq(-3,3,length.out = n_len),(i.stn-n_stn/2)/4, 1+i.stn/10)
    }
    
    # population density per station
    dens <- rbeta(n_stn, 3,5)
    dens_mat <- dens * len_comp
    
    # Catchability
    s_A <- len_list/60 * 1.2
    q_A <- s_A/sum(s_A)
    
    s_B <-  len_list/60 * 1.4 + sqrt(len_list/60)
    q_B <- s_B/sum(s_B)
    
    #  true conversion: 
    mu <- q_A/(q_A+q_B)
    rho <- q_A/q_B
    
    # Catch at Length: area for scaling catch numbers
    area <- area
    N_A <- N_B <- matrix(NA, n_stn, n_len)
    
    # Poisson sampling variation
    for(i.stn in 1:n_stn){
        N_A[i.stn,] <- rpois(n = n_len, lambda = q_A*dens_mat[i.stn,]*area*rlnorm(1,-0.18,0.6))
        N_B[i.stn,] <- rpois(n = n_len, lambda = q_B*dens_mat[i.stn,]*area*rlnorm(1,-0.18,0.6))
    }
    
    simu <- list(N_A = N_A, N_B = N_B, len_list = len_list, mu = mu, rho = rho)
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

dir.create("res/output/",recursive = T)

# run repeated simulations
n_seed <- 100
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
for(i.seed in 1:n_seed){
    set.seed(i.seed)
    simu <- simu_data(area = 1000)
    save(simu, file = paste0("res/output/data-",i.seed,".rda"))
    for(i.model in model_vec){
        res <- try(fit_model(i.model, simu))
        if(!is.null(res)){if(!inherits(res, "try-error")){
            save(res, file = paste0("res/output/res-",i.seed,"-",i.model,".rda"))
        }}
    }
}


# organize results: AIC table

aic_mat <- matrix(NA, n_seed, length(model_vec))
colnames(aic_mat) <- model_vec

for(i.seed in 1:n_seed){
    for(i.model in model_vec){
        res_file <- paste0("res/output/res-",i.seed,"-",i.model,".rda")
        if(file.exists(res_file)){
            load(res_file)
                aic_mat[i.seed, i.model] <- res$aic
            rm("res", "res_file")
        }
    }
}

t(round(aic_mat, digits = 0))


round(aic_mat - apply(aic_mat,1,function(x){min(x,na.rm = T)}), digits = 0) %>%
    write.csv(file = "res/aic_table.csv")


aic.summary <- data.frame(
    nconv = as.integer(apply(aic_mat, 2, function(x){sum(!is.na(x))})),
    nbest = as.integer(colSums(aic_mat==apply(aic_mat,1,function(x){min(x,na.rm = T)}), na.rm = T)),
    row.names = model_vec)
    

xtable::xtable(t(aic.summary))
           


# organize results: est.mu

for(i.seed in 1:n_seed){
    
    load(paste0("res/output/data-",i.seed, ".rda"))
    
    jpeg(paste0("res/est_mu-",i.seed,".jpg"),res = 600, width = 10, height = 8, units = "in")
    matplot(simu$len_list, t(simu$N_A/(simu$N_A + simu$N_B)), type = "p", cex = 0.2, col = "black", xlab = "len", ylab = NA)
    lines(simu$len_list, colMeans(simu$N_A/(simu$N_A + simu$N_B), na.rm = T), col = "orange")
    lines(simu$len_list, simu$mu, col = "black")
    for(i.model in model_vec){
        res_file <- paste0("res/output/res-",i.seed,"-",i.model,".rda")
        if(file.exists(res_file)){
            load(res_file)
            est.mean_mu <- res$obj$report()$mean_mu
            lines(simu$len_list, est.mean_mu, type = "l",  col = "blue")
            text(x=simu$len_list[1]-1, y=est.mean_mu[1], labels = i.model,
                 col = ifelse(i.model == model_vec[which.min(aic_mat[i.seed,])], "red", "blue"))
            rm("res", "res_file", "est.mean_mu")
        }
    }
    legend("bottomright", pch = c(1,NA,NA), lty = c(0,1,1), 
           col = c("black","orange","black"),
           legend = c("obs mu/stn", "obs mu/mean", "true mu"))    
    dev.off()
    
    rm("simu")
}



load("res/output/data-1.rda")
len_list <- simu$len_list
rm("simu")
resid_mat <- matrix(NA, n_seed, length(len_list))
for(i.seed in 1:n_seed){
    aic_mat <- read.csv("res/aic_table.csv")[,-1]
    dat_file <- paste0("res/output/data-",i.seed, ".rda")
    res_file <- paste0("res/output/res-",i.seed,"-",model_vec[which.min(aic_mat[i.seed,])],".rda")
    load(dat_file)
    load(res_file)
    resid_mat[i.seed,] <- res$obj$report()$mean_mu - simu$mu
    rm("simu","res")
}



resid_mat %>%
    as.data.frame() %>%
    `colnames<-`(len_list) %>%
    gather() %>%
    `colnames<-`(c("len","resid")) %>%
    mutate(len = as.integer(len)) %>%
    ggplot(aes(x = len, y = resid)) +
    geom_point( color = "blue", size = 0.2, alpha = 0.2) +
    geom_hline(yintercept = 0, color = "black") +
    stat_summary(fun.y = median, geom="point", color = "blue") +
    theme_bw()
ggsave(filename = "res/resid_median-100simu.jpg", width = 10, height = 6)





