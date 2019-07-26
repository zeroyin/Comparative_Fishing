# ###############################################
# Simulate data for model testing
# ###############################################


rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(ggplot2)


# ----- Presets -----#

# strata list
n_strat <- 20
strat_list <- seq(1, n_strat)

# station list
n_stn <- 10
# stn_strat <- 1 + rmultinom(1, n_station-n_strat, rep(1/n_strat,n_strat)) # random number of stations per stratum
n_stn_strat <- rep(1, n_stn) # fixed number of stations per stratum
stn_list <- seq(1, n_stn)

# length list
n_len <- 50
len_list <- seq(11, 10 + n_len, length.out = n_len)


# ----- Population Density -----#

# population density same within stratum
## use flat beta disribution
## use spatial distribution inferred from survey

dens <- rbeta(n_strat, 3,5)

# plot
hist(dens)

# length composition per stratum
## flat uniform distribution
## simulate a smooth curve
## use smoothed curve from survey

# len_comp <- rep(1,n_len)
# len_comp <- dbeta(len_list/70,1.8,2)

len_comp <- matrix(NA, n_stn,n_len)
for(i.stn in 1:n_stn){
    len_comp[i.stn, ] <- dnorm(seq(-3,3,length.out = n_len),(i.stn-n_stn/2)/4, 0.4 + i.stn/10)
}
    
# plot
matplot(len_list, t(len_comp), type = "l")


dens_mat <- dens * len_comp

matplot(len_list, t(dens_mat), type = "l")


# ----- Gear Catchability -----#

# cachability is more than gear selectivity 
## parametric selectivity models: logistic, lognormal etc
## designed selectivity 

s_A <-  dlnorm(len_list, meanlog = log(40), sdlog = 0.3)
# s_A <-  cos(len_list/10+20)+1
q_A <- s_A/sum(s_A)

s_B <-  dlnorm(len_list, meanlog = log(30), sdlog = 0.3)
q_B <- s_B/sum(s_B)

#  true conversion: 
rho <- q_A/q_B

# plot
matplot(len_list, cbind(q_A, q_B), type = "l", lty = c(1,2), col = "black")
legend("topleft", lty = c(1,2), legend = c("A", "B"))


plot(len_list, q_A/(q_A+q_B), type = "l")

plot(len_list, rho, type = "l")


# ----- Catch at Length -----#

# area for scaling catch numbers
area <- 1000
N_A <- N_B <- matrix(NA, n_stn, n_len)

# Poisson sampling variation
for(i.stn in 1:n_stn){
    N_A[i.stn,] <- rpois(n = n_len, lambda = q_A*dens_mat[i.stn,]*area*rlnorm(1,-0.18,0.6))
    N_B[i.stn,] <- rpois(n = n_len, lambda = q_B*dens_mat[i.stn,]*area*rlnorm(1,-0.18,0.6))
}

# combined catch 
N <- N_A + N_B


# plot

jpeg(paste(sep = "-","res/truth_catch.jpg"),
     res = 300, width = 8, height = 10, units = "in")
par(mfrow = c(2,1))
matplot(len_list, t(N_A), type = "l", lty = "dashed", col = "black")
lines(len_list, colMeans(N_A), col = "orange")
matplot(len_list, t(N_B), type = "l", col = "black")
lines(len_list, colMeans(N_B), col = "orange")
dev.off()

jpeg(paste(sep = "-","res/truth_mu.jpg"),
     res = 300, width = 8, height = 6, units = "in")
matplot(len_list, t(N_A/N), type = "p", cex = 0.2, col = "black")
lines(len_list, colMeans(N_A/N, na.rm = T), col = "orange")
lines(len_list, boot::inv.logit(log(rho)), col = "black")
legend("bottomright", pch = c(1,NA,NA), lty = c(0,1,1), col = c("black","orange","black"), c("mu by stn", "mean obs mu", "true mu"))    
dev.off()



# ########################################
# Fit models
# ########################################



# load TMB model
library(TMB)
version <- "binom_1"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))



# function to fit tmb model, given species and model name
fit_model <- function(i.model, N_A, N_B, lenseq){
    
    # model related: nll index and df index
    ind_nll<- 1 + switch(
        i.model,
        "ZB2"=c(5,12), "ZB3"=c(0,5,12), 
        "BI0"=10,"BI1"=c(4,10),"BI2"=c(0,10),"BI3"=c(0,4,10),"BI4"=c(0,2,3,10),
        "BB0"=11,"BB1"=c(4,11),"BB2"=c(0,11),"BB3"=c(0,1,11),"BB4"=c(0,4,11),"BB5"=c(0,1,4,11),"BB6"=c(0,2,3,11),"BB7"=c(0,1,2,3,11)
    )
    ind_df <- switch(
        i.model,
        "ZB2"=4,"ZB3"=5,
        "BI0"=1,"BI1"=2,"BI2"=3,"BI3"=4,"BI4"=7,
        "BB0"=2,"BB1"=3,"BB2"=4,"BB3"=6,"BB4"=5,"BB5"=7,"BB6"=8,"BB7"=10
    )
    
    N <- N_A + N_B
    
    # basis and penalty matrices for cubic spline: default to 10 knots
    library(mgcv)
    cs <- smooth.construct(
        object = s(len, bs = "cr"),
        data = cbind(len = seq(1, 70, 1), catch = 2) %>% as.data.frame(),
        # data = cbind(len = len_list, catch = colSums(N)) %>% as.data.frame(),
        knots = data.frame(knots = len_list)
    )
    
    n_f <- 2
    n_r <- cs$df - n_f
    eigende <- eigen(cs$S[[1]])
    
    # input for TMB
    nlen = ncol(N)
    nstation = nrow(N)
    data = list(
        A = N_A,
        B = N_B,
        offset = matrix(0, nrow = nrow(N), ncol = ncol(N)),
        Xf = cs$X %*% eigende$vectors[,1:n_f+n_r],
        Xr = cs$X %*% eigende$vectors[,1:n_r],
        d = eigende$value[1:n_r],
        idist = switch(substr(i.model, 1, 2), "BI"=0, "BB"=1, "ZB"=2)
    )

    parameters = list(
        beta = rep(0, n_f),
        b = rep(0, n_r),
        log_s_b = log(0.1),
        gamma = rep(0, n_f),
        g = rep(0, n_r),
        log_s_g = log(10),
        delta = matrix(0, nstation, n_f),
        chol_delta = c(1,1,1), # use chol decomp in vector form
        epsilon = matrix(0, nstation, n_r),
        log_s_epsilon = log(10),
        beta_0 = 0,
        gamma_0 = 0,
        delta_0 = rep(0, nstation),
        log_sigma_delta_0 = 0,
        p = matrix(0.5, nstation, nlen),
        log_p_s1 = log(2),
        log_p_s2 = log(2)
    )
    
    # model specifications
    if(i.model == "BB7"){
        map = list(
            beta_0 = factor(NA),
            gamma_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB6"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            beta_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB5"){
        map = list(
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            beta_0 = factor(NA),
            gamma_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB4"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            beta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB3"){
        map = list(
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            beta_0 = factor(NA),
            gamma_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB2"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            beta_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB1"){
        map = list(
            beta = factor(rep(NA, n_f)),
            b = factor(rep(NA, n_r)),
            log_s_b = factor(NA),
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BB0"){
        map = list(
            beta = factor(rep(NA, n_f)),
            b = factor(rep(NA, n_r)),
            log_s_b = factor(NA),
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BI4"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            beta_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            gamma_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BI3"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            beta_0 = factor(NA),
            gamma_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BI2"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            beta_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            gamma_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BI1"){
        map = list(
            beta = factor(rep(NA, n_f)),
            b = factor(rep(NA, n_r)),
            log_s_b = factor(NA),
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            gamma_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "BI0"){
        map = list(
            beta = factor(rep(NA, n_f)),
            b = factor(rep(NA, n_r)),
            log_s_b = factor(NA),
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            gamma_0 = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            log_sigma_delta_0 = factor(NA),
            p = factor(matrix(NA, nstation, nlen)),
            log_p_s1 = factor(NA),
            log_p_s2 = factor(NA)
        )
    }else if(i.model == "ZB3"){
        map = list(
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            gamma_0 = factor(NA),
            log_sigma_delta_0 = factor(NA)
        )
    }else if(i.model == "ZB2"){
        map = list(
            beta = factor(rep(NA, n_f)),
            b = factor(rep(NA, n_r)),
            log_s_b = factor(NA),
            gamma = factor(rep(NA, n_f)),
            g = factor(rep(NA, n_r)),
            log_s_g = factor(NA),
            delta = factor(matrix(NA, nstation, n_f)),
            chol_delta = factor(c(NA,NA,NA)),
            epsilon = factor(matrix(NA, nstation, n_r)),
            log_s_epsilon = factor(NA),
            delta_0 = factor(rep(NA, nstation)),
            gamma_0 = factor(NA),
            log_sigma_delta_0 = factor(NA)
        )
    }
    
    
    obj = MakeADFun(data=data,
                    parameters=parameters,
                    map = map,
                    DLL=version,
                    random = c("b", "g", "delta", "epsilon", "delta_0", "p"),
                    silent = T)
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    
    # check convergence, maximum gradient and positive definite
    if(exists("opt")&!opt$convergence){
        gra <- obj$gr()
        hes <- eigen(optimHess(par=opt$par, fn=obj$fn, gr=obj$gr))$values
        if(max(abs(gra)) < 0.1 & min(hes) > -0.1){
            aic <- 2*sum(obj$report()$nll[ind_nll]) + 2*ind_df
            rep <- try(sdreport(obj))
            res <- list(obj = obj, opt = opt, rep = rep, aic = aic, gra = gra, hes = hes)
            save(res, file = paste0("res/", i.model,".rda"))
            
            # estimate and std of mu and phi and rho
            est <- summary(rep, "report")[,"Estimate"]
            std <-  summary(rep, "report")[,"Std. Error"]
            est.mean_mu <- est[names(est) == "mean_mu"]
            std.mean_mu <- std[names(std) == "mean_mu"]
            
            jpeg(paste(sep = "-","res/estimates",i.model,"CI95_zscore.jpg"),
                 res = 300, width = 8, height = 6, units = "in")
            plot(lenseq, est.mean_mu, ylim = c(0,1), type = "l", col = "orange")
            lines(lenseq, boot::inv.logit(log(rho)), col = "black")
            for(i in 1:nstation){lines(lenseq, obj$report()$mu[i,], col = "gray")}
            lines(lenseq, est.mean_mu + 1.96*std.mean_mu, col = "blue", lty = "dashed")
            lines(lenseq, est.mean_mu - 1.96*std.mean_mu, col = "blue", lty = "dashed")
            dev.off()
            
            return(res)
        }
    }
}



# ---------------------------------------------------------
# run a suite of models given species and model name
# ---------------------------------------------------------

# Note : false convergence filtering is removed

model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4),paste0("ZB", 2:3), "GB")
for(i.model in model_vec){
        res <- fit_model(i.model, N_A = N_A, N_B = N_B, lenseq = len_list)
}


# AIC table
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4),paste0("ZB", 2:3))
aic_mat <- rep(NA, length(model_vec))
names(aic_mat) <- model_vec
for(i_model in 1:length(model_vec)){
    res_file <- paste0("res/",model_vec[i_model],".rda")
    if(file.exists(res_file)){
        load(res_file)
        aic_mat[i_model] <- res$aic
        rm("res", "res_file")
    }
}

t(round(aic_mat, digits = 0)) %>%
    write.csv(file = "res/aic_table.csv")



