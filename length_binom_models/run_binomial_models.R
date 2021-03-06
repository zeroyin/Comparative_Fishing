# ###############################################
# run a suite of binomial and beta binomial models:
# summarize and compare results
# note: NED2005: same gear (9) different vessel
# ###############################################


rm(list = ls())
setwd("C://Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/length_binom_models/")

library(dplyr)
library(tidyr)
library(ggplot2)

# load data
load("../read_data/data-NED2005.RData")

# load TMB model
library(TMB)
version <- "binom"
compile(paste0(version,".cpp"))
dyn.load(dynlib(version))


# function to fit tmb model, given species and model name
fit_model <- function(i.model, i.species, b.len = 1, d.length){
    
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
    
    # organize data
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
    lenseq <- unique(d$len)
    
    # data for offset: same within a tow, use zeros for now
    d.offset <- d %>% 
        distinct(station) %>%
        mutate(offset = 0)
    
    # basis and penalty matrices for cubic spline: default to 10 knots
    library(mgcv)
    cs <- smooth.construct(
        object = s(len, bs = "cr"),
        data = data.frame(len = seq(min(lenseq)-10,max(lenseq)+10,1)),
        knots = data.frame(knots = lenseq)
    )
    
    n_f <- 2
    n_r <- cs$df - n_f
    eigende <- eigen(cs$S[[1]])
    
    # input for TMB
    nlen = length(lenseq)
    nstation = nlevels(d$station)
    data = list(
        A = d %>% filter(vessel == levels(d$vessel)[1]) %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
        B = d %>% filter(vessel == levels(d$vessel)[2]) %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
        offset = outer(d.offset$offset,rep(1,length(lenseq))),
        Xf = cs$X %*% eigende$vectors[,1:n_f+n_r],
        Xr = cs$X %*% eigende$vectors[,1:n_r],
        d = eigende$value[1:n_r],
        idist = switch(substr(i.model, 1, 2), "BI"=0, "BB"=1, "ZB"=2)
    )
    parameters = list(
        beta = rep(0, n_f),
        b = rep(0, n_r),
        log_s_b = log(10),
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
    if(exists("opt")){
        if(!opt$convergence){
            gra <- obj$gr()
            hes <- eigen(optimHess(par=opt$par, fn=obj$fn, gr=obj$gr))$values
            if(max(abs(gra)) < 0.1 & min(hes) > 0){
                aic <- 2*sum(obj$report()$nll[ind_nll]) + 2*ind_df
                rep <- try(sdreport(obj))
                res <- list(obj = obj, opt = opt, rep = rep, aic = aic, gra = gra, hes = hes)
                save(res, file = paste0("res/", i.species, "-",i.model,".rda"))
                
                # estimate and std of mu and phi and rho
                est <- summary(rep, "report")[,"Estimate"]
                std <-  summary(rep, "report")[,"Std. Error"]
                est.mean_mu <- est[names(est) == "mean_mu"]
                std.mean_mu <- std[names(std) == "mean_mu"]
                est.mean_phi <- est[names(est) == "mean_phi"]
                std.mean_phi <- std[names(std) == "mean_phi"]
                est.mean_log_rho <- est[names(est) == "mean_log_rho"]
                std.mean_log_rho <- std[names(std) == "mean_log_rho"]
                
                jpeg(paste(sep = "-","res/estimates",i.model,"species",i.species,"CI95_zscore.jpg"),
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
                
                return(res)
            }
        }
    }
    
}



# ---------------------------------------------------------
# run a suite of models given species and model name
# ---------------------------------------------------------

species_vec <- c(10, 11, 23, 14, 201, 204)
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4), paste0("ZB", 2:3))
for(i.model in model_vec){
    for(i.species in species_vec){
        res <- fit_model(i.model, i.species, d.length = d.length)
    }
}


# AIC table
species_vec <- c(10, 11, 23, 14, 201, 204)
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))#, paste0("ZB", 2:3), "GB")
aic_mat <- matrix(NA, length(species_vec), length(model_vec), dimnames = list(species_vec, model_vec))
for(i.species in 1:length(species_vec)){
    for(i.model in 1:length(model_vec)){
        res_file <- paste0("res/", species_vec[i.species], "-",model_vec[i.model],".rda")
        if(file.exists(res_file)){
            load(res_file)
            aic_mat[i.species, i.model] <- res$aic
            rm("res", "res_file")
        }
    }
}

t(round(aic_mat - apply(aic_mat, MARGIN = 1, FUN = function(x) min(x, na.rm = T)), digits = 0)) %>%
    write.csv(file = "res/aic_table.csv")


t(aic_mat)

