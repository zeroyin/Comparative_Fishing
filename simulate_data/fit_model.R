# function to fit model, given model name


fit_model <- function(i.model, simu){
    
    # attach data to environment
    N_A <- simu$N_A
    N_B <- simu$N_B
    len_list <- simu$len_list
    
    # # remove trivial stations and trivial length
    i <- rowSums(N_A+N_B)>0
    # j <- colSums(N_A+N_B)>0
    # 
    N_A <- N_A[i,]
    N_B <- N_B[i,]
    # len_list <- len_list[j]
    
    
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
    
    # total catch
    N <- N_A + N_B
    
    # basis and penalty matrices for cubic spline: default to 10 knots
    library(mgcv)
    cs <- smooth.construct(
        object = s(len, bs = "cr"),
        data = cbind(len = seq(min(len_list)-5, max(len_list)+5, 1), catch = 2) %>% as.data.frame(),
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
        chol_delta = c(1,1,1), 
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
    
    # model specifications: mapping
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
    
    # run model
    library(TMB)
    obj = MakeADFun(data=data,
                    parameters=parameters,
                    map = map,
                    DLL="binom",
                    random = c("b", "g", "delta", "epsilon", "delta_0", "p"),
                    silent = T)
    opt <- try(nlminb(obj$par,obj$fn,obj$gr))
    
    # check convergence, maximum gradient and positive definite
    if(exists("opt")){
        if(!inherits(opt, "try-error")){
            # if(!opt$convergence){
            gra <- obj$gr()
            if(max(abs(gra)) < 0.1){
                hes <- eigen(optimHess(par=opt$par, fn=obj$fn, gr=obj$gr))$values
                if(min(hes) > 0){
                    aic <- 2*sum(obj$report()$nll[ind_nll]) + 2*ind_df
                    if(aic < 10^10){
                        # rep <- try(sdreport(obj))
                        res <- list(obj = obj, opt = opt, aic = aic)
                        return(res)
                    }}}}}
}

