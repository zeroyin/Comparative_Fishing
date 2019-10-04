

rm(list = ls())
setwd("/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

library(TMB)
dyn.load("binom.dll")


# ---------------------------------
# mu and rho

jpeg(filename = "resampling/mu_rho_assumption.jpg",8,8,units = "in",res=200)
par(mfrow=c(3,2), mar=c(2,4,1,1))

len <- 3:100

# ratio of two logistic catchability functions
logistic <- function(len,L,k,x0) L/(1+exp(-k*(len-x0)))
q_A <- logistic(len, 1, 0.2, len[floor(length(len)*1/3)])
q_B <- logistic(len, 1, 0.1, len[floor(length(len)*1/2)])
mu <-  q_B/(q_A+q_B)
rho <- exp(boot::logit(mu))

plot(len, mu, type = "l", ylim = c(0,1),xlab = "Length")
lines(len, q_A, col="grey")
lines(len, q_B, col="grey",lty="dashed")
legend(80,0.3,legend = c("q_A","q_B"),col="grey",lty = c("solid","dashed"))

plot(len, rho, type = "l", xlab = "Length", ylim = c(0,max(rho)))
abline(h = 1, col = "grey")

# exponential or linear: a + b*len; a+b*exp(-c*len)
mu <- 1 - 1/(exp(0.7653+2.4198*exp(-0.2607*len))+1)
rho <- exp(boot::logit(mu))

plot(len, mu, type = "l", ylim = c(0,1),xlab = "Length")
plot(len, rho, type = "l", xlab = "Length", ylim = c(0,max(rho)))
abline(h = 1, col = "grey")

# station-level variation of shape
nstn <- 90
nlen <- length(len)
mu <- rho <- matrix(NA, nstn, nlen)
dimnames(mu)[[2]] <- dimnames(rho)[[2]] <- len
for(i.stn in 1:nstn){
    mu[i.stn,] <- 1 - 1/(exp(0.03*rlnorm(1,0-0.3^2/2,0.3)*len-2*rlnorm(1,0-0.3^2/2,0.3))+1)
    rho[i.stn,] <- exp(boot::logit(mu[i.stn,]))
}

matplot(len, t(mu), type = "l", col="black",lty="dashed", ylim = c(0,1), ylab = "mu", xlab = "Length")
matplot(len, t(rho), type = "l", col="black",lty="dashed", ylab = "rho", xlab = "Length")
abline(h = 1, col = "grey")

dev.off()


# ------------------------------
# phi

jpeg(filename = "resampling/phi_assumption.jpg",8,6,units = "in",res=200)

len <- 3:100
nlen <- length(len)
phi1 <- rep(0.4, nlen)  
phi_len <- seq(0:(max(len)+3))
phi2 <- dbeta((phi_len-min(phi_len))/(max(phi_len)-min(phi_len)),0.5,0.5)[len]

plot(len, phi1, type = "l", lty="dashed", ylim = c(0,3),xlab = "Length",ylab = "phi")
lines(len, phi2)
legend(80,3,legend = c("Scenario 1","Scenario 2"), lty = c("dashed","solid"))

dev.off()


# ------------------------------
# underlying density

mat2data <- function(x) 
    cbind(x, stn = 1:nrow(x)) %>%
    as.data.frame() %>%
    gather(len, n,-stn) %>%
    mutate(len = as.numeric(len))


species_list <- c(10,11,14,23,201,204)
mu_method = 2
phi_method = 1

for (i.species in species_list){
    i.scenario <- paste0(i.species,"-",mu_method,"-",phi_method)
    load(paste0("~/workspace/simulation/resampling/res-",i.scenario,"/output/data-1.rda"))
    
    bind_rows(
        mat2data(simu$N_A) %>% mutate(type="A"),
        mat2data(simu$N_B) %>% mutate(type="B")
    ) %>%
        filter(n > 0) %>%
        ggplot() +
        geom_tile(aes(as.factor(stn), len, fill = n)) +
        scale_fill_continuous(low = "white", high = "red", trans = "log10", na.value = "white", limits = c(1, NA)) +
        labs(fill = paste(i.scenario, ": Catch")) +
        theme(axis.text.x = element_blank(),
              axis.title = element_blank(),
              legend.position = "bottom",
              panel.border = element_rect(fill = NA),
              panel.background = element_blank(),
              panel.grid= element_blank()) +
        facet_wrap(~type, ncol = 2) 
    ggsave(filename = paste0("resampling/catch-",i.scenario,".jpg"),width = 12,height = 6)
    
}



# ---------------------------------
# gather plot residuals: targeted stations

model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
species_list <- c(10,11,14,23,201,204)
mu_method = 1
phi_method = 1

# residual plot
plotresid <- function(vec,mat1,mat2,lim,name){
    matplot(vec-0.1, t(mat1), pch=".", col="blue", ylim = c(-lim,lim), main = name, xlab = "Length", ylab = "Residual")
    matplot(vec+0.1, t(mat2), pch=".", col="red ", add = T)
    abline(h=0)
    points(vec-0.1, apply(mat1,2,function(x) median(x,na.rm = T)), col = "blue")
    points(vec+0.1, apply(mat2,2,function(x) median(x,na.rm = T)), col = "red")
    segments(vec-0.4,apply(mat1,2,function(x)quantile(x, 0.75, na.rm = T)),
             vec+0.2,apply(mat1,2,function(x)quantile(x, 0.75, na.rm = T)),
             col = "blue")
    segments(vec-0.4,apply(mat1,2,function(x)quantile(x, 0.25, na.rm = T)),
             vec+0.2,apply(mat1,2,function(x)quantile(x, 0.25, na.rm = T)),
             col = "blue")
    segments(vec-0.2,apply(mat2,2,function(x)quantile(x, 0.75, na.rm = T)),
             vec+0.4,apply(mat2,2,function(x)quantile(x, 0.75, na.rm = T)),
             col = "red ")
    segments(vec-0.2,apply(mat2,2,function(x)quantile(x, 0.25, na.rm = T)),
             vec+0.4,apply(mat2,2,function(x)quantile(x, 0.25, na.rm = T)),
             col = "red ")
}

for(i.species in species_list){
    
    scenario1 <- paste0(i.species,"-",mu_method,"-",phi_method)
    scenario2 <- paste0(scenario1, "-2+STN")
    
    load(paste0("~/workspace/simulation/resampling/res-",scenario1,"/output/data-1.rda"))
    nlen <- length(simu$len_list)
    len <- simu$len_list
    
    mu <- simu$m_mu
    rho <- exp(boot::logit(mu))
    
    aic_mat1 <- read.csv(paste0("resampling/res/",scenario1, "/aic_table.csv"))[,-1]
    aic_mat2 <- read.csv(paste0("resampling/res/",scenario2, "/aic_table.csv"))[,-1]
    resid2.rho <- resid1.rho <- resid2.mu <- resid1.mu <- matrix(NA,100,nlen) 
    for(i.seed in 1:100){
        aic.which1 <- model_vec[which.min(aic_mat1[i.seed,])]
        f1 <- paste0("~/workspace/simulation/resampling/res-",scenario1,"/output/res-",i.seed,"-",aic.which1,".rda")
        if(file.exists(f1)) {
            load(f1)
            resid1.mu[i.seed,] <- res$obj$report()$mean_mu-mu
            resid1.rho[i.seed,] <- res$obj$report()$mean_log_rho-log(rho)
        }
        rm(list = c("res","f1"))
        aic.which2 <- model_vec[which.min(aic_mat2[i.seed,])]
        f2 <- paste0("~/workspace/simulation/resampling/res-",scenario2,"/output/res-",i.seed,"-",aic.which2,".rda")
        if(file.exists(f2)){
            load(f2)
            resid2.mu[i.seed,] <-  res$obj$report()$mean_mu-mu
            resid2.rho[i.seed,] <- res$obj$report()$mean_log_rho-log(rho)
        }
        rm(list = c("res","f2"))
    }
    
    jpeg(paste0("resampling/resid-",scenario1,".jpg"), width = 10, height = 6, units = "in", res = 300)
    # par(mfrow = c(1,2),mar = c(2,2,2,1))
    # plotresid(simu$len_list, resid1.mu, resid2.mu, 0.25,name="mu")
    par(mar = c(2,4,2,1))
    plotresid(simu$len_list, resid1.rho, resid2.rho, log(2),name="log(rho)")
    abline(h=c(-log(1.5),log(1.5),-log(1.2),log(1.2)), col="gray")
    dev.off()
    
}

# ---------------------------------
# compare residuals: targeted stations

species_list <- c(10,11,14,23,201,204)
resid.mat <- matrix(NA,3*6,4)
for(i in 1:length(species_list)){
    i.species = species_list[i]
    for(j in 1:3){
        if(j==1){mu_method = 3; phi_method = 2}else if(j==2){mu_method = 2; phi_method = 1}else if(j==3){mu_method = 1; phi_method = 1}
        scenario1 <- paste0(i.species,"-",mu_method,"-",phi_method)
        scenario2 <- paste0(scenario1, "-2+STN")
        f1 <- read.csv(paste0("resampling/res/",scenario1,"/wmresid.csv"), row.names = 1)
        f2 <- read.csv(paste0("resampling/res/",scenario2,"/wmresid.csv"), row.names = 1)
        resid.mat[i+(j-1)*length(species_list),] <- as.numeric(cbind(f1["AIC",c(3,6)], f2["AIC",c(3,6)]))
        rm(list = c("f1","f2"))
    }
}

cbind(rbind(expand.grid(spec=species_list,mu=1,phi=1),
            expand.grid(spec=species_list,mu=2,phi=1),
            expand.grid(spec=species_list,mu=3,phi=2)),
      resid.mat[,c(1,3)], 
      log_rho=resid.mat[,3]-resid.mat[,1],
      resid.mat[,c(2,4)],
      sd_log_rho=resid.mat[,4]-resid.mat[,2]) %>%
  `colnames<-`(c("spec", "mu", "phi", 1,2,3,4, "log_rho", "sd_log_rho")) %>%
    mutate(spec=as.integer(spec),mu=as.integer(mu),phi=as.integer(phi)) %>%
    arrange(spec) %>%
    xtable::xtable(digits = 3) %>%
    print(include.rownames=FALSE)


# ---------------------------------
# 3 station scenarios
species_list <- c(10,201)
mu_list <- c(1,3)
phi_list <- c(1,2)
s_list <- expand.grid(spec = species_list, mu = mu_list, phi = phi_list)
resid.mat <- matrix(NA,nrow(s_list),4)
for(i in 1:nrow(s_list)){
        scenario1 <- paste0(s_list[i,1],"-",s_list[i,2],"-",s_list[i,3])
        scenario2 <- paste0(scenario1, "-3STN")
        f1 <- read.csv(paste0("resampling/res/",scenario1,"/wmresid.csv"),row.names = 1)
        f2 <- read.csv(paste0("resampling/res/",scenario2,"/wmresid.csv"),row.names = 1)
        resid.mat[i,] <- as.numeric(cbind(f1["AIC",c(3,6)], f2["AIC",c(3,6)]))
        rm(list = c("f1","f2"))
}


cbind(s_list,
      resid.mat[,c(1,3)], 
      log_rho=resid.mat[,3]-resid.mat[,1],
      resid.mat[,c(2,4)],
      p.20=resid.mat[,4]-resid.mat[,2]) %>%
  `colnames<-`(c("spec", "mu", "phi", 1,2,3,4, "log_rho", "p.20")) %>%
  mutate(spec=as.integer(spec),mu=as.integer(mu),phi=as.integer(phi)) %>%
  arrange(spec) %>%
  xtable::xtable(digits = 4) %>%
  print(include.rownames=FALSE)


# ---------------------------------
# tuncated size range scenarios

load("../read_data/MARsurvey.rda")
jpeg(filename = "truncate_size.jpg", 10,8,units = "in",res=300)
par(mfrow=c(2,1),mar=c(3,3,2,2))
plot(table(d[d$species==11,]$len),main = "11: nobs MAR survey")
abline(v=c(15,60),col="red")
plot(table(d[d$species==14,]$len),xlab="Length",main = "14: nobs MAR survey")
abline(v=c(10,50),col="red")
dev.off()


species_list <- c(11,14)
mu_list <- c(3)
phi_list <- c(1,2)
s_list <- expand.grid(spec = species_list, mu = mu_list, phi = phi_list)
resid.mat <- matrix(NA,nrow(s_list),4)
for(i in 1:nrow(s_list)){
  scenario1 <- paste0(s_list[i,1],"-",s_list[i,2],"-",s_list[i,3])
  scenario2 <- paste0(scenario1, "-TLen")
  f1 <- read.csv(paste0("resampling/res/",scenario1,"/wmresid.csv"),row.names = 1)
  f2 <- read.csv(paste0("resampling/res/",scenario2,"/wmresid.csv"),row.names = 1)
  resid.mat[i,] <- as.numeric(cbind(f1["AIC",c(3,6)], f2["AIC",c(3,6)]))
  rm(list = c("f1","f2"))
}

cbind(s_list,
      resid.mat[,c(1,3)], 
      log_rho=resid.mat[,3]-resid.mat[,1],
      resid.mat[,c(2,4)],
      p.20=resid.mat[,4]-resid.mat[,2]) %>%
  `colnames<-`(c("spec", "mu", "phi", 1,2,3,4, "log_rho", "p.20")) %>%
  mutate(spec=as.integer(spec),mu=as.integer(mu),phi=as.integer(phi)) %>%
  arrange(spec) %>%
  xtable::xtable(digits = 4) %>%
  print(include.rownames=FALSE)



# # ---------------------------------
# # Calculate additional statistics on residuals
# for(i.species in c(10,11,14,23,201,204)){
#     for(mu_method in 1:3){
#         for(phi_method in 1:2){
#             
#             save.i.scenario <- paste0(i.species,"-",mu_method,"-",phi_method)
#             save.dir.scenario <- paste0("resampling/res/",save.i.scenario,"/")
#             save.dir.output <- paste0("C:/Users/yinyi/Documents/workspace/simulation/resampling/res-",save.i.scenario,"/output/")
#             
#             aic_mat <- read.csv(paste0(save.dir.scenario, "aic_table.csv"))[,-1]
#             
#             load(paste0(save.dir.output, "data-1.rda"))
#             len <- simu$len_list
#             
#             model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4),"AIC")
#             wm_abs_mu <- rep(NA, length(model_vec)); names(wm_abs_mu) <- model_vec
#             p50_log_rho <- p20_log_rho <- wsd_log_rho <- wsd_mu <- wm_log_rho <- wm_mu <- wm_abs_mu
#             
#             load(paste0(save.dir.output, "enc.rda"))
#             
#             for(i.model in model_vec){
#                 
#                 load(paste0(save.dir.output,i.model,"-resid.rda"))
#                 
#                 if(sum(!is.na(rowMeans(residmu_mat)))>40){
#                     wm_abs_mu[i.model] <- mean(colMeans(residabsmu_mat*p_enc_mat/mean(p_enc_mat,na.rm = T), na.rm = T), na.rm = T)
#                     wm_mu[i.model] <- mean(abs(residmu_mat*p_enc_mat/mean(p_enc_mat,na.rm = T)), na.rm = T)
#                     wm_log_rho[i.model] <- mean(abs(residlogrho_mat*p_enc_mat/mean(p_enc_mat,na.rm = T)), na.rm = T)
#                     wsd_mu[i.model] <- mean(apply(residmu_mat,2,function(x)sd(x,na.rm = T))*colMeans(p_enc_mat)/mean(p_enc_mat,na.rm = T))
#                     wsd_log_rho[i.model] <- mean(apply(residlogrho_mat,2,function(x)sd(x,na.rm = T))*colMeans(p_enc_mat)/mean(p_enc_mat,na.rm = T))
#                     p20_log_rho[i.model] <- sum(abs(residlogrho_mat)<log(1.2))/sum(!is.na(residlogrho_mat))
#                     p50_log_rho[i.model] <- sum(abs(residlogrho_mat)<log(1.5))/sum(!is.na(residlogrho_mat))
#                     
#                 }
#                 
#                 rm("residmu_mat","residabsmu_mat","residlogrho_mat")
#                 
#             }
#             
#             
#             cbind(wm_abs_mu,wm_mu,wm_log_rho,wsd_mu,wsd_log_rho,p20_log_rho,p50_log_rho) %>%
#                 write.csv(file = paste0(save.dir.scenario, "wmresid.csv"))
#         }}}



# ---------------------------------
# summary table

s_list <- expand.grid(i.species=c(10,11,14,23,201,204), mu_method=1:3, phimethod=1:2)
res_resid <- res_summary <- list()
for(i in 1:nrow(s_list)){
    summary.path <- paste0("resampling/res/", paste(s_list[i,1],s_list[i,2],s_list[i,3],sep="-"),"/summary_stats.csv")
    resid.path <- paste0("resampling/res/", paste(s_list[i,1],s_list[i,2],s_list[i,3],sep="-"),"/wmresid.csv")
    res_summary[[i]] <- read.table(summary.path, row.names = 1, header = T, sep = ",")
    res_resid[[i]] <- read.table(resid.path, row.names = 1, header = T, sep = ",")["AIC",]
}

cbind(s_list,t(bind_cols(res_summary)),bind_rows(res_resid)) %>% 
    arrange(i.species) %>% 
    mutate(i.species=as.integer(i.species)) %>%
    xtable::xtable(digits = 3)


# ---------------------------------
# model comparison
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
s_list <- expand.grid(spec=c(10,11,14,23,201,204), mu=1:3, phi=1:2)
nbest_summary <- nconv_summary <- list()
for(i in 1:nrow(s_list)){
    summary.path <- paste0("resampling/res/", paste(s_list[i,1],s_list[i,2],s_list[i,3],sep="-"),"/aic_summary.csv")
    summary.d <- read.table(summary.path, row.names = 1, header = T, sep = ",")
    nconv_summary[[i]] <- summary.d[,1]
    nbest_summary[[i]] <- summary.d[,2]
}

model.summary <- cbind(s_list,t(rbind(bind_cols(nbest_summary))))
colnames(model.summary)[4:ncol(model.summary)] <- model_vec
model.summary %>% 
    arrange(spec) %>% 
    mutate(spec=as.integer(spec)) %>%
    xtable::xtable() %>%
    print(include.rownames=FALSE)






# # ---------------------------------
# # percentage of discrepency in the total number
# 
# model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))
# 
# i.species <- 10
# scenario <- "10-1-1"
# 
# load("../read_data/MARsurvey.rda")
# 
# load(paste0("~/workspace/simulation/resampling/res-",scenario,"/output/data-1.rda"))
# nlen <- length(simu$len_list)
# 
# # size-frenquency over martime during 2018
# n_vec <- d %>% filter(species == i.species) %>%
#     filter(year == 2018) %>%
#     group_by(len) %>%
#     summarise(n = sum(catch)/length(unique(d[d$year==2018,]$ID))) %>%
#     ungroup() %>%
#     complete(len=full_seq(simu$len_list,1), fill = list(n = 0)) %>%
#     mutate(n = na.fill(rollmean(n, 5, fill = NA), "extend"))
# 
# n_vec %>%
#     ggplot() +
#     geom_point(aes(len,n))
# 
# load("~/workspace/simulation/resampling/res-10-1-1/output/AIC-resid.rda")
# plot(exp(residlogrho_mat[1,]) * n_vec$n)




# ---------------------------------
# model comparison

res_spec <- expand.grid(i.species = c(10,11,14,23,201,204), mu_method = 3, phi_method = 2)
res_list <- list()
for(i in 1:nrow(res_spec)){
  resid.path <- paste0("resampling/res-",res_spec[i,1],"-",res_spec[i,2],"-",res_spec[i,3],"-2+STN/wmresid.csv")
  res_list[[i]] <- read.table(resid.path, row.names = 1, header = T, sep = ",")
}

res.summary <- cbind(res_spec,t(bind_cols(res_list)))
colnames(res.summary)[4:ncol(res.summary)] <- rownames(res_list[[1]])
res.summary[,c(1:3,ncol(res.summary)-3:0,4:(ncol(res.summary)-4))]

res.summary[,1:(ncol(res.summary)-4)] %>%
  as.data.frame() %>%
  mutate(scenario=paste0(i.species,"-",mu_method,"-",phi_method)) %>%
  select(-i.species,-mu_method,-phi_method) %>%
  gather(model, value,-scenario) %>%
  ggplot() +
  geom_tile(aes(model, scenario, fill= value)) +
  scale_fill_continuous(limits = c(0,NA)) +
  theme_bw()


res.summary[,1:(ncol(res.summary)-4)] %>%
  as.data.frame() %>% 
  gather(model, value,-i.species,-mu_method,-phi_method) %>%
  ggplot() +
  geom_point(aes(model, value, color = as.factor(i.species))) +
  geom_line(aes(model, value, color = as.factor(i.species), group = as.factor(i.species))) +
  facet_grid(mu_method~phi_method, scales = "free_y") +
  theme_bw()








# ---------------------------------
# NED2013


rm(list = ls())
setwd("/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/length_binom_models/")

dyn.load("../simulate_data/binom.dll")

species_list <- c(10,11,14,23,201,204)
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))

load("../read_data/data-NED2013.RData")


d.length %>% 
    filter(species %in% species_list) %>%
    group_by(species) %>%
    summarize(n=n_distinct(station))
        
read.csv("../length_binom_models/NED2013/aic_table.csv", header = F) %>%
    xtable::xtable() %>%
    print(include.rownames = F, include.colnames=F)


for(i in 1:length(species_list)){
    
    i.species <- species_list[i]
    b.len <- 1
    
    # organize data
    d <- d.length %>% 
        filter(species == i.species) %>%
        transmute(
            station = factor(station),
            vessel = factor(gear), 
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
    len <- unique(d$len)
    nlen = length(len)
    nstation = nlevels(d$station)
    data = list(
        N_A = d %>% filter(vessel == levels(d$vessel)[1]) %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
        N_B = d %>% filter(vessel == levels(d$vessel)[2]) %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix())
    
    aic_mat <- read.csv("NED2013/aic_table.csv")[1:13,]
    
    jpeg(paste0("NED2013/",i.species, "-est-mu.jpg"),res = 600, width = 10, height = 8, units = "in")
    matplot(len, t(data$N_A/(data$N_A + data$N_B)), type = "p", cex = 0.05, pch =20, col = "black", xlab = "len", ylab = NA)
    lines(len, colMeans(data$N_A/(data$N_A + data$N_B), na.rm = T), col = "orange")
    for(i.model in model_vec){
        res_file <- paste0("NED2013/",i.species,"-",i.model,".rda")
        if(file.exists(res_file)){
            load(res_file)
            est.mean_mu <- res$obj$report()$mean_mu
            col = ifelse(i.model == model_vec[which.min(aic_mat[,i+1])], "red", "blue")
            lines(len, est.mean_mu, type = "l",  col = col)
            text(x=len[1]-1, y=est.mean_mu[1], labels = i.model, col = col)
            rm("res", "res_file", "est.mean_mu", "col")
        }
    }
    legend("bottom", pch = c(1,NA,NA), lty = c(0,1,1),
           col = c("black","orange","black"),
           legend = c("obs mu/stn", "obs mu/mean", "true mu"))
    dev.off()
    
}



# ---------------------------------
# NED2005

rm(list = ls())
setwd("/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/length_binom_models/")

dyn.load("../simulate_data/binom.dll")

species_list <- c(10,11,14,23,201,204)
model_vec <- c(paste0("BB", 0:7),paste0("BI", 0:4))

load("../read_data/data-NED2005.RData")


d.length %>% 
  filter(species %in% species_list) %>%
  group_by(species) %>%
  summarize(n=n_distinct(station))


read.csv("../length_binom_models/NED2005/aic_table.csv", header = F) %>%
  xtable::xtable() %>%
  print(include.rownames = F, include.colnames=F)


for(i in 1:length(species_list)){
  
  i.species <- species_list[i]
  b.len <- 1
  
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
  len <- unique(d$len)
  nlen = length(len)
  nstation = nlevels(d$station)
  data = list(
    N_A = d %>% filter(vessel == levels(d$vessel)[1]) %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix(),
    N_B = d %>% filter(vessel == levels(d$vessel)[2]) %>% spread(len, catch) %>% select(-station, -vessel) %>% as.matrix())
  
  aic_mat <- read.csv("NED2005/aic_table.csv")[1:13,]
  
  jpeg(paste0("NED2005/",i.species, "-est-mu.jpg"),res = 600, width = 10, height = 8, units = "in")
  matplot(len, t(data$N_A/(data$N_A + data$N_B)), type = "p", cex = 0.05, pch =20, col = "black", xlab = "len", ylab = NA)
  lines(len, colMeans(data$N_A/(data$N_A + data$N_B), na.rm = T), col = "orange")
  for(i.model in model_vec){
    res_file <- paste0("NED2005/",i.species,"-",i.model,".rda")
    if(file.exists(res_file)){
      load(res_file)
      est.mean_mu <- res$obj$report()$mean_mu
      col = ifelse(i.model == model_vec[which.min(aic_mat[,i+1])], "red", "blue")
      lines(len, est.mean_mu, type = "l",  col = col)
      text(x=len[1]-1, y=est.mean_mu[1], labels = i.model, col = col)
      rm("res", "res_file", "est.mean_mu", "col")
    }
  }
  legend("bottom", pch = c(1,NA,NA), lty = c(0,1,1),
         col = c("black","orange","black"),
         legend = c("obs mu/stn", "obs mu/mean", "true mu"))
  dev.off()
  
}


