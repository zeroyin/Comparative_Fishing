# #########################################
# test beta binomial model
# #########################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)
library(ggplot2)

load("data-NED2016.RData")

# -----------------------------------------
# data for model 

d.length %>% count(species, name) %>% print(n=Inf)

i.species <- 11

d <- d.length %>% 
    filter(species == i.species) %>%
    transmute(station, gear, len, c) %>%
    spread(gear, c, fill = 0) %>%
    mutate(len = (floor(len/5))) %>% # length grouping
    group_by(len) %>% 
    summarise(A = (sum(`9`)), B = (sum(`15`)))%>%
    mutate(N = A+B) %>%
    complete(len = full_seq(len, 1), fill = list(A=0,B=0,N=0)) %>%
    mutate(A = A/min(N[N!=0])*2,
           B = B/min(N[N!=0])*2,
           N = N/min(N[N!=0])*2) %>% # for gamlss computation, blow up min number to 2
    print(n=Inf)


# -----------------------------------------
# fit model
# overdispersion (of a beta distribution) only exists when trials > 1
# higher degree of freedom may not be satified given "bad" data

library(gamlss)
rm("fit")
fit <- gamlss(
    formula = cbind(A, B) ~ cs(len),
    sigma.formula = ~ 1,
    family = BB(mu.link = "logit", sigma.link = "log"),
    data = d)

summary(fit)


# -----------------------------------------
# results plot

jpeg(filename = paste0("beta-binom/estimates_species-", i.species, ".jpg"), 
     units = "in", res = 200, width = 4, height = 6)
fittedPlot(fit, x = d$len)
dev.off()


# Proportion of A
jpeg(filename = paste0("beta-binom/mu_CI95(fast_simu)_species-", i.species, ".jpg"),
     units = "in", res = 200, width = 6, height = 5)
plot(d$len, d$A/d$N, ylim = c(0, 1), ylab = "Prop. of Gear 9")
lines(d$len, fit$mu.fv)
# computation of CI is based on simulation of estimated BB dist
# median more robust than mean
lines(d$len, qBB(0.5, mu = fit$mu.fv, sigma = fit$sigma.fv, bd = 1000)/1000, col = "blue")
lines(d$len, qBB(0.05, mu = fit$mu.fv, sigma = fit$sigma.fv, bd = 1000)/1000, col = "blue", lty = "dashed")
lines(d$len, qBB(0.95, mu = fit$mu.fv, sigma = fit$sigma.fv, bd = 1000)/1000, col = "blue", lty = "dashed")
dev.off()

# conversion factor: rho
jpeg(filename = paste0("beta-binom/rho_CI95(fast_simu)_species-", i.species, ".jpg"),
     units = "in", res = 200, width = 6, height = 5)
plot(d$len, d$A/d$B, ylab = "rho: Gear 9 / Gear 15")
lines(d$len, fit$mu.fv/(1-fit$mu.fv))
abline(a = 1, b = 0)
# CI approximated by lognormal given variance of logit(P)
lines(d$len, fit$mu.fv/(1-fit$mu.fv)-1.96*sqrt(fit$mu.var) * (fit$mu.fv/(1-fit$mu.fv)), col = "blue", lty = "dashed")
lines(d$len, fit$mu.fv/(1-fit$mu.fv)+1.96*sqrt(fit$mu.var) * (fit$mu.fv/(1-fit$mu.fv)), col = "blue", lty = "dashed")
dev.off()

# -----------------------------------------
# spatial exploratory
# map setting
base_map <- ggplot() +
    borders("world", colour="gray50", fill = "gray90", xlim = c(-67,-65), ylim = c(43,45)) +
    coord_map(xlim = c(-63,-57), ylim = c(43, 46)) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          legend.position = "bottom")

base_map +
    geom_point(data = d.length %>%
                   filter(species == i.species) %>% 
                   mutate(len = as.factor(floor(len/10))),
               aes(lon, lat, color = gear, size = c), shape = 21) +
    facet_wrap(~len)

