# ###############################################
# Simulate data for model testing
# test:
# interpolate density from survey
# Depreciated
# ###############################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/simulate_data//")

library(dplyr)
library(tidyr)
library(ggplot2)


# --------------------------------
# Simulate density


# (1) resampling within stratum: assume same density within stratum; 

load("../read_data/MARsurvey.rda")
d.cod <- d %>% 
    filter(species == 10) %>% 
    group_by(strat,len) %>%
    summarise(catch = sum(catch)) %>%
    ungroup %>%
    complete(strat, len=full_seq(len,1), fill = list(catch = 0))
    


# (2) size-spatial interpolation:



# (3) size-indepedent spatial interpolation:
library(sp)

load("../read_data/MARsurvey.rda")
d.catch <- d %>%
    filter(species == 11, len==35) %>% 
    group_by(ID) %>%
    mutate(catch = sum(catch)) %>% slice(1) %>%
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    PBSmapping::convUL() %>%
    `coordinates<-`(~ X + Y)

load("../read_data/MARgrid_10km.rda")
MARgrid <- MARgrid %>%
    filter(strat %in% c(440:498)) %>%#,paste0("5Z",1:9)
    mutate(X=lon, Y=lat) %>%
    `attr<-`("projection", "LL") %>%
    PBSmapping::convUL() %>%
    `coordinates<-`(~ X + Y)

library(gstat)

# use log(catch) instead of catch
# empirical distribution of the data is skewed 
hist(d.catch$catch)
# then the kriging estimators are sensitive to a few large data values
moments::skewness(log(d.catch$catch))
# not close to normal distribution
# transforamtion to stablize variance
shapiro.test(log(d.catch$catch))

# anisotropic
v.map <- variogram(log(catch) ~ 1, d.catch, map = TRUE,cutoff=600, width=600/18)
plot(v.map, col.regions = bpy.colors(64))

v.samp <- variogram(log(catch)~1, d.catch)
plot(v.samp)


krig <- automap::autoKrige(log(catch)~1, d.catch, MARgrid)


plot(krig)


ggplot() + 
    geom_point(data = as.data.frame(krig$krige_output),
               aes(x=X, y=Y, color=exp(var1.pred)))+
    scale_color_continuous(low = "white",high = "red", trans="log10", limits=c(1,NA)) + 
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())

ggplot() + 
    geom_point(data = as.data.frame(d.catch),
               aes(x=lon, y=lat, color=catch))+
    scale_color_continuous(low = "white",high = "red", trans="log10", limits=c(1,NA)) + 
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())













