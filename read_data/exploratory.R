# ####################################################
# Species similarity
# ####################################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/read_data/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)

load("data-NED2005.RData")
str(d.length)

d <- d.length %>% 
    filter(species %in% c(10,11,14,23,201,204)) %>%
    transmute(
        station = station,
        vessel = factor(vessel), 
        ilen = (floor(len/5)*5), 
        species = factor(species),
        catch = catch
    ) %>%
    mutate(ilen = ifelse(ilen>60,60,ilen)) %>%
    group_by(species, station, vessel, ilen) %>%
    summarise(catch = sum(catch)) %>%
    ungroup() %>%
    complete(species, ilen, station, vessel, fill = list(catch = 0))

ggplot(d %>% spread(vessel, catch)) +
    geom_point(aes(NED,TEL,color = ilen, shape = species)) +
    scale_x_log10() + scale_y_log10() +
    theme_bw()


# x <- d %>% 
#     spread(vessel, catch) %>% 
#     nest(-c(species,ilen)) %>% 
#     mutate(fit = map(data, ~ lm(NED ~ TEL, data = .)),
#            intercept = map_dbl(fit, function(x) x$coefficients[[1]]),
#            slope = map_dbl(fit, function(x) x$coefficients[[2]]),
#            results = map(fit, augment)) %>% 
#     unnest(results)

d %>% 
    spread(vessel, catch) %>%
    ggplot(aes(x = TEL, y = NED)) +
    geom_point(aes(color = station)) +
    geom_smooth(method = "lm", formula = y ~ x + 0) +
    geom_abline(slope = 1, color = "gray") +
    facet_grid(species ~ ilen) +
    scale_x_log10(limits = c(1,1000)) + scale_y_log10(limits = c(1,1000)) +
    theme_bw()+
    theme(panel.grid = element_blank())
