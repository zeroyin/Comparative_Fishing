# #######################################
# Read data: NED2013
# Note: catch already scaled by tow dist
# #######################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)

# read tables
library("xlsx")
data.path <- "../../Data/20190510-Don-Survey/compare sumr2013 december 2013dc.xlsx"

options(stringsAsFactors = F)
data.catch <- read.xlsx(data.path, sheetIndex = 4, header = T, endRow = 15287) 
data.set <- read.xlsx(data.path, sheetIndex = 1, header = T, startRow = 8) 


# paired by station

d <- data.set %>% 
    filter(TYPE == 5) %>%
    transmute(mission = MISSION, 
              setno = as.integer(SETNO),
              station = as.integer(STAtion),
              gear = as.integer(GEAR)) %>%
    inner_join(
        data.catch %>% transmute(
            mission, 
            setno = as.integer(setno),
            station = as.integer(station),
            species = as.integer(spec),
            name = comm,
            len = flen,
            c = clen),
        by = c("mission", "setno", "station")
    )

save(d, file = "data-NED2013.RData")


# Check data
library(ggplot2)

# Table: number of obs by species
d %>% 
    count(species, name) %>%
    print(n=Inf)


# figure: catch by length by species
d %>%
    ggplot() +
    geom_line(aes(len, c,group = as.factor(station)), color = "gray") +
    stat_summary(aes(len, c), fun.y=mean, geom="line", colour="orange", size = 1) +
    facet_grid(species~gear, scales = "free")
ggsave("data_description/NED2013/catch_compare_by_len_species-all_species.pdf", width = 6, height = 100, limitsize = FALSE)


# sample rho
d %>% 
    select(species, station, gear, len, c) %>%
    spread(gear, c, fill = NA) %>%
    mutate(ratio = `15`/`9`) %>%
    ggplot(aes(len, log10(ratio), group = station)) +
    geom_line(color = "gray") +
    geom_point(size = 0.1) +
    geom_abline(intercept = 0, slope = 0) +
    theme_bw() +
    facet_wrap(~species, scales = "free", nrow = 10)
ggsave("data_description/NED2013/sample_rho-all_species.pdf", width = 40, height = 30)



# species 11
d %>% 
    filter(species == 11) %>%
    ggplot()+
    geom_line(aes(len, c, group = station), color = "gray") +
    stat_summary(aes(len, c), fun.y=median, geom="line", colour="orange") +
    stat_summary(aes(len, c), fun.y=mean, geom="line", colour="blue") +
    facet_wrap(~gear, ncol = 1) +
    theme_bw() 




