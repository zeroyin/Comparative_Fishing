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
options(stringsAsFactors = F)
data.set <- read.xlsx(
    file = "../../Data/20190510-Don-Survey/compare sumr2013 december 2013dc.xlsx",
    sheetIndex = 1, header = T, startRow = 8) 

# catch details table updated 20190613
data.catch <- read.xlsx(
    file = "../../Data/20190510-Don-Survey/catch_details_NED2013.xlsx",
    sheetIndex = 1, header = F) 


# paired by station

d <- data.set %>% 
    filter(TYPE == 5) %>%
    transmute(mission = MISSION, 
              setno = as.integer(SETNO),
              station = as.integer(STAtion),
              gear = as.integer(GEAR)) %>%
    inner_join(
        data.catch %>% 
            transmute(
                mission = X1, 
                setno = as.integer(X2),
                species = as.integer(X3),
                name = X4,
                size = as.factor(X5),
                len = as.integer(X6),
                c = X7) %>%
            group_by(mission, setno, species, name, len) %>%
            summarise(catch = sum(c)) %>%
            ungroup(),
        by = c("mission", "setno")
    )

save(d, file = "data-NED2013.RData")


# Check data
library(ggplot2)

# Table: number of obs by species
d %>% 
    count(species, name) %>%
    print(n=Inf)


# figure: catch by length by species
p <- d %>% 
    mutate(species = paste(species,name,sep = "-")) %>%
    ggplot() +
    geom_line(aes(len, catch, group = as.factor(station)), color = "gray") +
    stat_summary(aes(len, catch), fun.y=mean, geom="line", colour="orange", size = 1) +
    facet_grid(species~gear, scales = "free_y") +
    theme_bw()
ggsave(filename = "data_description/NED2013/catch_compare_by_len_species-all_species.pdf", 
       plot = p, width = 10, height = 150, limitsize = FALSE)


# sample rho
p <- d %>% 
    mutate(species = paste(species,name,sep = "-")) %>%
    select(species, station, gear, len, catch) %>%
    spread(gear, catch, fill = NA) %>%
    mutate(ratio = `15`/`9`) %>%
    ggplot(aes(len, log10(ratio), group = station)) +
    geom_line(color = "gray") +
    geom_point(size = 0.1) +
    geom_abline(intercept = 0, slope = 0) +
    facet_wrap(~species, scales = "free", ncol = 1) +
    theme_bw()
ggsave(filename = "data_description/NED2013/sample_rho-all_species.pdf",
       plot = p, width = 5, height = 150, limitsize = FALSE)



# species 11
i.species <- 11

p <- d %>% 
    filter(species == i.species) %>%
    ggplot()+
    geom_line(aes(len, catch, group = station), color = "gray") +
    stat_summary(aes(len, catch), fun.y=median, geom="line", colour="orange") +
    stat_summary(aes(len, catch), fun.y=mean, geom="line", colour="blue") +
    facet_wrap(~gear, ncol = 2) +
    theme_bw() 
ggsave(filename = paste0("data_description/NED2013/paired_catch-species_",i.species,".pdf"),
       plot = p, width = 8, height = 6)




