# #######################################
# Read data: NED2016
# #######################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)

# read tables
library("xlsx")
data.path <- "../../Data/20190510-Don-Survey/NED2016116.xlsx"
data.set <- read.xlsx(data.path, sheetIndex = 1, header = T) 
data.catch <- read.xlsx(data.path, sheetIndex = 2, header = T) 
data.length <- read.xlsx(data.path, sheetIndex = 3, header = T, endRow = 3163)

# map settings
library(maps)
library(ggplot2)
base_map <- ggplot() +
    borders("world", colour="gray50", fill = "gray90", xlim = c(-67,-65), ylim = c(43,45)) +
    coord_fixed(ratio = 1.2, xlim = c(-64,-56), ylim = c(43, 46)) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        legend.position = "bottom"
    )

# #########################################
# set info
# #########################################

# check data

ggplot(data.set) +
    geom_histogram(aes(x = DIST)) +
    facet_wrap(~as.factor(GEAR), ncol = 1)

ggplot() +
    geom_point(data = data.set, aes(as.factor(STATION),DUR, color = as.factor(GEAR))) +
    geom_text(data = data.set %>% filter(TYPE!=5), aes(as.factor(STATION),DUR+1, label = TYPE, color = as.factor(GEAR))) 

# organize

d.set <- data.set %>%
    filter(
        TYPE == 5,# successful compartative tows
        !is.na(DIST)# valid distance data
    )%>% 
    transmute(
        mission = MISSION,
        gear = as.factor(GEAR),
        stratum = as.integer(STR), 
        station = as.integer(STATION),
        setno = as.integer(SETNO),
        dmin = DMIN,
        duration = DUR,
        distance = DIST, 
        lon = - (floor(SLONG/100) + SLONG%%100/60), # convert lat-long
        lat = floor(SLAT/100) + SLAT%%100/60
    ) %>%
    filter(
        station!=55 # remove one tow: short distance/duration
    )


# summarize

d.set %>%
    group_by(gear) %>%
    summarise(
        m_dist = mean(distance), 
        n_station = n()
    )
    

# visualize

base_map +
    geom_point(data = d.set, aes(lon, lat, shape = gear, color = gear)) +
    geom_text(data = d.set %>% group_by(station) %>% slice(1), aes(lon, lat+0.1, label = station), size = 2)
ggsave("tow_locations_map.pdf", width = 8, height = 5)


# #########################################
# catch
# #########################################


# check data

data.catch %>% 
    count(spec, comm) %>%
    ggplot() +
    geom_histogram(aes(n))


# organize

d.catch <- data.catch %>%
    transmute(
        mission,
        gear = as.factor(gear),
        stratum = as.integer(strat), 
        station = as.integer(station),
        species = as.factor(spec),
        name = comm,
        n = totno,
        w = totwgt
    )


# visualize: catch by tow, by species, by gear

d.catch %>%
    ggplot() +
    geom_point(aes(as.factor(station), n, color = gear, shape = gear)) +
    facet_wrap(~paste(species, name), scales = "free", ncol = 4) +
    theme_bw() + theme(legend.position = "top")
ggsave(filename = "paired_catch_n_per_species.pdf", dpi = 300, width = 10, height = 40)



# #########################################
# Join tables with length frequency
# #########################################


d.length <- data.length %>% 
    transmute(
        mission,
        station = as.integer(station),
        setno = as.integer(setno),
        species = as.integer(spec),
        name = as.factor(comm),
        flen, clen
        ) %>%
    inner_join(d.set, by = c("mission", "station", "setno")) %>%
    mutate(
        len = flen,
        catch = clen/distance # adjust catch for tow distance
    )


# visualize: length frenquency by tow, by species, by gear

ggplot(d.length) +
    geom_line(aes(len, c, group = station), color = "gray") +
    stat_summary(aes(len, c), fun.y=mean, geom="line", colour="red", alpha = 0.3, size = 1.2) +
    facet_grid(paste(species, substr(name,1,12))~gear, scales = "free_y", switch = "y") +
    theme_bw()
ggsave(filename = "catch_length_per_species.pdf", dpi = 300, width = 8, height = 80, limitsize = F)



# #########################################
# Save data frame
# #########################################

save(file = "data-NED2016.RData", d.length)


