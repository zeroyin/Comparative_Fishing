# #############################################
# Read data: NED2005
# 20190628: obsolete as Don indicates that
# one of the four paired trips is unsuccessful
# and updates me with station information
# #############################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/test_size_gam/")

library(dplyr)
library(tidyr)
library(purrr)

# read tables
library(readxl)
data.path <- "../../Data/20190604-Jamie-NED2005/"

# -----------------------------------------------
# Mission table: set information, and for pairing tows

path.set <- paste0(data.path,"2005_2006 Comparative Missions (All Sets).xlsx")
sheet.set <- excel_sheets(path.set)[1:8]
data.set <- map_df(sheet.set, ~ read_excel(
    path.set, sheet = .x,
    col_names = T,
    range = cell_cols("A:L"), # removed column "LABEL"
    col_types = c("text",rep("numeric",3),"text","numeric","numeric","skip",rep("numeric",4)))) %>%
    filter(
        TYPE %in% c(1,5), # filter by 1 and 5: comparative tows
        DUR >= 20 # retain tows of duration > 20min
    ) %>%
    transmute(mission = as.character(MISSION), 
              gear = as.integer(GEAR), # all is gear 9
              type = as.integer(TYPE),
              setno = as.integer(SETNO),
              strat = as.character(STRAT),
              depth = DMIN, # depth
              lon = - (floor(SLONG/100) + SLONG%%100/60), 
              lat = floor(SLAT/100) + SLAT%%100/60) %>%
    mutate(vessel = as.factor(substr(mission, 1, 3)))
              

# ----------------------------------------------------------
# pair tows by distance/lat-long

# For each TEL comparative tow, find the nearest NED tow
data.TEL <- data.set %>% 
    filter(vessel == "TEL", type == 5) %>%
    mutate(X=lon, Y=lat) %>% 
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", 19) %>%
    PBSmapping::convUL()
data.NED <- data.set %>% 
    filter(vessel == "NED") %>%
    mutate(X=lon, Y=lat) %>% 
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", 19) %>%
    PBSmapping::convUL()

# pair: within 1km
tmp.dist <- RANN::nn2(data.NED[, c("X","Y")], data.TEL[c("X","Y")], 1) 

data.paired <- bind_rows(
    data.TEL %>%
        mutate(station = row_number(),
               dist = tmp.dist$nn.dists[,1],
               trip = factor(mission, labels = c(1:4), levels = c("TEL2005545","TEL2005605","TEL2005633","TEL2006614"))),
    data.NED[tmp.dist$nn.idx,] %>% 
        mutate(station = row_number(),
               dist = tmp.dist$nn.dists[,1],
               trip = factor(mission, labels = c(1:4), levels = c("NED2005001","NED2005027","NED2005034","NED2006001"))))%>%
    mutate(filter.dist = dist < 2) %>% # side by side within 2km
    mutate(trip = as.integer(trip)) %>% 
    group_by(station) %>%
    mutate(filter.trip = all(trip == mean(trip))) %>%
    ungroup # trips are paired

# check setno: one special trip
data.paired %>% 
    filter(filter.dist,filter.trip) %>%
    transmute(station, mission, setno) %>%
    spread(mission, setno) %>%
    print(n = Inf)

# check strat: four pairs have different strat
data.paired %>% 
    filter(filter.dist,filter.trip) %>%
    transmute(station, vessel, strat) %>%
    spread(vessel, strat) %>%
    filter(NED!=TEL) %>%
    left_join(data.paired) %>% 
    select(-gear,-type,-filter.dist,-filter.trip,-NED,-TEL)

# chech map: different colors for various filtering conditions
library(maps)
library(ggplot2)
base_map <- ggplot() +
    borders("world", colour="gray50", fill = "gray90", xlim = c(-70,-56), ylim = c(40, 48)) +
    coord_fixed(ratio = 1.2, xlim = c(-70,-56), ylim = c(40, 48)) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        legend.position = "bottom"
    )

base_map +
    geom_point(data=data.set, aes(lon, lat, color = mission)) +
    theme_bw()

base_map +
    geom_point(data=data.set, aes(lon, lat, color = vessel), shape = 21) +
    geom_point(data=data.paired %>% filter(filter.dist,filter.trip), aes(lon, lat), shape = 19, size = 0.2) +
    geom_point(data=data.paired %>% filter(!filter.dist), aes(lon, lat), shape = 3, color = "blue") +
    geom_point(data=data.paired %>% filter(!filter.trip), aes(lon, lat), shape = 4, color = "yellow") +
    theme_bw()
ggsave(filename = "data_description/NED2005/distance-paired_NED2005_map.jpg", width = 12, height = 8, dpi = 600)


# ---------------------------------------------------
# Details table: catch at lengh
# need to confirm: sum together duplicating rows?

path.detail <- paste0(data.path,"2005_2006 Comparative Details (All Sets).xlsx")
sheet.detail <- excel_sheets(path.detail)[1:8]
data.detail <- map_df(
    sheet.detail, ~ read_excel(
        path.detail, sheet = .x, 
        col_names = c("mission","setno","species","size","len","catch"),
        range = cell_limits(c(2, 1), c(NA, 7)),
        col_types = c("text", "numeric", "skip", "numeric", "text", "numeric", "numeric"))) %>%
    select(-size) %>% # need to confirm if size matters, though size = 1 for all
    transmute(mission = as.character(mission),
              setno = as.integer(setno),
              species = as.integer(species),
              len = as.integer(len),
              catch = catch) %>%
    group_by(mission, setno, species, len) %>%
    summarise(catch = sum(catch)) %>% # sum together duplicating rows
    ungroup()


# ---------------------------------------------------
# Catch table:catch by tow, here only for extracting species names

path.catch <- paste0(data.path,"2005_2006 Comparative Catch (All Sets).xlsx")
sheet.catch <- excel_sheets(path.catch)[1:8]
data.catch <- map_df(
    sheet.catch, ~ read_excel(
        path.catch, sheet = .x, 
        col_names = c("species","name"),
        range = cell_limits(c(2,6),c(NA,7)),
        col_types = c("numeric", "text"))) %>%
    distinct(species, name) %>%
    mutate(species = as.integer(species))


# ---------------------------------------------------
# combine tables
d.length <- data.paired %>%
    filter(filter.dist,filter.trip) %>%
    select(mission, setno, station, depth, lon, lat, vessel) %>%
    inner_join(data.detail, by = c("mission", "setno")) %>%
    left_join(data.catch, by = "species") %>%
    as.data.frame()


save(file = "data-NED2005.RData", d.length)


