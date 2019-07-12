# #############################################
# Read data: NED2005
# Two trips paired with Don's station table
# One trip paired by mission/setno
# #############################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/read_data/")

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
# pair by new table from Don
path.stn <- paste0(data.path,"compare sets 4vwx05.xls")
excel_sheets(path.stn)

data.set.p1 <- bind_rows(
    read_excel(
        path.stn,
        sheet = "ned num",
        range = cell_limits(c(3,1),c(NA,5)),
        col_names = c("mission","station","strat","id","setno"),
        col_types = c("text","numeric","text","numeric","numeric")),
    read_excel(
        path.stn,
        sheet = "tel num",
        range = cell_limits(c(3,1),c(NA,5)),
        col_names = c("mission","strat","station","id","setno"),
        col_types = c("text","text","numeric","numeric","numeric"))
) %>%
    select(-strat,-id) %>%
    filter(!is.na(station),!is.na(mission),!is.na(setno)) %>%
    inner_join(data.set, by = c("mission","setno")) %>%
    group_by(station) %>%
    filter(n() == 2) %>%
    ungroup() %>% # filter paired station
    mutate(station = as.integer(station),
           setno = as.integer(setno))


range(data.set.p1$station)

# ----------------------------------------------------------
# Update: pair "TEL2005633","NED2005034" by setno
data.set.p2 <-
    data.set %>% 
    filter(mission %in% c("TEL2005633","NED2005034")) %>%
    group_by(setno) %>%
    filter(n() == 2) %>%
    ungroup() %>% # filter paired setno
    mutate(station = as.integer(setno + 1000)) # assign station number
    

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
# Combine the three paired trips
# Then combine tables

d.length <- bind_rows(data.set.p1, data.set.p2) %>%
    group_by(station) %>%
    mutate(strat = min(strat)) %>%
    ungroup() %>% # uniform strat for same pair
    inner_join(data.detail, by = c("mission", "setno")) %>%
    left_join(data.catch, by = "species") %>%
    filter(station!=1014) %>% # additional exclusions suggested by Don
    as.data.frame()


save(file = "data-NED2005.RData", d.length)



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

# color by mission
base_map +
    geom_point(data=data.set, aes(lon, lat, color = mission)) +
    theme_bw()
ggsave(filename = "data_description/NED2005/mission_map_NED2005.jpg", width = 12, height = 8, dpi = 600)

# successful pairs
base_map +
    geom_point(data=data.set, aes(lon, lat, color = vessel), shape = 21) +
    geom_point(data=data.set.p1, aes(lon, lat), shape = 3) +
    geom_point(data=data.set.p2, aes(lon, lat), shape = 4) +
    theme_bw() +
    ggtitle(paste0("Successful Pairs = ", n_distinct(d.length$station)))
ggsave(filename = "data_description/NED2005/paired_NED2005_map.jpg", width = 12, height = 8, dpi = 600)


# check distanceL: one pair are 4km apart
x <- bind_rows(data.set.p1, data.set.p2)

y <- x %>% 
    mutate(X=lon, Y=lat) %>% 
    `attr<-`("projection", "LL") %>%
    `attr<-`("zone", 19) %>%
    PBSmapping::convUL() %>%
    group_by(station) %>%
    summarise(dist = sqrt(diff(X)^2 + diff(Y)^2)) 

y %>%
    filter(dist >2)%>%
    left_join(x) %>% 
    select(-gear,-type)

ggplot(y) +
    geom_histogram(aes(dist), fill = "white", color = "gray") +
    theme_bw()
ggsave(filename = "data_description/NED2005/paired_NED2005_distance_hist.jpg", width = 8, height = 8, dpi = 300)


# check strat: four pairs have different strat
x %>% 
    transmute(station, vessel, strat) %>%
    spread(vessel, strat) %>%
    filter(NED!=TEL) %>%
    left_join(x) %>% 
    select(-gear,-type,-NED,-TEL)


# ----------------------------------------------------
# data description

rm(list = ls())
load("data-NED2005.RData")


# Table: number of obs by species
d.length %>% 
    count(species, name) %>%
    print(n=Inf)


# figure: catch by length by species
p <- d.length %>% 
    mutate(species = paste(species,name,sep = "-")) %>%
    ggplot() +
    geom_line(aes(len, catch, group = as.factor(station)), color = "gray") +
    stat_summary(aes(len, catch), fun.y=mean, geom="line", colour="orange", size = 1) +
    facet_grid(species~vessel, scales = "free_y") +
    theme_bw()
ggsave(filename = "data_description/NED2005/catch_compare_by_len_species-all_species.pdf", 
       plot = p, width = 10, height = 200, limitsize = FALSE)


# sample mu
p <- d.length %>% 
    mutate(species = paste(species,name,sep = "-")) %>%
    select(species, station, vessel, len, catch) %>%
    spread(as.integer(vessel), catch, fill = 0) %>%
    mutate(mu = NED/(TEL+NED)) %>%
    ggplot(aes(len, mu, group = station)) +
    geom_line(color = "gray") +
    geom_point(size = 0.1) +
    ylim(c(0,1)) +
    facet_wrap(~species, scales = "free", ncol = 1) +
    theme_bw()
ggsave(filename = "data_description/NED2005/sample_mu-all_species.pdf",
       plot = p, width = 5, height = 200, limitsize = FALSE)



# sample rho
p <- d.length %>% 
    mutate(species = paste(species,name,sep = "-")) %>%
    select(species, station, vessel, len, catch) %>%
    spread(vessel, catch, fill = 0) %>%
    mutate(ratio = NED/TEL) %>%
    ggplot(aes(len, log10(ratio), group = station)) +
    geom_line(color = "gray") +
    geom_point(size = 0.1) +
    geom_abline(intercept = 0, slope = 0, color = "orange") +
    facet_wrap(~species, scales = "free", ncol = 1) +
    theme_bw()
ggsave(filename = "data_description/NED2005/sample_log10_rho-all_species.pdf",
       plot = p, width = 5, height = 200, limitsize = FALSE)



# select species
i.species <- 10

# paired catch
p <- d.length %>% 
    filter(species == i.species) %>%
    ggplot()+
    geom_line(aes(len, catch, group = station), color = "gray") +
    stat_summary(aes(len, catch), fun.y=median, geom="line", colour="orange") +
    stat_summary(aes(len, catch), fun.y=mean, geom="line", colour="blue") +
    facet_wrap(~vessel, ncol = 2) +
    theme_bw() +
    xlim(c(0,60))
ggsave(filename = paste0("data_description/NED2005/paired_catch-species_",i.species,".pdf"),
       plot = p, width = 8, height = 6)


# paired catch by station: length spectrum
p <- d.length %>%
    filter(species == i.species) %>%
    ggplot() +
    geom_tile(aes(as.factor(station), len, fill = catch)) +
    scale_fill_continuous(low = "gray", high = "black", na.value = "white", trans = "log10") +
    facet_wrap(~vessel, nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave(filename = paste0("data_description/NED2005/paired_catch_spectrum-species_",i.species,".pdf"),
       plot = p, width = 15, height = 10)



# -----------------------------------------------
# Group by strata/depth
# -----------------------------------------------


d.length %>%
    ggplot(aes(x = strat, y = depth)) +
    geom_boxplot(width = 0.5, fill = "white", ) +
    geom_jitter(aes(color = vessel), width = 0.1, size = 0.1) +
    theme_bw()




