# ############################################################
# For simulating realistic data
# Density: MAR Survey (2015-2018)
# ############################################################

rm(list = ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/read_data/")

library(dplyr)
library(tidyr)
library(ggplot2)


# Maritime groundfish survey data

d.inf <- readRDS("../../Data/20190801-Adam-SurveyDensity/Groundfish.gsinf.rds")
d.det <- readRDS("../../Data/20190801-Adam-SurveyDensity/Groundfish.gsdet.rds")
d.cat <- readRDS("../../Data/20190801-Adam-SurveyDensity/Groundfish.gscat.rds")

# organize data
d <- d.det %>% 
    group_by(MISSION, SETNO, SPEC, FLEN) %>%
    summarise(catch = sum(CLEN)) %>%
    ungroup() %>%
    transmute(ID = paste(MISSION, SETNO),
              species = SPEC, 
              len = FLEN,
              catch = catch) %>%
    left_join(d.cat %>%
                  filter(SAMPWGT > 0, TOTWGT > 0) %>%
                  transmute(ID = paste(MISSION, SETNO),
                            species = SPEC, 
                            scale = TOTWGT/SAMPWGT),
              by = c("ID","species")) %>%
    mutate(catch = catch * scale) %>%
    filter(!is.na(catch)) %>%
    select(-scale) %>%
    left_join(d.inf %>%
                  transmute(ID = paste(MISSION, SETNO),
                            year = format(SDATE, format = "%Y"), 
                            strat = STRAT,
                            lon = - (floor(SLONG/100) + SLONG%%100/60),
                            lat = floor(SLAT/100) + SLAT%%100/60,
                            area = AREA,
                            depth = START_DEPTH),
              by = "ID") %>%
    mutate(catch = catch/area*1000) %>%
    select(-area)

save(d, file = "MARsurvey.rda")


# divide into size classes and visualize catch on map
for(i.species in c(10,11,14,23,201,204)){

    p.length.dist <- d %>%
        filter(species == i.species) %>%
        mutate(isize = floor(3*(len - min(len))/(max(len)-min(len)))) %>%
        mutate(isize = if_else(isize == 3, 2, isize)) %>%
        mutate(isizelabel = case_when(
            isize == 0 ~ paste0(min(len)," - ", floor(min(len)+(max(len)-min(len))/3)),
            isize == 1 ~ paste0(min(len)+floor((max(len)-min(len))/3)," - ", min(len)+floor((max(len)-min(len))/3*2)),
            isize == 2 ~ paste0(min(len)+floor((max(len)-min(len))/3*2)," - ", min(len)+floor((max(len)-min(len))/3*3)))) %>%
        ggplot() +
        geom_point(aes(x = lon, y = lat, color = catch)) +
        scale_color_continuous(high = "red", low = "white", trans="log10", limits = c(1, NA)) +
        geom_polygon(data = map_data("world"), aes(long, lat, group = group), fill = "white", color = "black") +
        coord_quickmap(xlim = c(-69,-57), ylim = c(41, 47)) +
        facet_grid(year~isize) +
        theme_bw()+
        theme(legend.position="bottom", axis.title = element_blank())
    ggsave(filename = paste0("survey_MAR/length_class_dist_map",i.species,".pdf"),
           plot = p.length.dist, width = 10, height = 10, limitsize = FALSE)


    p.length.strat <- d %>%
        filter(species == i.species) %>%
        group_by(strat, len) %>%
        summarise(catch = sum(catch)) %>%
        ggplot() +
        geom_tile(aes(as.factor(strat), len, fill = catch)) +
        scale_fill_continuous(low = "white", high = "red", trans = "log10", na.value = "white", limits = c(1, NA)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90)) +
        ylim(c(0,60))
    ggsave(filename = paste0("survey_MAR/length_strata_spectrum",i.species,".jpg"),
           plot = p.length.strat, width = 10, height = 6, limitsize = FALSE)
}


# ---------------------------------------------------

# Area shape file:
shp.mar <- rgdal::readOGR("../../Data/20190801-Adam-SurveyDensity/MaritimesRegionEcosystemAssessmentStrata(2014-)/MaritimesRegionEcosystemAssessmentStrata(2014-).shp")

# plot map
ggplot() +
    geom_polygon(data = shp.mar, 
                 aes(x = long, y = lat, group = group, fill = group),
                 show.legend = F,
                 color = "black") +
    geom_point(data = d.inf %>%
                   transmute(
                       lon = - (floor(SLONG/100) + SLONG%%100/60),
                       lat = floor(SLAT/100) + SLAT%%100/60),
               aes(x = lon, y = lat),
               color = "blue", size = 0.1) +
    theme_bw()


# Make grid 1km*1km: lon, lat, area
source("functions/getGrid.R")
space <- 1000
shape <- spTransform(shp.mar, CRS("+proj=utm +zone=19 +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
myGrid <- getGrid(space, shape, proj4string(shape))
MARgrid <- data.frame(
    X = myGrid$gridData$x, Y = myGrid$gridData$y,
    lon = myGrid$gridData$lon-360, lat = myGrid$gridData$lat,
    Area_km2 = (space/1000)^2,
    strat =  data.frame(X = myGrid$gridData$x, Y = myGrid$gridData$y) %>%
        SpatialPoints() %>%
        `proj4string<-`(proj4string(shape)) %>%
        over(shape) %>% 
        .$StrataID)
save(MARgrid, file = "MARgrid_1km.rda")


ggplot(MARgrid) +
    geom_histogram(aes(strat), stat = "count", fill = "gray") +
    theme_bw() +
    ylab("Area_km^2") +
    theme(axis.text.x = element_text(angle = 90),
          panel.grid = element_blank())



