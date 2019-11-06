rm(list=ls())
setwd("C:/Users/yinyi/Dropbox/BIO/Comparative_Fishing/Workspace/spatial_model/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(INLA)

# function to transform to UTM from LL
deg2m <- function(x) x %>%
  `attr<-`("projection","LL") %>% 
  `attr<-`("zone","19") %>%
  `attr<-`("datum","WGS84") %>%
  PBSmapping::convUL(southern = F) 

# base map for ggplot
basemap <- ggplot() +
  borders("world", colour="gray50", fill = "gray90", xlim = c(-68,-64), ylim = c(42,46)) +
  coord_map(xlim = c(-68,-64), ylim = c(42.5,45.8)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        legend.position = "bottom")

# --------------------------------------
# Organize data
# --------------------------------------

# load data and grid information
grid.infor = read.csv("data/34_38grids.dat",sep=',')
grid.nb = read.csv("data/adjacency_matrix.dat",sep=',')
data.survey <- read.csv("data/AllSurveyR1new.csv", header = T)


# mesh and SPDE: extend Adam's grid 
loc <- grid.infor %>% transmute(X=(x1+x2+x3+x4)/4,Y=(y1+y2+y3+y4)/4) %>% deg2m()
bnd <- inla.nonconvex.hull(as.matrix(loc), convex = 15)
mesh <- inla.mesh.create(loc = loc, boundary = bnd, refine = F, extend = F)
spde <- inla.spde2.matern(mesh, alpha = 1)

# assign tows to a knot: use Adam's grid cells
for(i in 1:nrow(grid.infor)){
  index <- sp::point.in.polygon(
    data.survey$X,data.survey$Y,
    t(grid.infor[i,c('x1','x2','x3','x4','x5')]),
    t(grid.infor[i,c('y1','y2','y3','y4','y5')]))
  data.survey$knotno[index==1]=spde$mesh$idx$loc[grid.infor$grid[i]]
}

# Pair 2016 comparative survey: Pair by nearest distance
d.Nest <- data.survey %>% filter(survey == "LobsterNest", year == 2016) %>% deg2m
d.Ball <- data.survey %>% filter(survey == "LobsterBalloon", year == 2016) %>% deg2m
pair.dist <- RANN::nn2(d.Nest[,c("X","Y")], d.Ball[, c("X","Y")], 1)
d.pair <- data.frame(
  EID = c(d.Ball[,"EID"],d.Nest[pair.dist$nn.idx,"EID"]), 
  pairno = 1:length(pair.dist$nn.idx))

# Organize survey data 
d <- data.survey %>% 
  filter(!is.na(knotno)) %>% # within the study area
  filter(year >= 2015 & year <=2018) %>% # select years 
  mutate(date = as.character(date)) %>% # change date data type
  filter(lubridate::month(date) %in% c(6,7,8)) %>% # select summer surveys
  left_join(d.pair, by = "EID") %>% # identify paired sets
  mutate(vessel = case_when( # identify vessel
    survey == "LobsterBalloon" ~ "CapeRosewayBalloon",
    survey == "LobsterNest" & year == 2016 ~ "CapeRosewayNest",
    survey == "LobsterNest" & year > 2016 ~ "JosiesPrideNest",
    survey == "DFOsummer" ~ "Needler",
    survey == "Scallop" ~ "Scallop",
    survey %in% c("NEFSCfall","NEFSCspring") ~ "Bigelow")) %>%
  mutate(survey = if_else(!is.na(pairno), "ComparativeLobster", as.character(survey))) %>% # identify survey
  mutate(siteno = if_else(is.na(pairno), row_number(), -pairno)) %>% # negative sites are comparative
  group_by(siteno) %>% # filter out comparative sites where sum catch is zero
  filter(!(siteno < 0 & sum(LobCatch) == 0)) %>%
  ungroup() 


# # Plot catch on map by year
# p <- basemap +
#   geom_point(data = d, aes(X,Y,color = LobCatch)) +
#   scale_color_continuous(low = "blue",high = "red",na.value = "gray", trans = "log10") +
#   facet_grid(vessel~year) +
#   theme_bw()
# ggsave(filename = "lobcatch_map.jpg", plot = p, width = 12, height = 15, units = "in")
# 
# # plot tow location distribution and knots
# plot(d$X, d$Y, pch = ".", col = "blue")
# for(i in 1:nrow(grid.infor)){
#   with(grid.infor[i,], polygon(c(x1,x4,x3,x2),c(y1,y4,y3,y2)), col = "gray")
# }
# 
# plot(mesh)

# --------------------------------------
# Make TMB input
# --------------------------------------
library(TMB)

# model options:
# (1) ST effect temporal dependece. 0: independent; 1: RW; 2: AR(1)
# (2) catch by site: 0: negative binomial; 1: zero-inflated NB
# (3) Paired catch. 0: binomial; 1: beta-binomial
# 
tmb.data = list(
  nobs = n_distinct(d$EID),
  nsite = n_distinct(d$siteno),
  nyear = n_distinct(d$year),
  nvessel = n_distinct(d$vessel),
  nsurvey = n_distinct(d$survey),
  knot = as.integer(as.factor(d$knotno)) - 1, # spatial knot number
  site = as.integer(as.factor(d$siteno)) - 1, # site number
  year = as.integer(as.factor(d$year)) - 1, # year 
  vessel = as.integer(as.factor(d$vessel)) - 1, # vessel number
  survey = as.integer(as.factor(d$survey)) -1, # survey number
  offset = log(d$AreaSwept), # offset for tow standardization
  C = d$LobCatch, # catch number
  spde = spde$param.inla[c("M0","M1","M2")], # spatial struct
  Options = c(1,1,1) # model options
)

tmb.param = list(
  log_rho = c(0,0,0,0), # relative catchability
  log_N = rep(1, tmb.data$nyear), # total abundance by year
  log_S = rep(0, mesh$n), # spatial effect
  log_ST = matrix(0, nrow = mesh$n, ncol = tmb.data$nyear), # spatiotemporal effect 
  log_kappa_s = 0, # range 
  log_tau_s = 0,  # variance
  log_kappa_st = 0, # range 
  log_tau_st = 0, # variance
  log_nbk = 10, # NB over-dispersion parameter
  log_zip = c(-10,-10,-10,-10, -1), # NB zero-inlfation probability for each survey
  log_phi = 0 # BB over-dispersion parameter
)

tmb.map <- NULL


dyn.unload("spatialmodel")
compile("spatialmodel.cpp")
dyn.load("spatialmodel")

obj <- MakeADFun(
  tmb.data, tmb.param, tmb.map,
  random = c("log_S", "log_ST"),
  DLL = "spatialmodel",
  inner.control = list(trace = FALSE)
)

system.time(
  opt<-nlminb(obj$par,obj$fn,obj$gr, control = list(iter.max=1000,eval.max=1000))
) 

exp(obj$report()$log_rho)


res <- data.frame(year = d$year,
                  vessel = d$vessel, 
                  catch = d$LobCatch,
                  mu = obj$report()$mu,
                  resid = d$LobCatch - obj$report()$mu)

ggplot(res) +
  geom_point(aes(catch, mu), color = "blue", size = 0.2) +
  facet_grid(vessel~year, scales = "free_y") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "gray")


ggplot(res) +
  geom_point(aes(mu, resid), color = "blue", size = 0.2) +
  facet_grid(vessel~year, scales = "free_y") +
  scale_x_log10() + 
  theme_bw() +
  geom_abline(intercept = 0, slope = 0, color = "gray")

