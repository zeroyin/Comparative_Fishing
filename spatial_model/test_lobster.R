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

m2deg <- function(x) x %>%
  `attr<-`("projection","UTM") %>% 
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
data.survey <- read.csv("data/AllSurveyR1new.csv", header = T)
grid.infor = read.csv("data/34_38grids.dat",sep=',')
grid.nb = read.csv("data/adjacency_matrix.dat",sep=',')

# mesh and SPDE: extend Adam's grid 
loc <- grid.infor %>% transmute(X=(x1+x2+x3+x4)/4,Y=(y1+y2+y3+y4)/4) %>% deg2m()
bnd <- inla.nonconvex.hull(as.matrix(loc), convex = 5)
mesh <- inla.mesh.create(loc = loc, boundary = bnd, refine = F, extend = F)
spde <- inla.spde2.matern(mesh, alpha = 2)

# mesh grid points as knots
d.knots <- cbind(1:mesh$n,mesh$loc[,1:2]) %>%
  `colnames<-`(c("knotno","X","Y")) %>%
  as.data.frame() %>%
  mutate(bnd = !(knotno %in% mesh$idx$loc)) %>%
  m2deg() 

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
  ungroup()  %>%
  mutate(date = as.character(date)) %>% # convert data type for date
  filter(lubridate::month(date) %in% c(6,7,8)) %>% # select summer surveys for model
  filter(vessel %in% c("CapeRosewayBalloon","CapeRosewayNest","JosiesPrideNest","Needler")) %>% # select these vessels for model
  filter(year >= 2014 & year <=2018) %>% # select years for model
  filter(!is.na(knotno)) # select data within the study area


update_grid = T
n_knots <- 250

# --------------------------------------
# update with a data-dependent grid
# aggregation by kmeans is random and not compatible with mesh cutoff
# aggregation by hclust can avoid knots too closeby using large nknots

if(update_grid){
  
  # make sure paried tows are assigned same knot
  knots.tow <- d %>% group_by(siteno) %>% summarise(X = mean(X), Y = mean(Y)) %>% ungroup
  
  # generate preliminary knots: kmeans and knn
  knots.fit <- knots.tow %>% select(X, Y) %>% deg2m() %>% stats::kmeans(centers = n_knots)
  
  # re-run for improved aggregation
  for(i.iter in 1:50){
    print(knots.fit$tot.withinss)
    knots.fit1 <- knots.tow %>% select(X, Y) %>% deg2m() %>% stats::kmeans(centers = n_knots)
    if(knots.fit1$tot.withinss < knots.fit$tot.withinss) knots.fit <- knots.fit1
  }
  
  # check distance between knots
  summary(RANN::nn2(data = knots.fit$centers,k = 2)$nn.dist[,2])
  
  # mesh and SPDE
  mesh <- inla.mesh.create(loc = knots.fit$centers, boundary = bnd, extend = F, refine = F)
  spde <- inla.spde2.matern(mesh, alpha = 2)
  
  # mesh grid points as knots
  d.knots <- cbind(1:mesh$n,mesh$loc[,1:2]) %>%
    `colnames<-`(c("knotno","X","Y")) %>%
    as.data.frame() %>%
    mutate(bnd = !(knotno %in% mesh$idx$loc)) %>%
    m2deg() 
  
  
  # attach updated knots
  d <- d %>%
    select(-knotno) %>% # remove existing knots
    left_join( # left_join: make sure order does not change
      knots.tow %>% mutate(knotno = RANN::nn2(
        # in case some mesh knots were NAs due to cutoff distance
        data = d.knots %>% select(X, Y) %>% deg2m(),
        query = knots.tow %>% select(X, Y) %>% deg2m(),
        k = 1)$nn.idx[,1]) %>%
        select(siteno, knotno),
      by = "siteno"
    )
  
  # check aggregation on map
  p <- basemap +
    geom_point(data = d, aes(X,Y,color = factor(knotno))) +
    geom_point(data = d.knots %>% filter(!bnd), aes(X,Y), shape = "+", size = 4) +
    geom_point(data = d.knots %>% filter(bnd), aes(X,Y), shape = "*", size = 1) +
    theme(legend.position = "none")
  ggsave(filename = paste0("knots_map-mesh_",n_knots,".jpg"), plot = p, width = 6, height = 8, units = "in")
  
}
# --------------------------------------

# # plot sample conversion
# p <- d %>%
#   mutate(c=LobCatch/AreaSwept) %>%
#   complete(knotno, year, vessel, fill=list(c=NA)) %>%
#   group_by(knotno, vessel, year) %>%
#   summarise(ave = mean(c,na.rm=T),sd = sd(c,na.rm=T)) %>%
#   ungroup() %>%
#   group_by(knotno)%>%
#   mutate(n.ave=sum(!is.na(ave)),m.ave=mean(ave,na.rm=T)) %>%
#   mutate(ave.order=if_else(n.ave>2,m.ave,m.ave*NA))%>%
#   ungroup()%>%
#   mutate(p=ave/ave.order,x=as.numeric(factor(ave.order))) %>%
#   ggplot(aes(x = x, y = ave, color = vessel)) +
#   geom_point(pch=19) +
#   geom_smooth(method = "lm", formula = y~x,se = F, size=0.1)+
#   scale_y_log10() +
#   facet_wrap(~year,nrow = 1,scales="free_y") +
#   theme_bw()+
#   theme(legend.position = "bottom")
# ggsave(filename = paste0("samp_catch_ordered-mesh_",n_knots,".jpg"), plot = p, width = 12, height = 8, units = "in")
# 
# # Plot catch on map by year
# p <- basemap +
#   geom_point(data = d, aes(X,Y,color = AreaSwept)) +
#   scale_color_gradientn(colours = plot3D::jet.col(3), na.value = "gray", trans = "log10") +
#   facet_grid(vessel~year) +
#   theme_bw()
# ggsave(filename = "lobcatch_map.jpg", plot = p, width = 12, height = 15, units = "in")
# 
# # plot Adam's grid
# plot(d$X, d$Y, pch = ".", col = "blue",
#      xlim = range(grid.infor$x1), ylim = range(grid.infor$y1))
# for(i in 1:nrow(grid.infor)){
#   with(grid.infor[i,], polygon(c(x1,x4,x3,x2),c(y1,y4,y3,y2)), col = "gray")
# }
# 
# # plot tow location distribution and knots
# plot(d %>% select(X, Y) %>% deg2m(),
#      xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]),
#      main = paste0("Dimension=",mesh$n), pch = ".", col = "blue")
# # points(mesh$loc)
# text(mesh$loc)
# plot(mesh, add = T)


# --------------------------------------
# Make TMB input
# --------------------------------------
library(TMB)

# model options:
# (1) ST effect temporal dependece. 0: independent; 1: RW; 2: AR(1)
# (2) catch by site: 0: negative binomial; 1: zero-inflated NB
# (3) Paired catch. 0: binomial; 1: beta-binomial
# (4) Depth effect. 0: fixed, covariate; 1

tmb.data = list(
  nobs = n_distinct(d$EID),
  nsite = n_distinct(d$siteno),
  nyear = n_distinct(d$year),
  nvessel = n_distinct(d$vessel),
  nsurvey = n_distinct(d$survey),
  knot = d$knotno-1, # spatial knot number
  site = as.integer(as.factor(d$siteno)) - 1, # site number
  year = as.integer(as.factor(d$year)) - 1, # year 
  vessel = as.integer(as.factor(d$vessel)) - 1, # vessel number
  survey = as.integer(as.factor(d$survey)) -1, # survey number
  offset = log(d$AreaSwept), # offset for tow standardization
  C = d$LobCatch, # catch number
  spde = spde$param.inla[c("M0","M1","M2")], # spatial struct
  Options = c(1) # model options
)

tmb.param = list(
  log_rho = rep(0, tmb.data$nvessel-1), # relative catchability
  log_N = rep(1, tmb.data$nyear), # total abundance by year
  log_S = rep(0, mesh$n), # spatial effect
  log_ST = matrix(0, nrow = mesh$n, ncol = tmb.data$nyear), # spatiotemporal effect 
  log_kappa_s = 0, # range for S
  log_tau_s = 0,  # variance for S
  log_kappa_st = 0, # range for ST
  log_tau_st = 0, # variance for ST
  beta = 0, # temporal correlation for ST
  log_nbk = 0, # NB over-dispersion parameter
  log_phi = 0 # BB over-dispersion parameter
)

tmb.map <- list()
if(tmb.data$Options[1] != 2) tmb.map$beta <- factor(NA)


dyn.unload("test_sp_model_3")
compile("test_sp_model_3.cpp")
dyn.load("test_sp_model_3")

obj <- MakeADFun(
  tmb.data, tmb.param, tmb.map,
  random = c("log_S", "log_ST"),
  DLL = "test_sp_model_3",
  inner.control = list(trace = FALSE)
)

system.time(
  opt<-nlminb(obj$par,obj$fn,obj$gr)
) 



# --------------------------------------
# Results
# --------------------------------------

# estimated conversion
est.rho <- data.frame(
  vessel = levels(as.factor(d$vessel)),
  rho = c(1, exp(obj$report()$log_rho)),
  stringsAsFactors = F
)

x <- c(1, exp(obj$report()$log_rho))
names(x) <- levels(as.factor(d$vessel))
x


rep <- obj$report()


# predicted spatial effect 

pred.dens <- cbind(
  1:mesh$n,
  exp(rep(1, nrow(rep$log_ST))%o%rep$log_N+rep$log_S%o%rep(1, ncol(rep$log_ST))+rep$log_ST)
) %>%
  `colnames<-`(c("knotno",as.numeric(levels(factor(d$year))))) %>%
  as.data.frame() %>%
  gather(year, dens, -knotno) %>%
  mutate(year = as.integer(year)) %>%
  inner_join(d.knots, by = "knotno")

obs.dens <- d %>% 
  inner_join(est.rho, by = "vessel") %>%
  mutate(catch=LobCatch/AreaSwept/rho) %>%
  group_by(knotno, year) %>%
  summarise(dens = mean(catch,na.rm=T)) %>%
  ungroup() %>%
  complete(knotno=d.knots$knotno,year=unique(d$year),fill = list(catch=NA)) %>%
  inner_join(d.knots,by = "knotno")


p <- basemap +
  geom_point(data = bind_rows(
    pred.dens %>% mutate(type="pred"),
    obs.dens %>% mutate(type= "obs")
  ) %>% filter(!bnd),
  aes(X,Y, fill = log(dens+1)),pch=21,color="gray") +
  scale_fill_gradientn(colors = plot3D::jet.col(5),na.value = "white") +
  facet_grid(type~year)
ggsave(filename = paste0("pred+obs+dens_mesh_","reg",".jpg"), plot = p, width = 12, height = 8, units = "in")


# residuals
res <- d %>% mutate(
  offset = tmb.data$offset,
  obs = LobCatch,
  pred = obj$report()$mu,
  resid = d$LobCatch - obj$report()$mu,
  log_resid = log(abs(resid))*sign(resid)
)

# residual hist
hist(res$log_resid, main = paste0("sd=", round(sd(res$resid), 2)))


# residual map
basemap +
  geom_point(data = res, aes(X,Y,color = log_resid), size = 0.5) +
  scale_color_gradientn(colours = plot3D::jet.col(3), na.value = "gray") +
  facet_grid(vessel~year) +
  theme_bw()


# residual scatter: obs-pred
ggplot(res) +
  geom_point(aes(obs, pred), color = "blue", size = 0.2) +
  facet_grid(vessel~year, scales = "free_y") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "gray")


# residual scatter: pred-resid
ggplot(res) +
  geom_point(aes(log(pred), resid), color = "blue", size = 0.2) +
  facet_grid(vessel~year, scales = "free_y") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 0, color = "gray")


#am(parse_date_time(date, orders = "%d%m%y %H%M"))
ggplot(res) +
  geom_point(aes(log(Z), log_resid)) + 
  facet_grid(vessel~year) +
  theme_bw()


