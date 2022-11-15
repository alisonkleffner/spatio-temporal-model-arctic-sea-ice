## Necessary Packages ----------------------------------------------------------

library(FNN) #Neearest neighbor function for simulation
library(expss) #for vlookup
library(reshape2) #melt
library(sf) #sf objects for spatial methods
#library(ggpubr) #used for figures with ggplot
library(tidyverse)
library(fields) #spatial stuff
library(crosstalk) #html widgets
library(factoextra) #get silhouette statistic
library(sp) #spatial methods
library(mvtnorm) #simulate multivariate normal
library(zoo) #cross validation
library(INLA) #install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(GpGp) #get covariance function for simulation

## Necessary Functions ---------------------------------------------------------

bbox_summary <- function(df) {
  # Function takes data frame with k, xmap, ymap
  # and returns a bounding box and trajectory
  df <- arrange(df, k) # sort by k
  df %>%
    summarize(
      x = mean(xmap), y = mean(ymap),
      xmin = min(xmap), xmax = max(xmap),
      ymin = min(ymap), ymax = max(ymap),
      xbox = xmax - xmin,
      ybox = ymax - ymin,
      dx = xmap[n()] - xmap[1],
      dy = ymap[n()] - ymap[1],
      angle = atan2(dy, dx),
      kmin = min(k),
      kmax = max(k),
      npts = n(),
    )
}

nn_df <- function(d, t, og_df){
  d_t <- filter(d, t==t)
  d_new<- rownames_to_column(d_t, "index")
  g <- knn(d_t[,1:2], og_df[,5:6], cl = d_t$cluster, k=1)
  indices <- attr(g, "nn.index")
  og_new <- data.frame(og_df,indices) #get 6th column with nn index
  og_new$ux<- vlookup(og_new$indices,d_new,6,1)
  og_new$uy<- vlookup(og_new$indices,d_new,7,1)
  og_new$xnew <- og_new$xnew+ og_new$ux
  og_new$ynew <- og_new$ynew+ og_new$uy
  og_new$t <- og_new$t+1
  final <- og_new[,-7]
  return(final)
}

week1 <- function(x_ini_disp, y_ini_disp, rho_x, rho_y, init_grid){
  q_t1 <- data.frame(xmap = x_ini_disp+cg$xmap, ymap = y_ini_disp+cg$ymap, t = cg$t+1, cluster = cg$cluster) #, gpid = 1:400)
  q_t2 <- data.frame(xmap = (rho_x*(q_t1$xmap-cg$x)) + q_t1$xmap, ymap = (rho_y*(q_t1$ymap-cg$y)) + q_t1$ymap, t = q_t1$t+1, cluster = cg$cluster) #, gpid = 1:400)
  q_t3 <- data.frame(xmap = (rho_x*(q_t2$xmap-q_t1$xmap)) + q_t2$xmap, ymap = (rho_y*(q_t2$ymap-q_t1$ymap)) + q_t2$ymap, t = q_t2$t+1, cluster = cg$cluster) #, gpid = 1:400)
  q_t4 <- data.frame(xmap = (rho_x*(q_t3$xmap-q_t2$xmap)) + q_t3$xmap, ymap = (rho_y*(q_t3$ymap-q_t2$ymap)) + q_t3$ymap, t = q_t3$t+1, cluster = cg$cluster) #, gpid = 1:400)
  q_t5 <- data.frame(xmap = (rho_x*(q_t4$xmap-q_t3$xmap)) + q_t4$xmap, ymap = (rho_y*(q_t4$ymap-q_t3$ymap)) + q_t4$ymap, t = q_t4$t+1, cluster = cg$cluster) #, gpid = 1:400)
  q_t6 <- data.frame(xmap = (rho_x*(q_t5$xmap-q_t4$xmap)) + q_t5$xmap, ymap = (rho_y*(q_t5$ymap-q_t4$ymap)) + q_t5$ymap, t = q_t5$t+1, cluster = cg$cluster) #, gpid = 1:400)
  
  all_w1 <- rbind(cg,q_t1,q_t2,q_t3,q_t4,q_t5,q_t6)
  return(all_w1)
}

polygon_sim <- function(tp, preds){
  feat_t1 <- filter(preds, t == tp)
  feat_t1_og <- rownames_to_column(feat_t1, "index")
  points <- st_as_sf(feat_t1, coords = c("xmap", "ymap"))
  q_smallest <- filter(feat_t1_og, clust == 1)
  ###1###
  coord <- q_smallest[,c(6,7)]
  ConvexHull = chull(coord)
  coord_1 <- data.frame(coord[ConvexHull,],1)
  my.sf.point <- st_as_sf(x = coord_1, 
                          coords = c("xmap", "ymap"))
  polys = st_sf(
    aggregate(
      my.sf.point$geometry,
      list(my.sf.point$X1),
      function(g){
        st_cast(st_combine(g),"POLYGON")
      }
    ))
  #Should be final data frame needed
  y1 <- st_cast(polys, "MULTIPOLYGON")
  polyg <- y1$geometry
  z <- st_intersects(polyg,points)
  z_idex <- data.frame(z[1])
  colnames(z_idex) <- "index"
  df <- merge(z_idex, feat_t1_og, by= "index")
  df2 <- subset(feat_t1_og,!(feat_t1_og$gpid%in%df$gpid))
  ###2###
  coord2 <- df2[,c(6,7)]
  ConvexHull2 = chull(coord2)
  coord_2 <- data.frame(coord2[ConvexHull2,],2)
  my.sf.point2 <- st_as_sf(x = coord_2, 
                           coords = c("xmap", "ymap"))
  polys2 = st_sf(
    aggregate(
      my.sf.point2$geometry,
      list(my.sf.point2$X2),
      function(g){
        st_cast(st_combine(g),"POLYGON")
      }
    ))
  #Should be final data frame needed
  y2 <- st_cast(polys2, "MULTIPOLYGON")
  polyg2 <- y2$geometry
  z2 <- st_intersects(polyg2,points)
  z_idex2 <- data.frame(z2[1])
  colnames(z_idex2) <- "index"
  df3 <- merge(z_idex2, df2, by= "index")
  df4 <- subset(df2,!(df2$gpid%in%df3$gpid))
  y = rbind(y1, y2)
  return(y)
}

grid_tp_sim <- function(frame, tp, m){
  time_og <- filter(frame, t==tp)
  pts <- st_as_sf(time_og, coords = c("xmap", "ymap")) 
  grid_50 <- st_make_grid(m, cellsize = c(5, 5)) %>% 
    st_sf(grid_id = 1:length(.))
  grid_lab <- st_centroid(grid_50) %>% cbind(st_coordinates(.))
  pts_grd <- pts %>% st_join(grid_50, join = st_intersects) %>% as.data.frame
  all_pts_grd <- left_join(pts_grd,grid_lab,by= "grid_id")
  all_pts_grd2 <- all_pts_grd %>% distinct(gpid, .keep_all= TRUE)
  g1 <- all_pts_grd2[,c(1,6,8,7,9,10,11)]
  g1 <- rownames_to_column(g1, "index")
  return(g1)
}

get_data <- function(df, tp, i){
  obs1 <- filter(df, t == tp & intersection == i)
  obs2 <- filter(df, t == tp-1 & !is.na(xmap) & intersection == i)
  obs3 <- filter(df, t == tp+1 & !is.na(xmap) & intersection == i)
  obs4 <- rbind(obs1, obs2, obs3)
  return(obs4)
}


data_inla <- function(f, tp, poly, w, n, int){
  grd <- grid_tp_sim(f, tp, poly)
  final_wgrd <- inner_join(w, grd, by = "gpid")[,c(1,2,3,4,5,6,8,13,14)]
  colnames(final_wgrd) <- c("gpid","t", "ux", "uy", "xmap", "ymap", "clust", "x_grid", "y_grid")
  withint <- inner_join(final_wgrd, n[,c(1,9)], by = "gpid")
  withint2<- withint %>% distinct(gpid, t, .keep_all= TRUE)
  with_act <- inner_join(withint2, f[,c(1,2,5,6)], by = c("gpid", "t"))
  colnames(with_act) <- c("gpid","t", "ux", "uy", "xmap", "ymap", "clust", "x_grid", "y_grid", "intersection", "x_actual", "y_actual")
  
  w1_known <- filter(with_act , !is.na(xmap))
  w1_miss <- filter(with_act , is.na(xmap))
  
  w1_miss$xval <- w1_miss$x_grid
  w1_miss$yval <- w1_miss$y_grid
  
  w1_known$xval <- w1_known$x_actual
  w1_known$yval <- w1_known$y_actual
  
  w_new <- rbind(w1_known, w1_miss)
  
  #Just Get Day Before and Day After
  d <- get_data(w_new, tp, int)
  #remove na grid
  d2 <- filter(d, !is.na(x_grid))
  return(d2)
  
}


## Simulate Bivariate Data -----------------------------------------------------

xmap <- seq(1,30,by = 1)
ymap <- seq(1,30,by = 1)

u_grid <- as.matrix(expand.grid(xmap,ymap)) # create grid

set.seed(4)
clust <- kmeans(u_grid,2) # clusters are straight lines
cluster <- clust$cluster

cg <- data.frame(u_grid,cluster)
cg = data.frame(xmap=cg$Var1,ymap=cg$Var2,t=0) # get data frame of inital grid with cluster

cg$cluster <- 1
cg$cluster[cg$ymap > 15.5] <- 2

#Move Grid
all_q2 <- week1(10, 10, 0.75, 0.8, cg)


covparms_c2 <- c(20,20,2,0)

covparms_c1 <- c(1,10,5,0)



d_c1 <- filter(all_q2, cluster == 1)
d_c2 <- filter(all_q2, cluster == 2)

set.seed(4)
sigma_clus1 <- exponential_spacetime(covparms_c1,as.matrix(d_c1[,1:3]))
set.seed(4)
sigma_clus2 <- exponential_spacetime(covparms_c2, as.matrix(d_c2[,1:3]))


#Playing with mean and covariance function -> see if can separate. 
set.seed(44)
foo1 = rmvnorm(2, rep(0.5,ncol(sigma_clus1)),sigma_clus1) #u_x,c1
set.seed(44)
foo2 = rmvnorm(2, rep(2, ncol(sigma_clus2)),sigma_clus2) #u_x,c2



d_c1$rex <- c(t(foo1)[,1])
d_c2$rex <- c(t(foo2)[,1])

d_c1$rey <- c(t(foo1)[,2])
d_c2$rey <- c(t(foo2)[,2])


d <- rbind(d_c1, d_c2)%>% arrange(t) #Final data frame for underlying distribution


## Observed Grid

xmap2 <- seq(5,25,by = 2) #half the number of cells as underlying
ymap2 <- seq(5,25,by = 2)

og_grid <- as.matrix(expand.grid(xmap2,ymap2)) # create grid

og <- data.frame(gpid = 1:nrow(og_grid), xmap = og_grid[,1], ymap = og_grid[,2], t = 0)

og_clust <- og
og_clust$cluster <- 1
og_clust$cluster[og_clust$ymap > 15.5] <- 2

#Get Data - match observed value to nearest neighbor in underlying distribution

#From t = 0 to t = 1
d_t1 <- filter(d, t==0)
d_new1<- rownames_to_column(d_t1, "index")

g <- knn(d_t1[,1:2], og[,2:3], cl = d_t1$cluster, k=1)
indices <- attr(g, "nn.index")

og_new <- data.frame(og,indices)
og_new$ux <- vlookup(og_new$indices,d_new1,6,1)
og_new$uy <- vlookup(og_new$indices,d_new1,7,1)
og_new$xnew <- og_new$xmap+ og_new$ux
og_new$ynew <- og_new$ymap+ og_new$uy
og_new$t <- og_new$t+1

t_1 <- og_new[c(1,4,6,7,8,9)]
t_2 <- nn_df(d, 1, t_1)
t_3 <- nn_df(d, 2, t_2)
t_4 <- nn_df(d, 3, t_3)
t_5 <- nn_df(d, 4, t_4)
t_6 <- nn_df(d, 5, t_5)
t_7 <- nn_df(d, 6, t_6)
t_8 <- nn_df(d, 7, t_7)
t_9 <- nn_df(d, 8, t_8)
t_10 <- nn_df(d, 9, t_9)
t_11 <- nn_df(d, 10, t_10)
t_12 <- nn_df(d, 11, t_11)
t_13 <- nn_df(d, 12, t_12)
t_14 <- nn_df(d, 13, t_13)

final_t2 <- rbind(t_1, t_2, t_3, t_4, t_5, t_6, t_7,t_8, t_9, t_10, t_11, t_12, t_13, t_14) #all days - final simulated dataset


#Plot of Full Trajectory (need crosstalk and plotly packages if want to run)
# final2_2 <- final_t2 %>% highlight_key(~gpid)
# 
# ggplot(final2_2, aes(x = xnew, y = ynew, group = gpid, 
#                      hoverinfo = NULL,
#                      color = factor(gpid %% 10))) + 
#   geom_path(arrow = arrow(length = unit(1, "mm")), alpha = .5) + 
#   scale_color_viridis_d() + ggtitle("Simulation 2") 

## Create Clusters -------------------------------------------------------------

d_new <- data.frame(final_t2, k=1)
colnames(d_new) <- c("gpid", "t", "ux", "uy", "xmap", "ymap", "k")


sim_bb <- d_new %>% 
  tidyr::nest(data = -gpid) %>% 
  mutate(summary = map(data, bbox_summary)) %>%
  unnest(summary)


data_scale = data.frame(scale(sim_bb[,3:13]))
names(data_scale) <- c("x_s1", "y_s1", "xmin_s1", "xmax_s1", "ymin_s1", "ymax_s1", "xbox_s1", "ybox_s1", "dx_s1", "dy_s1", "angle_s1")
sim_bb_feat <- data.frame(sim_bb, data_scale)

set.seed(4)
km= kmeans(select(sim_bb_feat, c(x_s1,y_s1:angle_s1)), 2, nstart=25)
clusters = as.factor(km$cluster)

feat1 = data.frame(sim_bb_feat, clusters)

feat_data1_2 <- d_new %>%
  left_join(select(feat1, gpid, clust = clusters)) %>%
  complete(crossing(gpid, t), fill = list(imputed = T)) %>%
  ungroup() %>%
  arrange(gpid, t) %>%
  group_by(gpid) %>%
  arrange(gpid)

f1 <- filter(feat_data1_2, t < 8)
f2 <- filter(feat_data1_2, t > 7)

y1 <- polygon_sim(7, f1)
y2 <- polygon_sim(8, f2)

## Get Intersection

intersections_pp <- st_intersection(y1, y2) %>% mutate(int_name = paste0(Group.1, "-", Group.1.1))
feat_t1 <- filter(f1, t == 7)
points <-st_as_sf(feat_t1, coords = c("xmap", "ymap"))
polyg2 <-intersections_pp$geometry # LIST OF INTERSECTION POLYGONS
z <- st_intersects(polyg2,points) #Return index of points within each intersection
FrameData <- lapply(z, function(x) as.data.frame(x))
new<- reshape2::melt(FrameData)[,-1] #gives list of index and intersection group
colnames(new) <- c("gpid", "intersection")
n2 <- merge(feat_t1, new, by = "gpid")#Data frame of all of the gpids with interactions

w2 <- f1
f_sim2 <- f1

p <- 0.1 #proportion of missing values
set.seed(5)
sel <- sample(nrow(w2), size = p*nrow(w2))
# Final Data Frame that has missing values (gpid, t, ux, uy, xmap, ymap, k, clust)
for(t in sel){
  w2[t,c(5,6)] <- NA 
}


## Start SPDE INLA Process -----------------------------------------------------

#https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html

#Just one tp and intersection as an example
val2 <- data_inla(f_sim2, 4, polyg2, w2, n2, 2) #Get Data set up properly (t=3, int = 1)
val2_known <- filter(val2, !is.na(xmap))

plot(val2_known$xmap, val2_known$ymap)


#Create Mesh - nonconvex (create around known values)
prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(100, 100))

prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(3, 3), cutoff = 0.35,
                        offset = c(-0.05, -0.05))

par(mfrow=c(1,2))
plot(prmesh1, asp = 1, main = "") #see what mesh looks like
plot(val2_known$xmap, val2_known$ymap)


val.spde <- inla.spde2.matern(mesh = prmesh1, alpha = 2) #create object for matern model
A.meuse <- inla.spde.make.A(mesh = prmesh1,
                      loc = cbind(val2_known$xval, val2_known$yval), group = val2_known$t) #inla wont run without grroup object

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = val.spde$n.spde) #a necessary list of named index vectors for the SPDE model


val.stack <- inla.stack(data  = list(x = val2_known$ux),
                        A = list(A.meuse, 1),
                        effects = list(c(s.index, list(Intercept = 1)),
                                       list(t = val2_known$t)
                        ),
                        tag = "val.data")


val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$xval, val2_pred2$yval), group = val2_pred2$t)

val.stack.pred <- inla.stack(data = list(x = NA),
                             A = list(A.pred, 1),
                             effects = list(c(s.index, list (Intercept = 1)),
                                            list(t = val2_pred2$t)
                             ),
                             tag = "val.pred")
#Join stack
join.stack <- inla.stack(val.stack, val.stack.pred)


#Fit model
form <- x ~ -1 + Intercept + f(spatial.field, model = spde) + f(t, model = "rw1")

m1 <- inla(form, data = inla.stack.data(join.stack, spde = val.spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE))

#Summary of results
summary(m1)

ndex.pred <- inla.stack.index(join.stack, "val.pred")$data

#There is one that is very off but the rest look close. 
#Much larger sd values if use grid value in A.meuse
#But still slightly larger sd values than before
spde.mn <- m1$summary.fitted.values[ndex.pred, "mean"]
spde.sd<- m1$summary.fitted.values[ndex.pred, "sd"]

## Attempt 2 -------------------------------------------------------------------

val2_known <- filter(val2, !is.na(xmap))


#Create Mesh - nonconvex (create around known values)
prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(100, 100))

prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(3, 3), cutoff = 0.35,
                        offset = c(-0.05, -0.05))

plot(prmesh1, asp = 1, main = "") #see what mesh looks like


spde <- inla.spde2.pcmatern(mesh = prmesh1, 
                            prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01


k <- max(val2_known$t)

iset <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = k)


A <- inla.spde.make.A(mesh = prmesh1,
                            loc = cbind(val2_known$xval, val2_known$yval), group = val2_known$t) 


sdat <- inla.stack(
  data = list(x = val2_known$ux), 
  A = list(A), 
  effects = list(c(iset)),
  tag = 'stdata') 



val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$xval, val2_pred2$yval), group = val2_pred2$t)
iset2 <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = max(val2_pred$t))

sdat.pred <- inla.stack(
  data = list(x = NA), 
  A = list(A.pred), 
  effects = list(c(iset2)),
  tag = 'stdata.pred') 

#Join stack
join.stack <- inla.stack(sdat, sdat.pred)


#Fit model

h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
# Model formula
formulae <- x ~ 0 + f(i, model = spde, group = i.group, 
                          control.group = list(model = 'ar1', hyper = h.spec)) 
# PC prior on the autoreg. param.
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
# Model fitting
res <- inla(formulae,  data = inla.stack.data(join.stack), 
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(join.stack)), 
            control.family = list(hyper = list(prec = prec.prior)), 
            control.fixed = list(expand.factor.strategy = 'inla'))

#Summary of results
summary(res)

ndex.pred <- join.stack$data$index$stdata.pred

#There is one that is very off but the rest look close. 
#Much larger sd values if use grid value in A.meuse
#But still slightly larger sd values than before
spde.mn2 <- res$summary.fitted.values[ndex.pred, "mean"]
spde.sd2<- res$summary.fitted.values[ndex.pred, "sd"]

#Both attempts are close to each other for most part

#Attempt 1
# 4.187057 18.387 13.389 3.898043e-20 8.389 18.3898 23.389 18.3898 8.38984

#Attempt 2
#4.1867 18.388 13.3898 2.386281e-05 8.3899 18.3899 23.389 18.389 8.389806

#Seem to be just returning things close to the grid values


## Attempt 3 -------------------------------------------------------------------

val2 <- data_inla(f_sim2, 3, polyg2, w2, n2, 2) #Get Data set up properly (t=3, int = 1)
val2_known <- filter(val2, !is.na(xmap))


#Create Mesh - nonconvex (create around known values)
prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(100, 100))

prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(3, 3), cutoff = 0.35,
                        offset = c(-0.05, -0.05))

plot(prmesh1, asp = 1, main = "") #see what mesh looks like


val.spde <- inla.spde2.matern(mesh = prmesh1, alpha = 2) #create object for matern model
A.meuse <- inla.spde.make.A(mesh = prmesh1,
                            loc = cbind(val2_known$x_grid, val2_known$y_grid)) #inla wont run without grroup object

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = val.spde$n.spde) #a necessary list of named index vectors for the SPDE model


val.stack <- inla.stack(data  = list(x = val2_known$xmap),
                        A = list(A.meuse, 1),
                        effects = list(c(s.index, list(Intercept = 1)),
                                       list(t = val2_known$t)
                        ),
                        tag = "val.data")


val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$x_grid, val2_pred2$y_grid))

val.stack.pred <- inla.stack(data = list(x = NA),
                             A = list(A.pred, 1),
                             effects = list(c(s.index, list (Intercept = 1)),
                                            list(t = val2_pred2$t)
                             ),
                             tag = "val.pred")
#Join stack
join.stack <- inla.stack(val.stack, val.stack.pred)


#Fit model
form <- x ~ -1 + Intercept + f(spatial.field, model = spde) + f(t, model = "rw1")

m2 <- inla(form, data = inla.stack.data(join.stack, spde = val.spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE))

#Summary of results
summary(m1)

ndex.pred <- inla.stack.index(join.stack, "val.pred")$data

#There is one that is very off but the rest look close. 
#Much larger sd values if use grid value in A.meuse
#But still slightly larger sd values than before
spde.mn3 <- m2$summary.fitted.values[ndex.pred, "mean"]
#5.74744 18.942 12.47619 -0.0001175  7.7235 18.787 23.7201 18.7871 9.1979

spde.sd3<- m2$summary.fitted.values[ndex.pred, "sd"]


## Attempt 4

k <- 4

iset <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = k)


A <- inla.spde.make.A(mesh = prmesh1,
                      loc = cbind(val2_known$x_grid, val2_known$y_grid), group = val2_known$t) 


sdat <- inla.stack(
  data = list(x = val2_known$xmap), 
  A = list(A), 
  effects = list(c(iset)),
  tag = 'stdata') 



val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$x_grid, val2_pred2$y_grid), group = val2_pred2$t)
iset2 <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = 3)

sdat.pred <- inla.stack(
  data = list(x = NA), 
  A = list(A.pred), 
  effects = list(c(iset2)),
  tag = 'stdata.pred') 

#Join stack
join.stack <- inla.stack(sdat, sdat.pred)


#Fit model

h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
# Model formula
formulae <- x ~ 0 + f(i, model = spde, group = i.group, 
                      control.group = list(model = 'ar1', hyper = h.spec)) 
# PC prior on the autoreg. param.
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
# Model fitting
res2 <- inla(formulae,  data = inla.stack.data(join.stack), 
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(join.stack)), 
            control.family = list(hyper = list(prec = prec.prior)), 
            control.fixed = list(expand.factor.strategy = 'inla'))
#Warning in inla.core(formula = formula, family = family, contrasts = contrasts,  :
#f(i, ...) :  There is no indices where group[]=1, this is *usually* a misspesification



ndex.pred <- join.stack$data$index$stdata.pred

#There is one that is very off but the rest look close. 
#Much larger sd values if use grid value in A.meuse
#But still slightly larger sd values than before
spde.mn4 <- res2$summary.fitted.values[ndex.pred, "mean"]
#5.793  18.853  12.507 -4.487449e-21  7.712  18.756 23.628  18.756  9.0777

spde.sd4<- res2$summary.fitted.values[ndex.pred, "sd"]

#5.793  18.853  12.507 -4.487449e-21  7.712  18.756 23.628  18.756  9.0777 (A4)

#5.74744 18.942 12.47619 -0.0001175  7.7235 18.787 23.7201 18.7871 9.1979 (A3)

#4.937 17.441 13.68 26.484 8.207 16.657 25.124 16.122 10.138

a <- c(4.937, 17.441, 13.68, 26.484, 8.207, 16.657, 25.124, 16.122, 10.138)

sqrt(sum((spde.mn4-a)^2)/9) # 8.947133
sqrt(sum((spde.mn3-a)^2)/9) # 8.947355

sqrt(sum((spde.mn-a)^2)/9) #8.92679
sqrt(sum((spde.mn2-a)^2)/9) #8.926787

grid_val <- val2_pred$x_grid

#Maybe the third one isn't doing as well because at edge of mesh? only two known points at that same grid value so maybe not encompassed well?



### Loop for Ux ----------------------------------------------------------------


spde_ux <- function(fsim, t, polyg, w, n, i){

val2 <- data_inla(fsim, t, polyg, w, n, i) #Get Data set up properly (t=3, int = 1)
val2_known <- filter(val2, !is.na(xmap))


#Create Mesh - nonconvex (create around known values)
prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(100, 100))

prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(3, 3), cutoff = 0.35,
                        offset = c(-0.05, -0.05))


val.spde <- inla.spde2.matern(mesh = prmesh1, alpha = 2) #create object for matern model
A.meuse <- inla.spde.make.A(mesh = prmesh1,
                            loc = cbind(val2_known$xval, val2_known$yval), group = val2_known$t) #inla wont run without grroup object

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = val.spde$n.spde) #a necessary list of named index vectors for the SPDE model


val.stack <- inla.stack(data  = list(x = val2_known$ux),
                        A = list(A.meuse, 1),
                        effects = list(c(s.index, list(Intercept = 1)),
                                       list(t = val2_known$t)
                        ),
                        tag = "val.data")


val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$xval, val2_pred2$yval), group = val2_pred2$t)

val.stack.pred <- inla.stack(data = list(x = NA),
                             A = list(A.pred, 1),
                             effects = list(c(s.index, list (Intercept = 1)),
                                            list(t = val2_pred2$t)
                             ),
                             tag = "val.pred")
#Join stack
join.stack <- inla.stack(val.stack, val.stack.pred)


#Fit model
form <- x ~ -1 + Intercept + f(spatial.field, model = spde) + f(t, model = "rw1")

m1 <- inla(form, data = inla.stack.data(join.stack, spde = val.spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE))

#Summary of results
summary(m1)

ndex.pred <- inla.stack.index(join.stack, "val.pred")$data

#There is one that is very off but the rest look close. 
#Much larger sd values if use grid value in A.meuse
#But still slightly larger sd values than before
spde.mn <- m1$summary.fitted.values[ndex.pred, "mean"]
spde.sd<- m1$summary.fitted.values[ndex.pred, "sd"]

#spde.df <- data.frame(spde.mn, spde.sd)

g <- val2_pred
g$uxped <- spde.mn
g$ux_sd <- spde.sd

return(g)
}

j <- spde_ux(f_sim2, 3, polyg2, w2, n2, 2)
k <- spde_ux(f_sim2, 4, polyg2, w2, n2, 2)


## loop

all_pred <- NULL
pred_cv <- NULL
all_inla_spde <- function(df, df2, int_list, int_df){
  for (r in 1:7){
    for (i in 1:4){
      tryCatch({
        pred = spde_ux(df,r,int_list,df2,int_df, i)
        pred2 = data.frame(pred)
        pred_cv = rbind.data.frame(pred_cv, pred2)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(pred_cv)
}

test1 <- all_inla_spde(f_sim2, w2, polyg2, n2)

error_spde <- sqrt(sum((test1$ux-test1$uxped)^2)/nrow(test1)) #1.083


## Second Way ------------------------------------------------------------------

#takes a lot longer

spde_ux2 <- function(fsim, t, polyg, w, n, i){
  
  val2 <- data_inla(fsim, t, polyg, w, n, i) #Get Data set up properly (t=3, int = 1)
  val2_known <- filter(val2, !is.na(xmap))
  
  
  #Create Mesh - nonconvex (create around known values)
  prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                  convex = -0.03, concave = -0.05,
                                  resolution = c(100, 100))
  
  prmesh1 <- inla.mesh.2d(boundary = prdomain,
                          max.edge = c(3, 3), cutoff = 0.35,
                          offset = c(-0.05, -0.05))
  
  plot(prmesh1, asp = 1, main = "") #see what mesh looks like
  
  
  spde <- inla.spde2.pcmatern(mesh = prmesh1, 
                              prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                              prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
  
  
  k <- max(val2_known$t)
  
  iset <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = k)
  
  
  A <- inla.spde.make.A(mesh = prmesh1,
                        loc = cbind(val2_known$xval, val2_known$yval), group = val2_known$t) 
  
  
  sdat <- inla.stack(
    data = list(x = val2_known$ux), 
    A = list(A), 
    effects = list(c(iset)),
    tag = 'stdata') 
  
  
  
  val2_pred <- filter(val2, is.na(xmap))
  val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)
  
  A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$xval, val2_pred2$yval), group = val2_pred2$t)
  iset2 <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = max(val2_pred$t))
  
  sdat.pred <- inla.stack(
    data = list(x = NA), 
    A = list(A.pred), 
    effects = list(c(iset2)),
    tag = 'stdata.pred') 
  
  #Join stack
  join.stack <- inla.stack(sdat, sdat.pred)
  
  
  #Fit model
  
  h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
  # Model formula
  formulae <- x ~ 0 + f(i, model = spde, group = i.group, 
                        control.group = list(model = 'ar1', hyper = h.spec)) 
  # PC prior on the autoreg. param.
  prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
  # Model fitting
  res <- inla(formulae,  data = inla.stack.data(join.stack), 
              control.predictor = list(compute = TRUE,
                                       A = inla.stack.A(join.stack)), 
              control.family = list(hyper = list(prec = prec.prior)), 
              control.fixed = list(expand.factor.strategy = 'inla'))
  
  #Summary of results
  summary(res)
  
  ndex.pred <- join.stack$data$index$stdata.pred
  
  #There is one that is very off but the rest look close. 
  #Much larger sd values if use grid value in A.meuse
  #But still slightly larger sd values than before
  spde.mn2 <- res$summary.fitted.values[ndex.pred, "mean"]
  spde.sd2<- res$summary.fitted.values[ndex.pred, "sd"]
  
  #spde.df <- data.frame(spde.mn, spde.sd)
  
  g <- val2_pred
  g$uxped <- spde.mn2
  g$ux_sd <- spde.sd2
  
  return(g)
}

p <- spde_ux2(f_sim2, 4, polyg2, w2, n2, 2)

all_pred <- NULL
pred_cv <- NULL
all_inla_spde2 <- function(df, df2, int_list, int_df){
  for (r in 1:7){
    for (i in 1:4){
      tryCatch({
        pred = spde_ux2(df,r,int_list,df2,int_df, i)
        pred2 = data.frame(pred)
        pred_cv = rbind.data.frame(pred_cv, pred2)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(pred_cv)
}

test2 <- all_inla_spde2(f_sim2, w2, polyg2, n2)

error_spde2 <- sqrt(sum((test2$ux-test2$uxped)^2)/nrow(test2)) #1.066 (better, but not enough for the longer calculation time?)
