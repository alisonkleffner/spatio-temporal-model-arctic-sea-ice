## Necessary Packages ----------------------------------------------------------

library(knitr)
library(FNN)
library(expss)
library(reshape2)
library(sf)
library(ggpubr)
library(tidyverse)
library(fields)
library(crosstalk)
library(plotly)
library(factoextra)
library(sp)
library(mvtnorm)
library(zoo)
library(INLA)
library(GpGp)

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


#Plot of Full Trajectory
final2_2 <- final_t2 %>% highlight_key(~gpid)

ggplot(final2_2, aes(x = xnew, y = ynew, group = gpid, 
                     hoverinfo = NULL,
                     color = factor(gpid %% 10))) + 
  geom_path(arrow = arrow(length = unit(1, "mm")), alpha = .5) + 
  scale_color_viridis_d() + ggtitle("Simulation 2") 


#Cluster
d_new <- data.frame(final_t2, k=1)
colnames(d_new) <- c("gpid", "t", "ux", "uy", "xmap", "ymap", "k")

d_new2 <- filter(d_new, t < 8)

sim_bb <- d_new2 %>% 
  tidyr::nest(data = -gpid) %>% 
  mutate(summary = map(data, bbox_summary)) %>%
  unnest(summary)

data_scale = data.frame(scale(sim_bb[,3:13]))
names(data_scale) <- c("x_s1", "y_s1", "xmin_s1", "xmax_s1", "ymin_s1", "ymax_s1", "xbox_s1", "ybox_s1", "dx_s1", "dy_s1", "angle_s1")
sim_bb_feat <- data.frame(sim_bb, data_scale)

set.seed(4)
km = kmeans(select(sim_bb_feat, c(x_s1,y_s1:angle_s1)), 2, nstart=25)
clusters = as.factor(km$cluster)

feat2 = data.frame(sim_bb_feat, clusters)

feat_data1_2 <- d_new2 %>%
  left_join(select(feat2, gpid, clust = clusters)) %>%
  complete(crossing(gpid, t), fill = list(imputed = T)) %>%
  ungroup() %>%
  arrange(gpid, t) %>%
  group_by(gpid) %>%
  arrange(gpid)

#Plot of Cluster Memberships
feat_data1_2 %>%
  ggplot(aes(x = xmap, y = ymap, group = gpid, ids = gpid,
             color = clust, fill = clust)) + 
  geom_point() +
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + ggtitle("Bivariate Simulation")

## Get Intersections and Missing Data ------------------------------------------


f1 <- filter(feat_data1_2, t < 8)
f2 <- filter(feat_data1_2, t > 7)


y1 <- polygon_sim(7, f1) #Get Polygons
y2 <- polygon_sim(8, f2)

## Get Intersection
intersections_pp <- st_intersection(y1, y2) %>% mutate(int_name = paste0(Group.1, "-", Group.1.1))
feat_t1 <- filter(f1, t == 7)
#feat_t1 <- filter(feat_data1_1, t == 7)
points <-st_as_sf(feat_t1, coords = c("xmap", "ymap"))
polyg1 <-intersections_pp$geometry # LIST OF INTERSECTION POLYGONS
z <- st_intersects(polyg1,points) #Return index of points within each intersection
FrameData <- lapply(z, function(x) as.data.frame(x))
new<- melt(FrameData)[,-1] #gives list of index and intersection group
colnames(new) <- c("gpid", "intersection")
n1 <- merge(feat_t1, new, by = "gpid") #Data frame of all of the gpids with interactions

w1 <- f1
f_sim1 <- f1

p <- 0.1 #proportion of missing values
set.seed(5)
sel <- sample( nrow(w1), size = p*nrow(w1))
# Final Data Frame that has missing values (gpid, t, ux, uy, xmap, ymap, k, clust)
for(t in sel){
  w1[t,c(5,6)] <- NA #w1[t,c(3,4)] <- NA
}




## Bivariate Set up ------------------------------------------------------------

val2 <- data_inla(f1, 3, polyg1, w1, n1, 1) #Get Data set up properly (t=3, int = 1)

n <- nrow(val2)  
XY <- matrix(NA, ncol =2, nrow= n*2)
XY[1:n,1] <- val2$xmap #likelihood 1
XY[n+1:n,2]<- val2$ymap #likelihood 2
I <- matrix(NA, nrow = 2*n, ncol=2)
I[1:n, 1] <- 1
I[n+1:n,2] <- 1

t.joint <- rep(val2$t,2)
idx.t <- rep(1:2, each = n)
x.joint <- rep(val2$xval, 2)
y.joint <- rep(val2$yval, 2)
idx.x <- rep(1:2, each = n)
idx.y <- rep(1:2, each = n)

data2 <- list(R = XY, I = I, xval = x.joint, yval = y.joint, idx.x = idx.x, idx.y = idx.y, t = t.joint,  idx.t = idx.t)

#Model?
form2 <- R ~ I + f(idx.x, xval, model="iid", n=2) + f(idx.y, yval, model="iid", n=2) + f(idx.t, t, model="rw1", n=2)


m0.2lik <- inla(form2, family = c("gaussian", "gaussian"),
                data = data2,
                control.predictor = list(compute = TRUE),
                control.compute=list(dic=TRUE)#, verbose = TRUE
)

summary(m0.2lik)

x.na <- which(is.na(val2$xmap))
y.na <- which(is.na(val2$ymap))

r1 <- m0.2lik$summary.fitted.values[x.na, c("mean", "sd")]
r2 <- m0.2lik$summary.fitted.values[n+ y.na, c("mean", "sd")]


m <- val2[x.na,]
m$xpred <- r1[,1]
m$x_sd <- r1[,2]
m$ypred <- r2[,1]
m$y_sd <- r2[,2]

View(m) #Seems to just be returning grid values