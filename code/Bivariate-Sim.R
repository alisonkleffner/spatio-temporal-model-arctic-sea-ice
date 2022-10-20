## Simulate Bivariate Data

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


final2_2 <- final_t2 %>% highlight_key(~gpid)

ggplot(final2_2, aes(x = xnew, y = ynew, group = gpid, 
                     hoverinfo = NULL,
                     color = factor(gpid %% 10))) + 
  geom_path(arrow = arrow(length = unit(1, "mm")), alpha = .5) + 
  scale_color_viridis_d() + ggtitle("Simulation 2") 

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


feat_data1_2 %>%
  ggplot(aes(x = xmap, y = ymap, group = gpid, ids = gpid,
             color = clust, fill = clust)) + 
  geom_point() +
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + ggtitle("Bivariate Simulation")