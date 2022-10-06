get_data <- function(df, tp, i){
  obs1 <- filter(df, t == tp & intersection == i)
  obs2 <- filter(df, t == tp-1 & !is.na(xmap) & intersection == i)
  obs3 <- filter(df, t == tp+1 & !is.na(xmap) & intersection == i)
  obs4 <- rbind(obs1, obs2, obs3)
  return(obs4)
}

tt <- filter(w_new, intersection == 5)
tt2 <- filter(tt , t == 55)

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


back_toget_new <- filter(back_toget, imputed == FALSE)
back_toget_new$cv <- FALSE
back_toget_new$cv[back_toget_new$missing2 == "yes"] <- TRUE

grd <- grid_tp(anim_plot_data1, 55, m3)

#f = anim_plot_data1
#w = back_toget

data_inla_ice(anim_plot_data1, 55, m3, back_toget, m1, 5)

data_inla_ice <- function(f, tp, poly, w, n, int){
  
  grd <- grid_tp(f, tp, poly)
  
  back_toget_new <- filter(w, imputed == FALSE)
  back_toget_new$cv <- FALSE
  back_toget_new$cv[back_toget_new$missing2 == "yes"] <- TRUE
  
  final_wgrd <- inner_join(back_toget_new, grd, by = "gpid")[,c(1,2,3,4,5,9,14,15)]
  
  colnames(final_wgrd) <- c("gpid","t", "cluster", "xmap", "ymap", "cv", "x_grid", "y_grid")
  
  withint <- inner_join(final_wgrd, n[,c(1,5)], by = "gpid")
  
  withint2<- withint %>% distinct(gpid, t, intersection, .keep_all= TRUE)
  
  with_act <- inner_join(withint2, f[,c(1,2,5,6,9)], by = c("gpid", "t"))
  colnames(with_act) <- c("gpid","t", "cluster1", "xmap", "ymap", "cv", "x_grid", "y_grid", "intersection", "x_actual", "y_actual", "cluster2")
  
  w1_known <- filter(with_act , cv == FALSE)
  w1_miss <- filter(with_act , cv == TRUE)
  
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

ice_inla(anim_plot_data3, 64, m4, back_toget3, m2, 10)

ice_inla <- function(f_new, tp_new, poly_new, w_new, n_new, int_new){
  val2 <- data_inla_ice(f_new, tp_new, poly_new, w_new, n_new, int_new)
  n <- nrow(val2) #Seems to be doing better when have more data 
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
  
  form2 <- R ~ I + f(idx.x, xval, model="iid", n=2) + f(idx.y, yval, model="iid", n=2) + f(idx.t, t, model="rw1", n=2)
  
  
  m0.2lik <- inla(form2, family = c("gaussian", "gaussian"),
                  data = data2,
                  control.predictor = list(compute = TRUE),
                  control.compute=list(dic=TRUE)#, verbose = TRUE
  )
  
  x.na <- which(is.na(val2$xmap))
  y.na <- which(is.na(val2$ymap))
  
  r1 <- m0.2lik$summary.fitted.values[x.na, c("mean", "sd")]
  r2 <- m0.2lik$summary.fitted.values[n+ y.na, c("mean", "sd")]
  
  
  m <- val2[x.na,]
  m$xpred <- r1[,1]
  m$x_sd <- r1[,2]
  m$ypred <- r2[,1]
  m$y_sd <- r2[,2]
  
  
  #m <- m[,c(1,9,15,16,10,17,18,11,12)] #giving us grid values as prediction
  
  return(m)
}

all_pred <- NULL
pred_cv <- NULL
all_inla_ice <- function(df, df2, int_list, int_df){
  for (r in min(df$t):max(df$t)){
   for (i in 1:max(int_df$intersection)){
      tryCatch({
        pred = ice_inla(df,r,int_list,df2,int_df, i)
        pred2 = data.frame(pred)
        df_p = data.frame(pred2, r)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(pred_cv)
}

82112

## Week 1 ----------------------------------------------------------------------
test1 <- all_inla_ice(anim_plot_data1, back_toget, m3,m1)

x_error1 <- sqrt(sum((test1$x_actual-test1$xpred)^2)/nrow(test1)) #3.162
y_error1 <- sqrt(sum((test1$y_actual-test1$ypred)^2)/nrow(test1)) #7.307

c_x1 <- test1 %>% group_by(cluster1) %>% summarise(x = sqrt((sum((xpred-x_actual)^2)/n()))) #1- 3.15, 2- 4.58, 3 - 2.85, 4 - 2.97, 5 - 4.00, 6 - 3.22
c_y1 <- test1 %>% group_by(cluster1) %>% summarise(y = sqrt((sum((ypred-y_actual)^2)/n()))) #1- 13.8, 2- 4.64, 3 - 9.83, 4 - 2.79, 5 - 8.81, 6 - 3.76

#Same pattern as with univariate. Generally performing worse, but able to get an estimate for each missing datapoint

## Week 3 ----------------------------------------------------------------------

test3 <- all_inla_ice(anim_plot_data3, back_toget3, m4,m2)

x_error3 <- sqrt(sum((test3$x_actual-test3$xpred)^2)/nrow(test3)) #3.131
y_error3 <- sqrt(sum((test3$y_actual-test3$ypred)^2)/nrow(test3)) #3.04

c_x3 <- test3 %>% group_by(cluster1) %>% summarise(x = sqrt((sum((xpred-x_actual)^2)/n()))) #1- 2.94, 2- 2.68, 3 - 3.70, 4 - 2.61, 5 - 3.11, 6 - 3.25
c_y3 <- test3 %>% group_by(cluster1) %>% summarise(y = sqrt((sum((ypred-y_actual)^2)/n()))) #1- 2.85, 2- 2.97, 3 - 3.61, 4 - 2.87, 5 - 2.99, 6 - 3.03

##Week 2 -----------------------------------------------------------------------


test2 <- all_inla_ice(anim_plot_data2, back_toget2, m6,m5)

x_error2 <- sqrt(sum((test2$x_actual-test2$xpred)^2)/nrow(test2)) #3.190
y_error2 <- sqrt(sum((test2$y_actual-test2$ypred)^2)/nrow(test2)) #3.104

c_x2 <- test2 %>% group_by(cluster1) %>% summarise(x = sqrt((sum((xpred-x_actual)^2)/n()))) #1- 3.06, 2- 3.04, 3 - 2.95, 4 - 3.38, 5 - 4.32
c_y2 <- test2 %>% group_by(cluster1) %>% summarise(y = sqrt((sum((ypred-y_actual)^2)/n()))) #1- 4.51, 2- 2.88, 3 - 2.90, 4 - 3.17, 5 - 4.28
