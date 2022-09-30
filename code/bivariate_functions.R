library(INLA)


get_data <- function(df, tp, i){
  obs1 <- filter(df, t == tp & intersection == i)
  obs2 <- filter(df, t == tp-1 & !is.na(xmap) & intersection == i)
  obs3 <- filter(df, t == tp+1 & !is.na(xmap) & intersection == i)
  obs4 <- rbind(obs1, obs2, obs3)
  return(obs4)
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

#val2 <- data_inla(f1, 3, polyg1, w1, n1, 1)


sim_inla <- function(f_new, tp_new, poly_new, w_new, n_new, int_new){
  val2 <- data_inla(f_new, tp_new, poly_new, w_new, n_new, int_new)
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

sim_inla(f_sim1, 3, polyg1, w1, n1, 1) 



all_pred <- NULL
pred_cv <- NULL
all_inla <- function(df, df2, int_list, int_df){
  for (r in 1:7){
    for (i in 1:4){
      tryCatch({
        pred = sim_inla(df,r,int_list,df2,int_df, i)
        pred2 = data.frame(pred)
        df_p = data.frame(pred2, r)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(pred_cv)
}


# Results-----------------------------------------------------------------------

guess1 <- all_inla(f_sim1, w1, polyg1, n1)

x_error1 <- sqrt(sum((guess1$x_actual-guess1$xpred)^2)/nrow(guess1)) #1.541
y_error1 <- sqrt(sum((guess1$y_actual-guess1$ypred)^2)/nrow(guess1)) #1.641

c_x1 <- guess1 %>% group_by(clust) %>% summarise(x = sqrt((sum((xpred-x_actual)^2)/n()))) #1- 1.42, 2- 1.63
c_y1 <- guess1 %>% group_by(clust) %>% summarise(y = sqrt((sum((ypred-y_actual)^2)/n()))) #1- 1.70, 2- 1.59

guess2 <- all_inla(f_sim2, w2, polyg2, n2)

x_error2 <- sqrt(sum((guess2$x_actual-guess2$xpred)^2)/nrow(guess2)) #1.627
y_error2 <- sqrt(sum((guess2$y_actual-guess2$ypred)^2)/nrow(guess2)) #1.58

c_x2 <- guess2 %>% group_by(clust) %>% summarise(x = sqrt((sum((xpred-x_actual)^2)/n()))) #1- 1.80, 2- 1.57
c_y2 <- guess2 %>% group_by(clust) %>% summarise(y = sqrt((sum((ypred-y_actual)^2)/n()))) #1- 1.72, 2- 1.54


guess3 <- all_inla(f_sim3, w3, polyg3, n3)

x_error3 <- sqrt(sum((guess3$x_actual-guess3$xpred)^2)/nrow(guess3)) #1.431
y_error3 <- sqrt(sum((guess3$y_actual-guess3$ypred)^2)/nrow(guess3)) #1.514

c_x3 <- guess3 %>% group_by(clust) %>% summarise(x = sqrt((sum((xpred-x_actual)^2)/n()))) #1- 1.49, 2- 1.36
c_y3 <- guess3 %>% group_by(clust) %>% summarise(y = sqrt((sum((ypred-y_actual)^2)/n()))) #1- 1.51, 2- 1.51




c_x1_u <- xy_join1 %>% group_by(clust) %>% summarise(x = sqrt((sum((xmap_pred-xmap)^2)/n()))) #1- 1.41, 2- 1.57
c_y1_u <- xy_join1 %>% group_by(clust) %>% summarise(y = sqrt((sum((ymap_pred-ymap)^2)/n()))) #1- 1.58, 2- 1.47

c_x2_u <- xy_join2 %>% group_by(clust) %>% summarise(x = sqrt((sum((xmap_pred-xmap)^2)/n()))) #1- 1.70, 2- 1.50
c_y2_u <- xy_join2 %>% group_by(clust) %>% summarise(y = sqrt((sum((ymap_pred-ymap)^2)/n()))) #1- 1.67, 2- 1.52

c_x3_u <- xy_join3 %>% group_by(clust) %>% summarise(x = sqrt((sum((xmap_pred-xmap)^2)/n()))) #1- 1.46, 2- 1.34
c_y3_u <- xy_join3 %>% group_by(clust) %>% summarise(y = sqrt((sum((ymap_pred-ymap)^2)/n()))) #1- 1.48, 2- 1.30
