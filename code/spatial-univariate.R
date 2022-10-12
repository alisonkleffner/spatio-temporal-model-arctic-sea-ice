###############################
## Try Univariate Functions 
###############################

#Simulated Data ----------------------------------------------------------------

lat_sim2 <- function(int, m, m_new, tp, frame, frame2){
  grid <- grid_tp_sim(frame2,tp,m)
  known1 <- sim_known(m_new,int,tp-1,grid,frame)
  known2 <- sim_known(m_new,int,tp,grid, frame)
  known3 <- sim_known(m_new,int,tp+1,grid, frame)
  known4 <- sim_known(m_new,int,tp+2,grid, frame)
  if(nrow(known3) ==0){known <- rbind(known1, known2, known4)}
  else(known <- rbind(known1, known2, known3))
  loc <- known[,c("xmap","ymap", "obs_time")] ## locations that I will use to fit the model with time (xmap, ymap,obs_time)
  locs<-as.matrix(loc)                                ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)
  ## intercept
  N <- nrow(known)
  Y <- known[,5]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #needs to go to number of rows - 1
  known_grd <- sim_estimate(m_new,int,tp,grid, frame)
  known_grd2 <- filter(known_grd, !is.na(gpid_id))
  locs_pred <- as.matrix(known_grd2[,c("xmap", "ymap","obs_time")]) # Would need to create obs_time when don't have them (round number?)
  colnames(locs_pred) <- c("xmap", "ymap", "obs_time")
  X_pred <- matrix(c(rep(1,nrow(locs_pred)),locs_pred[,3]), ncol = 2, byrow=FALSE)
  ncondsim <- 30
  sims <- cond_sim(fit = fit, locs_pred = locs_pred[,1:2], X_pred = X_pred, covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                   ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE, nsims = ncondsim)
  ##for prediction
  pred_grid<-predictions(fit,locs_pred[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, 
                         y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  
  b <- list(sims,pred_grid)
  return(b)
}

long_sim2 <- function(int, m, m_new, tp, frame, frame2){
  grid <- grid_tp_sim(frame2,tp,m)
  known1 <- sim_known(m_new,int,tp-1,grid,frame)
  known2 <- sim_known(m_new,int,tp,grid, frame)
  known3 <- sim_known(m_new,int,tp+1,grid, frame)
  known4 <- sim_known(m_new,int,tp+2,grid, frame)
  if(nrow(known3) ==0){known <- rbind(known1, known2, known4)}
  else(known <- rbind(known1, known2, known3))
  loc <- known[,c("xmap","ymap", "obs_time")] ## locations that I will use to fit the model with time (xmap, ymap,obs_time)
  locs<-as.matrix(loc)                                ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)
  ## intercept
  N <- nrow(known)
  Y <- known[,6]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #needs to go to number of rows - 1
  known_grd <- sim_estimate(m_new,int,tp,grid, frame)
  known_grd2 <- filter(known_grd, !is.na(gpid_id))
  locs_pred <- as.matrix(known_grd2[,c("xmap", "ymap","obs_time")]) # Would need to create obs_time when don't have them (round number?)
  colnames(locs_pred) <- c("xmap", "ymap", "obs_time")
  X_pred <- matrix(c(rep(1,nrow(locs_pred)),locs_pred[,3]), ncol = 2, byrow=FALSE)
  ncondsim <- 30
  sims <- cond_sim(fit = fit, locs_pred = locs_pred[,1:2], X_pred = X_pred, covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                   ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE, nsims = ncondsim)
  ##for prediction
  pred_grid<-predictions(fit,locs_pred[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, 
                         y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  b <- list(sims,pred_grid)
  return(b)
}

all_pred <- NULL
pred_cv <- NULL
all_lat2 <- function(df, df2, int_list, int_df){
  for (r in 1:7){
    for (i in 1:4){
      tryCatch({
        pred = lat_sim2(i,int_list,int_df,r,df, df2)
        sd_df <- sd_pred_sim(pred[[1]],pred[[2]])
        grid_t <- grid_tp_sim(df2,r,int_list)
        known2 <- sim_estimate(int_df,i,r,grid_t, df)
        known3 <- filter(known2, !is.na(gpid_id))
        df_p = data.frame(sd_df, i, r, known3$gpid)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(pred_cv)
}

all_long2 <- function(df, df2, int_list, int_df){
  for (r in 1:7){
    for (i in 1:4){
      tryCatch({
        pred = long_sim2(i,int_list,int_df,r,df, df2)
        sd_df <- sd_pred_sim(pred[[1]],pred[[2]])
        grid_t <- grid_tp_sim(df2,r,int_list)
        known2 <- sim_estimate(int_df,i,r,grid_t, df)
        known3 <- filter(known2, !is.na(gpid_id))
        df_p = data.frame(sd_df, i, r, known3$gpid)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(pred_cv)
}


sd_pred_sim <- function(sim, pred_df){
  v <- NULL
  pred_df2 <- data.frame(pred_df)
  sum_squares <- rep(0,nrow(sim))
  for(j in 1:30){
    v <- sim[,j] #fill in Known with estimated
    sum_squares <- sum_squares + (pred_df - v)^2
  }
  pred_rmse <- sqrt(1/30*sum_squares)
  pred_df2$sd <- pred_rmse
  pred_df2$lc <- pred_df2$pred_df - (2*pred_rmse)
  pred_df2$uc <- pred_df2$pred_df + (2*pred_rmse)
  return(pred_df2)
}

#### Results -------------------------------------------------------------------

# Simulation 1

int_1xn = all_lat2(w1, f_sim1,polyg1,n1)


#Gives Location for Y
int_1yn = all_long2(w1,f_sim1,polyg1,n1)

xy1n <- data.frame(int_1xn$known3.gpid, int_1xn$r, int_1xn$i, int_1xn$pred_df, int_1xn$sd, int_1xn$lc, int_1xn$uc, int_1yn$pred_df, int_1yn$sd, int_1yn$lc, int_1yn$uc)
colnames(xy1n) <- c("gpid", "t", "intersection", "xmap_pred", "x_sd", "x_lc", "x_uc", "ymap_pred", "y_sd", "y_lc", "y_uc")


xy_join1n <- right_join(f_sim1, xy1n, by = c("gpid", "t"))
xy_join1n <- xy_join1n %>% distinct(gpid,t,.keep_all= TRUE)


error_x1n <- sqrt(sum((xy_join1n$xmap-xy_join1n$xmap_pred)^2)/nrow(xy_join1n)) #1.5
error_y1n <- sqrt(sum((xy_join1n$ymap-xy_join1n$ymap_pred)^2)/nrow(xy_join1n)) #1.506


xy_join1n %>% group_by(clust) %>% summarise(x = sqrt((sum((xmap-xmap_pred)^2)/n())))
xy_join1n %>% group_by(clust) %>% summarise(y = sqrt((sum((ymap-ymap_pred)^2)/n())))

# Simulation 2

#Gives Location for X
int_2xn = all_lat2(w2, f_sim2,polyg2,n2)


#Gives Location for Y
int_2yn = all_long2(w2,f_sim2,polyg2,n2)

xy2n <- data.frame(int_2xn$known3.gpid, int_2xn$r, int_2xn$i, int_2xn$pred_df, int_2xn$sd, int_2xn$lc, int_2xn$uc, int_2yn$pred_df, int_2yn$sd, int_2yn$lc, int_2yn$uc)
colnames(xy2n) <- c("gpid", "t", "intersection", "xmap_pred", "x_sd", "x_lc", "x_uc", "ymap_pred", "y_sd", "y_lc", "y_uc")


xy_join2n <- right_join(f_sim2, xy2n, by = c("gpid", "t"))
xy_join2n <- xy_join2n %>% distinct(gpid,t,.keep_all= TRUE)


error_x2n <- sqrt(sum((xy_join2n$xmap-xy_join2n$xmap_pred)^2)/nrow(xy_join2n)) #1.535
error_y2n <- sqrt(sum((xy_join2n$ymap-xy_join2n$ymap_pred)^2)/nrow(xy_join2n)) #1.544

xy_join2n %>% group_by(clust) %>% summarise(x = sqrt((sum((xmap-xmap_pred)^2)/n())))
xy_join2n %>% group_by(clust) %>% summarise(y = sqrt((sum((ymap-ymap_pred)^2)/n())))

# Simulation 3

#Gives Location for X
int_3xn = all_lat2(w3, f_sim3,polyg3,n3)

#Gives Location for Y
int_3yn = all_long2(w3,f_sim3,polyg3,n3)

xy3n <- data.frame(int_3xn$known3.gpid, int_3xn$r, int_3xn$i, int_3xn$pred_df, int_3xn$sd, int_3xn$lc, int_3xn$uc, int_3yn$pred_df, int_3yn$sd, int_3yn$lc, int_3yn$uc)
colnames(xy3n) <- c("gpid", "t", "intersection", "xmap_pred", "x_sd", "x_lc", "x_uc", "ymap_pred", "y_sd", "y_lc", "y_uc")

xy_join3n <- right_join(f_sim3, xy3n, by = c("gpid", "t"))
xy_join3n <- xy_join3n %>% distinct(gpid,t, .keep_all= TRUE)

error_x3n <- sqrt(sum((xy_join3n$xmap-xy_join3n$xmap_pred)^2)/nrow(xy_join3n)) #1.409
error_y3n <- sqrt(sum((xy_join3n$ymap-xy_join3n$ymap_pred)^2)/nrow(xy_join3n)) #1.394

xy_join3n %>% group_by(clust) %>% summarise(x = sqrt((sum((xmap-xmap_pred)^2)/n())))
xy_join3n %>% group_by(clust) %>% summarise(y = sqrt((sum((ymap-ymap_pred)^2)/n())))

### Ice Data -------------------------------------------------------------------

xmap_10CV_pred2 <- function(m_new, m, int, tp, frame){
  grid <- grid_tp(frame, tp , m)
  obs1 <- known_data(m_new,int,tp-1,grid, frame)
  obs2 <- known_data(m_new,int,tp,grid, frame)
  obs3 <- known_data(m_new,int,tp+1,grid, frame)
  obs4 <- known_data(m_new,int,tp+2,grid, frame)
  set.seed(4)
  rand_df <- obs2[sample(nrow(obs2), size=(round(0.9*nrow(obs2)))),]
  if(nrow(obs3) ==0){obs <- rbind(obs1, rand_df, obs4)}else{obs <- rbind(obs1, rand_df, obs3)}
  loc <- obs[,c("xmap","ymap", "obs_time")]
  locs<-as.matrix(loc)
  N <- nrow(locs)
  ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)
  Y<- obs[,5]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #Silent option doesn't show the iterations. Maybe nice so get list of errors
  predicted_locs <- known_data_grid(m_new,int,tp,grid, frame)
  est_df <- subset(predicted_locs,!(predicted_locs$gpid%in%rand_df$gpid))
  pred_loc <- as.matrix(est_df[ ,c("xmap","ymap", "obs_time")])
  X_pred <- matrix(c(rep(1,nrow(pred_loc)),pred_loc[,3]), ncol = 2, byrow=FALSE)
  ncondsim <- 30
  sims <- cond_sim(fit = fit, locs_pred = pred_loc[,1:2], X_pred = X_pred, covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                   ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE, nsims = ncondsim)
  ## for prediction
  pred_grid<-predictions(fit,pred_loc[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  preds <- data.frame(est_df$gpid,pred_grid)
  b <- list(sims,preds)
  return(b)
}

ymap_10CV_pred2 <- function(m_new, m, int, tp, frame){
  grid <- grid_tp(frame, tp , m)
  obs1 <- known_data(m_new,int,tp-1,grid, frame)
  obs2 <- known_data(m_new,int,tp,grid, frame)
  obs3 <- known_data(m_new,int,tp+1,grid, frame)
  obs4 <- known_data(m_new,int,tp+2,grid, frame)
  set.seed(4)
  rand_df <- obs2[sample(nrow(obs2), size=(round(0.9*nrow(obs2)))),]
  if(nrow(obs3) ==0){obs <- rbind(obs1, rand_df, obs4)}else{obs <- rbind(obs1, rand_df, obs3)}
  loc <- obs[,c("xmap","ymap", "obs_time")]
  locs<-as.matrix(loc)
  N <- nrow(locs)
  ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)
  Y<- obs[,6]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #Silent option doesn't show the iterations. Maybe nice so get list of errors
  predicted_locs <- known_data_grid(m_new,int,tp,grid, frame)
  est_df <- subset(predicted_locs,!(predicted_locs$gpid%in%rand_df$gpid))
  pred_loc <- as.matrix(est_df[ ,c("xmap","ymap", "obs_time")])
  X_pred <- matrix(c(rep(1,nrow(pred_loc)),pred_loc[,3]), ncol = 2, byrow=FALSE)
  ncondsim <- 30
  sims <- cond_sim(fit = fit, locs_pred = pred_loc[,1:2], X_pred = X_pred, covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                   ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE, nsims = ncondsim)
  ## for prediction
  pred_grid<-predictions(fit,pred_loc[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  preds <- data.frame(est_df$gpid,pred_grid)
  b <- list(sims,preds)
  return(b)
}

sd_pred <- function(sim, pred_df){
  v <- NULL
  sum_squares <- rep(0,nrow(sim))
  for(j in 1:30){
    v <- sim[,j] #fill in Known with estimated
    sum_squares <- sum_squares + (pred_df$pred_grid - v)^2
  }
  pred_rmse <- sqrt(1/30*sum_squares)
  pred_df$sd <- pred_rmse
  pred_df$lc <- pred_df$pred_grid - pred_rmse
  pred_df$uc <- pred_df$pred_grid + pred_rmse
  return(pred_df)
}


all_cv_x_pred2 <- function(int_df,int_list,df){
  for (p in min(df$t):max(df$t)){
    for (q in 1:max(int_df$intersection)){
      tryCatch({
        pred = xmap_10CV_pred2(int_df, int_list, q, p, df)
        sd_df <- sd_pred(pred[[1]],pred[[2]])
        df_p = data.frame(sd_df, q, p)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  all_pred <- rbind.data.frame(all_pred,pred_cv)
  return(all_pred)
}

all_cv_y_pred2 <- function(int_df,int_list,df){
  for (p in min(df$t):max(df$t)){
    for (q in 1:max(int_df$intersection)){
      tryCatch({
        pred = ymap_10CV_pred2(int_df, int_list, q, p, df)
        sd_df <- sd_pred(pred[[1]],pred[[2]])
        df_p = data.frame(sd_df, q, p)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    #all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  all_pred <- rbind.data.frame(all_pred,pred_cv)
  return(all_pred)
}

xmap_10CV_pred(m1,m3, 6, 55, anim_plot_data1) #20, 9: 1,764.0655 4,926.8697 58,799,533.8787 0 
xmap_10CV_pred2(m1,m3, 6, 55, anim_plot_data1) #17, 8: 1,688.3852 4,734.5232 0


xmap_10CV_pred(m1,m3, 10, 55, anim_plot_data1) #27, 8: 39,066.2894 457,079.7451 1,430,608,083.1768 0   
xmap_10CV_pred2(m1,m3, 10, 55, anim_plot_data1) #21, 7: 39,224.2646 450,245.4961 0  

xmap_10CV_pred(m2,m4, 8, 66, anim_plot_data3) #18, 7: 12,004.8512 66,120.0122 21,796,526.4702 0   
xmap_10CV_pred2(m2,m4, 8, 66, anim_plot_data3) #16, 7: 11,950.8747 66,468.6363 0

xmap_10CV_pred(m2,m4, 15, 66, anim_plot_data3) #23, 7: 12,365.2133 92,621.5645 163,000,834.3901 0     
xmap_10CV_pred2(m2,m4, 15, 66, anim_plot_data3) #19, 7: 12,013.4649 92,671.5933 0 


#get this warning a lot: system is singular; attempting approx solution

w1_xcvn <- all_cv_x_pred2(m1,m3,anim_plot_data1) #seems to be having more issues with linearly depedent columns (og= 2154, new-2143, diff-9)

w2_xcvn <- all_cv_x_pred2(m5,m6,anim_plot_data2) #og - 3873, new - 3859

w3_xcvn <- all_cv_x_pred2(m2,m4,anim_plot_data3) #og - 5008, new - 5005

w1_ycvn <- all_cv_y_pred2(m1,m3,anim_plot_data1) #og - 2128, new - 2150

w2_ycvn <- all_cv_y_pred2(m5,m6,anim_plot_data2) #both with 3859

w3_ycvn <- all_cv_y_pred2(m2,m4,anim_plot_data3) #og - 5008, new - 5005

e_w1_xn <- cv_error_x_p(m1, anim_plot_data1, w1_xcvn) #10.575
e_w2_xn <- cv_error_x_p(m5, anim_plot_data2, w2_xcvn) #9.882
e_w3_xn <- cv_error_x_p(m2, anim_plot_data3, w3_xcvn) #10.078

e_w1_yn <- cv_error_y_p(m1, anim_plot_data1, w1_ycvn) #10.276
e_w2_yn <- cv_error_y_p(m5, anim_plot_data2, w2_ycvn) #9.377408
e_w3_yn <- cv_error_y_p(m2, anim_plot_data3, w3_ycvn) #9.559


e_w1_df_xn <- cv_df_p(m1, anim_plot_data1,w1_xcvn) 
e_w2_df_xn <- cv_df_p(m5, anim_plot_data2,w2_xcvn) 
e_w3_df_xn <- cv_df_p(m2, anim_plot_data3,w3_xcvn)

e_w1_df_yn <- cv_df_p(m1, anim_plot_data1,w1_ycvn) 
e_w2_df_yn <- cv_df_p(m5, anim_plot_data2,w2_ycvn)
e_w3_df_yn <- cv_df_p(m2, anim_plot_data3,w3_ycvn)

c_w1_xn <- e_w1_df_xn %>% group_by(clust) %>% summarise(x = sqrt((sum((pred-xmap)^2)/n()))) # 1-3.07, 2-4.36, 3-3.08, 4-3.1, 5 - 3.91, 6 - 3.1
c_w1_yn <- e_w1_df_yn %>% group_by(clust) %>% summarise(y = sqrt((sum((pred-ymap)^2)/n()))) # 1-4.1, 2-4.3, 3-3.1, 4-2.85, 5 - 3.67, 6 - 3.06
c_w2_xn <- e_w2_df_xn %>% group_by(clust) %>% summarise(x = sqrt((sum((pred-xmap)^2)/n()))) # 1-3, 2-2.94, 3-2.94, 4-3.26, 5 - 4.23
c_w2_yn <- e_w2_df_yn %>% group_by(clust) %>% summarise(y = sqrt((sum((pred-ymap)^2)/n()))) # 1-2.93, 2-2.89, 3-2.90, 4-3.13, 5 - 4.16
c_w3_xn <- e_w3_df_xn %>% group_by(clust) %>% summarise(x = sqrt((sum((pred-xmap)^2)/n()))) # 1-2.85, 2-2.81, 3-3.65, 4-2.67, 5 - 3.16, 6 - 3.21
c_w3_yn <- e_w3_df_yn %>% group_by(clust) %>% summarise(y = sqrt((sum((pred-ymap)^2)/n()))) # 1-2.85, 2-2.95, 3-3.59, 4-2.78, 5 - 3.04, 6 - 3.01


# Coverage --------------

sd_mean_w1xn <- mean(e_w1_df_xn$sd) #0.788
sd_mean_w1yn <- mean(e_w1_df_yn$sd) #0.834

e_w1_df_xn$low2 <- e_w1_df_xn$pred - (2*e_w1_df_xn$sd)
e_w1_df_xn$up2 <- e_w1_df_xn$pred + (2*e_w1_df_xn$sd)

e_w1_df_yn$low2 <- e_w1_df_yn$pred - (2*e_w1_df_yn$sd)
e_w1_df_yn$up2 <- e_w1_df_yn$pred + (2*e_w1_df_yn$sd)

e_w1_df_xn$inint2 <- 0
e_w1_df_xn$inint2[e_w1_df_xn$low2 < e_w1_df_xn$xmap & e_w1_df_xn$xmap < e_w1_df_xn$up2] <- 1

e_w1_df_yn$inint2 <- 0
e_w1_df_yn$inint2[e_w1_df_yn$low2 < e_w1_df_yn$ymap & e_w1_df_yn$ymap < e_w1_df_yn$up2] <- 1

w1_cp_x2n <- sum(e_w1_df_xn$inint2)/nrow(e_w1_df_xn) #coverage probability - 0.256
w1_cp_y2n <- sum(e_w1_df_yn$inint2)/nrow(e_w1_df_yn) #coverage probability - 0.293

##

sd_mean_w2xn <- mean(e_w2_df_xn$sd) #0.7509988
sd_mean_w2yn <- mean(e_w2_df_yn$sd) #0.8104493

e_w2_df_xn$low2 <- e_w2_df_xn$pred - (2*e_w2_df_xn$sd)
e_w2_df_xn$up2 <- e_w2_df_xn$pred + (2*e_w2_df_xn$sd)

e_w2_df_yn$low2 <- e_w2_df_yn$pred - (2*e_w2_df_yn$sd)
e_w2_df_yn$up2 <- e_w2_df_yn$pred + (2*e_w2_df_yn$sd)

e_w2_df_xn$inint2 <- 0
e_w2_df_xn$inint2[e_w2_df_xn$low2 < e_w2_df_xn$xmap & e_w2_df_xn$xmap < e_w2_df_xn$up2] <- 1

e_w2_df_yn$inint2 <- 0
e_w2_df_yn$inint2[e_w2_df_yn$low2 < e_w2_df_yn$ymap & e_w2_df_yn$ymap < e_w2_df_yn$up2] <- 1

w2_cp_x2n <- sum(e_w2_df_xn$inint2)/nrow(e_w2_df_xn) #coverage probability - 0.2630683
w2_cp_y2n <- sum(e_w2_df_yn$inint2)/nrow(e_w2_df_yn) #coverage probability - 0.2980793

##

sd_mean_w3xn <- mean(e_w3_df_xn$sd) #0.7600099
sd_mean_w3yn <- mean(e_w3_df_yn$sd) #0.7960184

e_w3_df_xn$low2 <- e_w3_df_xn$pred - (2*e_w3_df_xn$sd)
e_w3_df_xn$up2 <- e_w3_df_xn$pred + (2*e_w3_df_xn$sd)

e_w3_df_yn$low2 <- e_w3_df_yn$pred - (2*e_w3_df_yn$sd)
e_w3_df_yn$up2 <- e_w3_df_yn$pred + (2*e_w3_df_yn$sd)

e_w3_df_xn$inint2 <- 0
e_w3_df_xn$inint2[e_w3_df_xn$low2 < e_w3_df_xn$xmap & e_w3_df_xn$xmap < e_w3_df_xn$up2] <- 1

e_w3_df_yn$inint2 <- 0
e_w3_df_yn$inint2[e_w3_df_yn$low2 < e_w3_df_yn$ymap & e_w3_df_yn$ymap < e_w3_df_yn$up2] <- 1

w3_cp_x2n <- sum(e_w3_df_xn$inint2)/nrow(e_w3_df_xn) #0.246385
w3_cp_y2n <- sum(e_w3_df_yn$inint2)/nrow(e_w3_df_yn) #0.2630986



##########################################
## No Intersection - Spatial Only 
##########################################

# Simulation Data --------------------------------------------------------------

lat_sim_noint <- function(tp, frame, frame2){
  p = polyg_noint(frame2, tp)
  grid <- grid_tp_sim_noint(frame2,tp,p)
  known1 <- sim_known_noint(tp-1,grid,frame)
  known2 <- sim_known_noint(tp,grid, frame)
  known3 <- sim_known_noint(tp+1,grid, frame)
  known4 <- sim_known_noint(tp+2,grid, frame)
  if(nrow(known3) ==0){known <- rbind(known1, known2, known4)}
  else(known <- rbind(known1, known2, known3))
  loc <- known[,c("xmap","ymap", "obs_time")] ## locations that I will use to fit the model with time (xmap, ymap,obs_time)
  locs<-as.matrix(loc)                                ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)                ## intercept
  N <- nrow(known)
  Y <- known[,4]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #needs to go to number of rows - 1
  known_grd <- sim_estimate_noint(tp,grid, frame)
  known_grd2 <- filter(known_grd, !is.na(gpid_id))
  locs_pred <- as.matrix(known_grd2[,c("xmap", "ymap","obs_time")]) # Would need to create obs_time when don't have them (round number?)
  colnames(locs_pred) <- c("xmap", "ymap", "obs_time")
  X_pred <- matrix(c(rep(1,nrow(locs_pred)),locs_pred[,3]), ncol = 2, byrow=FALSE)
  pred_grid<-predictions(fit,locs_pred[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, 
                         y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  return(pred_grid) 
}

long_sim_noint <- function(tp, frame, frame2){
  p = polyg_noint(frame2, tp)
  grid <- grid_tp_sim_noint(frame2,tp,p)
  known1 <- sim_known_noint(tp-1,grid,frame)
  known2 <- sim_known_noint(tp,grid, frame)
  known3 <- sim_known_noint(tp+1,grid, frame)
  known4 <- sim_known_noint(tp+2,grid, frame)
  if(nrow(known3) ==0){known <- rbind(known1, known2, known4)}
  else(known <- rbind(known1, known2, known3))
  loc <- known[,c("xmap","ymap", "obs_time")] ## locations that I will use to fit the model with time (xmap, ymap,obs_time)
  locs<-as.matrix(loc)                                ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)                ## intercept
  N <- nrow(known)
  Y <- known[,5]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #needs to go to number of rows - 1
  known_grd <- sim_estimate_noint(tp,grid, frame)
  known_grd2 <- filter(known_grd, !is.na(gpid_id))
  locs_pred <- as.matrix(known_grd2[,c("xmap", "ymap","obs_time")]) # Would need to create obs_time when don't have them (round number?)
  colnames(locs_pred) <- c("xmap", "ymap", "obs_time")
  X_pred <- matrix(c(rep(1,nrow(locs_pred)),locs_pred[,3]), ncol = 2, byrow=FALSE)
  pred_grid<-predictions(fit,locs_pred[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, 
                         y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  return(pred_grid) 
}

#Loop Through Each Time Point

pred_all <- NULL
all_sim_lat_noint <- function(df, df2){
  for (r in min(df$t):max(df$t)){
    tryCatch({
      pred = lat_sim_noint(r,df, df2)
      pred_df = data.frame(pred,r)
      p = polyg_noint(df2, r)
      grid_t <- grid_tp_sim_noint(df2,r,p)
      known2 <- sim_estimate_noint(r,grid_t, df)
      pred_df = data.frame(x = pred, t=r, gpid = known2$gpid)
      pred_all = rbind(pred_all, pred_df)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(pred_all)
}

all_sim_long_noint <- function(df, df2){
  for (r in min(df$t):max(df$t)){
    tryCatch({
      pred = long_sim_noint(r,df, df2)
      pred_df = data.frame(pred,r)
      p = polyg_noint(df2, r)
      grid_t <- grid_tp_sim_noint(df2,r,p)
      known2 <- sim_estimate_noint(r,grid_t, df)
      pred_df = data.frame(y = pred, t=r, gpid = known2$gpid)
      pred_all = rbind(pred_all, pred_df)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(pred_all)
}


## Simulation 1

x_sim1_noint <- all_sim_lat_noint(w1, f_sim1)
y_sim1_noint <- all_sim_long_noint(w1,f_sim1)


xy_noint1 <- data.frame(x_sim1_noint$gpid, x_sim1_noint$t, x_sim1_noint$x, y_sim1_noint$y)
colnames(xy_noint1) <- c("gpid", "t", "xmap_pred", "ymap_pred")
xy_join_noint1 <- right_join(f_sim1, xy_noint1, by = c("gpid", "t"))

noint_x_sim1 <- sqrt(sum((xy_join_noint1$xmap-xy_join_noint1$xmap_pred)^2)/nrow(xy_join_noint1)) #1.435
noint_y_sim1 <- sqrt(sum((xy_join_noint1$ymap-xy_join_noint1$ymap_pred)^2)/nrow(xy_join_noint1)) #1.295

x1_c1_noint <- filter(xy_join_noint1, clust == 1)
x1_c2_noint <- filter(xy_join_noint1, clust == 2)

sim1_x_c1_noi <- sqrt(sum((x1_c1_noint$xmap-x1_c1_noint$xmap_pred)^2)/nrow(x1_c1_noint)) #1.536
sim1_y_c1_noi <- sqrt(sum((x1_c1_noint$ymap-x1_c1_noint$ymap_pred)^2)/nrow(x1_c1_noint)) #1.2904
sim1_x_c2_noi <- sqrt(sum((x1_c2_noint$xmap-x1_c2_noint$xmap_pred)^2)/nrow(x1_c2_noint)) #1.354
sim1_y_c2_noi <- sqrt(sum((x1_c2_noint$ymap-x1_c2_noint$ymap_pred)^2)/nrow(x1_c2_noint)) #1.298


#Simulation 2

x_sim2_noint <- all_sim_lat_noint(w2, f_sim2)
y_sim2_noint <- all_sim_long_noint(w2,f_sim2)

xy_noint2 <- data.frame(x_sim2_noint$gpid, x_sim2_noint$t, x_sim2_noint$x, y_sim2_noint$y)
colnames(xy_noint2) <- c("gpid", "t", "xmap_pred", "ymap_pred")
xy_join_noint2 <- right_join(f_sim2, xy_noint2, by = c("gpid", "t"))

noint_x_sim2 <- sqrt(sum((xy_join_noint2$xmap-xy_join_noint2$xmap_pred)^2)/nrow(xy_join_noint2)) #11.263
noint_y_sim2 <- sqrt(sum((xy_join_noint2$ymap-xy_join_noint2$ymap_pred)^2)/nrow(xy_join_noint2)) #1.488


x2_c1_noint <- filter(xy_join_noint2, clust == 1)
x2_c2_noint <- filter(xy_join_noint2, clust == 2)


sim2_x_c1_noi <- sqrt(sum((x2_c1_noint$xmap-x2_c1_noint$xmap_pred)^2)/nrow(x2_c1_noint)) #16.361
sim2_y_c1_noi <- sqrt(sum((x2_c1_noint$ymap-x2_c1_noint$ymap_pred)^2)/nrow(x2_c1_noint)) #1.596
sim2_x_c2_noi <- sqrt(sum((x2_c2_noint$xmap-x2_c2_noint$xmap_pred)^2)/nrow(x2_c2_noint)) #8.943
sim2_y_c2_noi <- sqrt(sum((x2_c2_noint$ymap-x2_c2_noint$ymap_pred)^2)/nrow(x2_c2_noint)) #1.451

#Simulation 3

x_sim3_noint <- all_sim_lat_noint(w3, f_sim3)
y_sim3_noint <- all_sim_long_noint(w3,f_sim3)

xy_noint3 <- data.frame(x_sim3_noint$gpid, x_sim3_noint$t, x_sim3_noint$x, y_sim3_noint$y)
colnames(xy_noint3) <- c("gpid", "t", "xmap_pred", "ymap_pred")
xy_join_noint3 <- right_join(f_sim3, xy_noint3, by = c("gpid", "t"))

noint_x_sim3 <- sqrt(sum((xy_join_noint3$xmap-xy_join_noint3$xmap_pred)^2)/nrow(xy_join_noint3)) #8.97
noint_y_sim3 <- sqrt(sum((xy_join_noint3$ymap-xy_join_noint3$ymap_pred)^2)/nrow(xy_join_noint3)) #1.489


x3_c1_noint <- filter(xy_join_noint3, clust == 1)
x3_c2_noint <- filter(xy_join_noint3, clust == 2)

sim3_x_c1_noi <- sqrt(sum((x3_c1_noint$xmap-x3_c1_noint$xmap_pred)^2)/nrow(x3_c1_noint)) #5.062
sim3_y_c1_noi <- sqrt(sum((x3_c1_noint$ymap-x3_c1_noint$ymap_pred)^2)/nrow(x3_c1_noint)) #1.405
sim3_x_c2_noi <- sqrt(sum((x3_c2_noint$xmap-x3_c2_noint$xmap_pred)^2)/nrow(x3_c2_noint)) #11.989
sim3_y_c2_noi <- sqrt(sum((x3_c2_noint$ymap-x3_c2_noint$ymap_pred)^2)/nrow(x3_c2_noint)) #1.581


## No Intersection----------------------------------------------------------------

one <- Sys.time()
b1 <- ymap_10CV_all(53, anim_plot_data1)
two <- Sys.time()

xmap_10CV_all <- function(tp, frame){
  p <- polyg_noint(frame, tp)
  grid <- grid_tp(frame, tp , p)
  obs1 <- known_all(frame, tp-1, grid)
  obs2 <- known_all(frame, tp, grid)
  obs3 <- known_all(frame, tp+1, grid)
  obs4 <- known_all(frame, tp+2, grid)
  set.seed(4)
  rand_df <- obs2[sample(nrow(obs2), size=(round(0.9*nrow(obs2)))),]
  if(nrow(obs3) ==0){obs <- rbind(obs1, rand_df, obs4)}else{obs <- rbind(obs1, rand_df, obs3)}
  loc <- obs[,c("xmap","ymap", "obs_time")]
  locs<-as.matrix(loc)
  N <- nrow(locs)
  ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)
  Y<- obs[,4]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = FALSE) #Silent option doesn't show the iterations. Maybe nice so get list of errors
  predicted_locs <- known_all_grid(frame, tp, grid)
  est_df <- subset(predicted_locs,!predicted_locs$gpid%in%rand_df$gpid)
  pred_loc <- as.matrix(est_df[ ,c("xmap","ymap", "obs_time")])
  X_pred <- matrix(c(rep(1,nrow(pred_loc)),pred_loc[,3]), ncol = 2, byrow=FALSE)
  ## for prediction
  pred_grid<-predictions(fit,pred_loc[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  preds <- data.frame(est_df$gpid,pred_grid)
  return(preds)
}

ymap_10CV_all <- function(tp, frame){
  p <- polyg_noint(frame, tp)
  grid <- grid_tp(frame, tp , p)
  obs1 <- known_all(frame, tp-1, grid)
  obs2 <- known_all(frame, tp, grid)
  obs3 <- known_all(frame, tp+1, grid)
  obs4 <- known_all(frame, tp+2, grid)
  set.seed(4)
  rand_df <- obs2[sample(nrow(obs2), size=(round(0.9*nrow(obs2)))),]
  if(nrow(obs3) ==0){obs <- rbind(obs1, rand_df, obs4)}else{obs <- rbind(obs1, rand_df, obs3)}
  loc <- obs[,c("xmap","ymap", "obs_time")]
  locs<-as.matrix(loc)
  N <- nrow(locs)
  ## convert them into matrix
  X <- matrix(c(rep(1,nrow(locs)),locs[,3]), ncol = 2, byrow=FALSE)
  Y<- obs[,5]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #Silent option doesn't show the iterations. Maybe nice so get list of errors
  predicted_locs <- known_all_grid(frame, tp, grid)
  est_df <- subset(predicted_locs,!predicted_locs$gpid%in%rand_df$gpid)
  pred_loc <- as.matrix(est_df[ ,c("xmap","ymap", "obs_time")])
  X_pred <- matrix(c(rep(1,nrow(pred_loc)),pred_loc[,3]), ncol = 2, byrow=FALSE)
  ## for prediction
  pred_grid<-predictions(fit,pred_loc[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  preds <- data.frame(est_df$gpid,pred_grid)
  return(preds)
}


all_pred_noi <- data.frame(NULL)
pred_cv_noi <- data.frame(NULL)


#xmap
all_cv_x_noint <- function(df){
  for (g in min(df$t):max(df$t)){
    #for (g in 51:52){
    tryCatch({
      pred = xmap_10CV_all(g, df)
      pred_df = data.frame(pred,g)
      pred_cv_noi = rbind.data.frame(pred_cv_noi, pred_df)}, error=function(e){return(NA)})
  }
  return(pred_cv_noi)
}

#ymap
all_cv_y_noint <- function(df){
  for (g in min(df$t):max(df$t)){
    tryCatch({
      pred = ymap_10CV_all(g, df)
      pred_df = data.frame(pred,g)
      pred_cv_noi = rbind.data.frame(pred_cv_noi, pred_df)}, error=function(e){return(NA)})
  }
  return(pred_cv_noi)
}

## Get Error
cv_error_x_noi <- function(frame, pred_df){
  final <- filter(frame, imputed =="FALSE")[,c(1,2,5,6,9)]
  colnames(pred_df) <- c("gpid", "pred", "t")
  pred_df$pred <- round(pred_df$pred, 3)
  all_data <- left_join(pred_df,final, by = c("gpid", "t"))
  all_data_distinct <- distinct(all_data)
  error <- sum((all_data_distinct$pred-all_data_distinct$xmap)^2)/nrow(all_data_distinct)
  return(error)
}

cv_error_y_noi <- function(frame, pred_df){
  final <- filter(frame, imputed =="FALSE")[,c(1,2,5,6,9)]
  colnames(pred_df) <- c("gpid", "pred", "t")
  pred_df$pred <- round(pred_df$pred, 3)
  all_data <- left_join(pred_df,final, by = c("gpid", "t"))
  all_data_distinct <- distinct(all_data)
  error <- sum((all_data_distinct$pred-all_data_distinct$ymap)^2)/nrow(all_data_distinct)
  return(error)
}

#Gives Dataframe with actual and predicted
cv_df_noi <- function(frame, pred_df){
  final <- filter(frame, imputed =="FALSE")[,c(1,2,5,6,9)]
  colnames(pred_df) <- c("gpid", "pred", "t")
  pred_df$pred <- round(pred_df$pred, 3)
  all_data <- left_join(pred_df,final, by = c("gpid", "t"))
  all_data_distinct <- distinct(all_data)
  return(all_data_distinct)
}

h1_x <- data.frame(xmap_10CV_all(57, anim_plot_data1), r = 57)

h2_x <- all_cv_x_noint(anim_plot_data2)
h3_x <- all_cv_x_noint(anim_plot_data3)

h1_y <- all_cv_y_noint(anim_plot_data1)
h2_y <- all_cv_y_noint(anim_plot_data2)
h3_y <- all_cv_y_noint(anim_plot_data3)

w1x <- cv_error_x_noi(anim_plot_data1,h1_x) #9.573 - 3.093
w2x <- cv_error_x_noi(anim_plot_data2,h2_x) #10.2287 - 3.198
w3x <- cv_error_x_noi(anim_plot_data3,h3_x) #10.175 - 3.1903

w1y <- cv_error_y_noi(anim_plot_data1,h1y) #8.287- 3.067 (55, 57)
w2y <- cv_error_y_noi(anim_plot_data2,h2y) #9.923 - 3.151
w3y <- cv_error_y_noi(anim_plot_data3,h3y) # 9.926 - 3.150238

