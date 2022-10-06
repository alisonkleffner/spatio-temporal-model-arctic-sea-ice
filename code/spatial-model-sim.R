##########Spatial Model with Time in mean function - Simulated Data ##############


#model error x1: 1.496
#model error y1: 1.517
#linear error x1: 1.042
#linear error y1: 1.226
#model error x2: 1.628
#model error y2: 1.579
#linear error x2: 1.455
#linear error y2: 1.539
#model error x3: 1.342
#model error y3: 1.338
#linear error x3: 0.95
#linear error y3: 0.92


### Function with spatial and time in mean

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
  X <- as.matrix(locs[,3])                  ## intercept
  N <- nrow(known)
  Y <- known[,5]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = FALSE) #needs to go to number of rows - 1
  known_grd <- sim_estimate(m_new,int,tp,grid, frame)
  locs_pred <- as.matrix(known_grd[,c("xmap", "ymap","obs_time")]) # Would need to create obs_time when don't have them (round number?)
  colnames(locs_pred) <- c("xmap", "ymap", "obs_time")
  X_pred <- as.matrix(locs_pred[,3])
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
  X <- as.matrix(locs[,3])                    ## intercept
  N <- nrow(known)
  Y <- known[,6]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #needs to go to number of rows - 1
  known_grd <- sim_estimate(m_new,int,tp,grid, frame)
  locs_pred <- as.matrix(known_grd[,c("xmap", "ymap","obs_time")]) # Would need to create obs_time when don't have them (round number?)
  colnames(locs_pred) <- c("xmap", "ymap", "obs_time")
  X_pred <- as.matrix(locs_pred[,3])
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
        df_p = data.frame(sd_df, i, r, known2$gpid)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(all_pred)
}

all_long2 <- function(df, df2, int_list, int_df){
  for (r in 1:7){
    for (i in 1:4){
      tryCatch({
        pred = long_sim2(i,int_list,int_df,r,df, df2)
        sd_df <- sd_pred_sim(pred[[1]],pred[[2]])
        grid_t <- grid_tp_sim(df2,r,int_list)
        known2 <- sim_estimate(int_df,i,r,grid_t, df)
        df_p = data.frame(sd_df, i, r, known2$gpid)
        pred_cv = rbind.data.frame(pred_cv, df_p)}, error=function(e){return(NA)})
    }
    all_pred <- rbind.data.frame(all_pred,pred_cv)
  }
  return(all_pred)
}

s <- lat_sim2(4, polyg3, n3, w3, f1)
s2 <- lat_sim(4, polyg3, n3, w3, f1)

##################
## Finding Errror for Simulation 1
##################


#Gives Location for X
int_1x = all_lat2(w1, f1,polyg1,n1)

#Gives Location for Y
int_1y = all_long2(w1,f1,polyg1,n1)

xy1 <- data.frame(int_1x$known2.gpid, int_1x$r, int_1x$i, int_1x$pred_df, int_1x$sd, int_1x$lc, int_1x$uc, int_1y$pred_df, int_1y$sd, int_1y$lc, int_1y$uc)
colnames(xy1) <- c("gpid", "t", "intersection", "xmap_pred", "x_sd", "x_lc", "x_uc", "ymap_pred", "y_sd", "y_lc", "y_uc")


xy_join1 <- right_join(f1, xy1, by = c("gpid", "t"))
xy_join1 <- xy_join1 %>% distinct(gpid,t,.keep_all= TRUE)


error_x_s1 <- sqrt(sum((xy_join1$xmap-xy_join1$xmap_pred)^2)/nrow(xy_join1)) #1.491
error_y_s1 <- sqrt(sum((xy_join1$ymap-xy_join1$ymap_pred)^2)/nrow(xy_join1)) #1.5

###################
## Finding Errror for Simulation 2
###################

#Gives Location for X
int_2x = all_lat2(w2, f1,polyg2,n2)


#Gives Location for Y
int_2y = all_long2(w2,f1,polyg2,n2)

xy2 <- data.frame(int_2x$known2.gpid, int_2x$r, int_2x$i, int_2x$pred_df, int_2x$sd, int_2x$lc, int_2x$uc, int_2y$pred_df, int_2y$sd, int_2y$lc, int_2y$uc)
colnames(xy2) <- c("gpid", "t", "intersection", "xmap_pred", "x_sd", "x_lc", "x_uc", "ymap_pred", "y_sd", "y_lc", "y_uc")


xy_join2 <- right_join(f1, xy2, by = c("gpid", "t"))
xy_join2 <- xy_join2 %>% distinct(gpid,t,.keep_all= TRUE)


error_x_s2 <- sqrt(sum((xy_join2$xmap-xy_join2$xmap_pred)^2)/nrow(xy_join2)) #1.620597
error_y_s2 <- sqrt(sum((xy_join2$ymap-xy_join2$ymap_pred)^2)/nrow(xy_join2)) #1.571889

###################
## Finding Errror for Simulation 3
###################

#Gives Location for X
int_3x = all_lat2(w3, f1,polyg3,n3)

#Gives Location for Y
int_3y = all_long2(w3,f1,polyg3,n3)

xy3 <- data.frame(int_3x$known2.gpid, int_3x$r, int_3x$i, int_3x$pred_df, int_3x$sd, int_3x$lc, int_3x$uc, int_3y$pred_df, int_3y$sd, int_3y$lc, int_3y$uc)
colnames(xy3) <- c("gpid", "t", "intersection", "xmap_pred", "x_sd", "x_lc", "x_uc", "ymap_pred", "y_sd", "y_lc", "y_uc")

xy_join3 <- right_join(f1, xy3, by = c("gpid", "t"))
xy_join3 <- xy_join3 %>% distinct(gpid,t, .keep_all= TRUE)

error_x_s3 <- sqrt(sum((xy_join3$xmap-xy_join3$xmap_pred)^2)/nrow(xy_join3)) #1.340666
error_y_s3 <- sqrt(sum((xy_join3$ymap-xy_join3$ymap_pred)^2)/nrow(xy_join3)) #1.331062


####################
## Check Coverage 
####################

## Simulation 1

xy_join1$inintx <- 0
xy_join1$inintx[xy_join1$x_lc < xy_join1$xmap & xy_join1$xmap < xy_join1$x_uc] <- 1
cp_sim1_x <- sum(xy_join1$inintx)/nrow(xy_join1) #0.298 (up)
sd_sim1_x <- mean(xy_join1$x_sd) #0.381 

xy_join1$ininty <- 0
xy_join1$ininty[xy_join1$y_lc < xy_join1$ymap & xy_join1$ymap < xy_join1$y_uc] <- 1
cp_sim1_y <- sum(xy_join1$ininty)/nrow(xy_join1) #0.281 (low)
sd_sim1_y <- mean(xy_join1$y_sd) #0.355

## Simulation 2

xy_join2$inintx <- 0
xy_join2$inintx[xy_join2$x_lc < xy_join2$xmap & xy_join2$xmap < xy_join2$x_uc] <- 1
cp_sim2_x <- sum(xy_join2$inintx)/nrow(xy_join2) #0.1875 (same)
sd_sim2_x <- mean(xy_join2$x_sd) #0.428

xy_join2$ininty <- 0
xy_join2$ininty[xy_join2$y_lc < xy_join2$ymap & xy_join2$ymap < xy_join2$y_uc] <- 1
cp_sim2_y <- sum(xy_join2$ininty)/nrow(xy_join2) #0.25 (same)
sd_sim2_y <- mean(xy_join2$y_sd) #0.439

## Simulation 3

xy_join3$inintx <- 0
xy_join3$inintx[xy_join3$x_lc < xy_join3$xmap & xy_join3$xmap < xy_join3$x_uc] <- 1
cp_sim3_x <- sum(xy_join3$inintx)/nrow(xy_join3) #0.3125 (up)
sd_sim3_x <- mean(xy_join3$x_sd) #0.376

xy_join3$ininty <- 0
xy_join3$ininty[xy_join3$y_lc < xy_join3$ymap & xy_join3$ymap < xy_join3$y_uc] <- 1
cp_sim3_y <- sum(xy_join3$ininty)/nrow(xy_join3) #0.396 (up)
sd_sim3_y <- mean(xy_join3$y_sd) #0.377

###################
##Notes
###################

##Generally only has slight improvements in both MSE and Coverage. MSE always better, but no by much. 
##Standard Deviations roughly the same.
