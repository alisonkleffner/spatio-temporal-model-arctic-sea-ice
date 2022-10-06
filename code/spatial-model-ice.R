##########Spatial Model with Time in mean function - Actual Data ##############


#Week1_x model error <- 3.156
#Week1_y model error <- 3.239
#Week1_x linear error <- 1.527
#Week1_y linear error <- 4.179

#Week2_x model error <- 3.142
#Week2_y model error <- 3.062
#Week2_x linear error <- 2
#Week2_y linear error <- 1.399

#Week3_x model error <- 3.178
#Week3_y model error <- 3.094
#Week3_x linear error <- 1.027
#Week3_y linear error <- 1.215


####################
#Updated Functions
####################

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
  X <- as.matrix(locs[,3])                  ## intercept
  N <- nrow(obs)
  Y <- obs[,5]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = FALSE) #needs to go to number of rows - 1  
  predicted_locs <- known_data_grid(m_new,int,tp,grid, frame)
  est_df <- subset(predicted_locs,!(predicted_locs$gpid%in%rand_df$gpid))
  pred_loc <- as.matrix(est_df[ ,c("xmap","ymap", "obs_time")])
  X_pred <- as.matrix(pred_loc[,3])
  ncondsim <- 30
  sims <- cond_sim(fit = fit, locs_pred = pred_loc[,1:2], X_pred = X_pred, covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                   ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE, nsims = ncondsim)
  ##for prediction
  pred_grid<-predictions(fit,pred_loc[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, 
                         y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  preds <- data.frame(est_df$gpid,pred_grid)
  b <- list(sims,preds)
  return(b)
}


#Taking a longer time to converge

startTime <- Sys.time()
v <- xmap_10CV_pred2(m2, m4, 8, 65, anim_plot_data3) #(323, 462)
endTime <- Sys.time()
#Time difference of 2.136734 mins


startTime2 <- Sys.time()
v2 <- xmap_10CV_pred(m2, m4, 8, 65, anim_plot_data3) #(original) (20,7)
endTime2 <- Sys.time()
#Time difference of 6.564464 secs

#72,208.863 405,585.4581 0 (spatial)
#12,097.9566 67,730.5095 42,698,712.9681 0 (spatio-temporal)


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
  X <- as.matrix(locs[,3])                  ## intercept
  N <- nrow(obs)
  Y <- obs[,6]
  fit <- fit_model(Y,locs[,1:2], X, "exponential_isotropic", max_iter = 1000, convtol = 1e-04,reorder = TRUE, m_seq = c(10,min(30,N-1)), silent = TRUE) #needs to go to number of rows - 1  
  predicted_locs <- known_data_grid(m_new,int,tp,grid, frame)
  est_df <- subset(predicted_locs,!(predicted_locs$gpid%in%rand_df$gpid))
  pred_loc <- as.matrix(est_df[ ,c("xmap","ymap", "obs_time")])
  X_pred <- as.matrix(pred_loc[,3])
  ncondsim <- 30
  sims <- cond_sim(fit = fit, locs_pred = pred_loc[,1:2], X_pred = X_pred, covparms = fit$covparms, covfun_name = fit$covfun_name, y_obs = fit$y,locs_obs = fit$locs
                   ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE, nsims = ncondsim)
  ##for prediction
  pred_grid<-predictions(fit,pred_loc[,1:2], X_pred,covparms = fit$covparms, covfun_name = fit$covfun_name, 
                         y_obs = fit$y,locs_obs = fit$locs
                         ,X_obs= fit$X, beta=fit$betahat, m = 10, reorder = TRUE)
  preds <- data.frame(est_df$gpid,pred_grid)
  b <- list(sims,preds)
  return(b)
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

###################
## Week 1
###################

startTime <- Sys.time()
w1_xcv <- all_cv_x_pred2(m1,m3,anim_plot_data1)
endTime <- Sys.time()

startTime - endTime

startTime <- Sys.time()
w1_ycv <- all_cv_y_pred(m1,m3,anim_plot_data1)
endTime <- Sys.time()

e_w1_x <- cv_error_x_p(m1, anim_plot_data1, w1_xcv) #
e_w1_y <- cv_error_y_p(m1, anim_plot_data1, w1_ycv) #


## Coverage Probability

sd_mean_w1x <- mean(e_w1_df_x$sd)
sd_mean_w1y <- mean(e_w1_df_y$sd)

e_w1_df_x <- cv_df_p(m1, anim_plot_data1,w1_xcv) 
e_w1_df_y <- cv_df_p(m1, anim_plot_data1,w1_ycv) 

e_w1_df_x$low2 <- e_w1_df_x$pred - (2*e_w1_df_x$sd)
e_w1_df_x$up2 <- e_w1_df_x$pred + (2*e_w1_df_x$sd)

e_w1_df_y$low2 <- e_w1_df_y$pred - (2*e_w1_df_y$sd)
e_w1_df_y$up2 <- e_w1_df_y$pred + (2*e_w1_df_y$sd)

e_w1_df_x$inint2 <- 0
e_w1_df_x$inint2[e_w1_df_x$low2 < e_w1_df_x$xmap & e_w1_df_x$xmap < e_w1_df_x$up2] <- 1

e_w1_df_y$inint2 <- 0
e_w1_df_y$inint2[e_w1_df_y$low2 < e_w1_df_y$ymap & e_w1_df_y$ymap < e_w1_df_y$up2] <- 1

sum(e_w1_df_x$inint2)/nrow(e_w1_df_x) #coverage probability - 
sum(e_w1_df_y$inint2)/nrow(e_w1_df_y) #coverage probability - 



###################
## Week 2
###################

startTime <- Sys.time()
w2_xcv <- all_cv_x_pred2(m5,m6,anim_plot_data2)
endTime <- Sys.time()

startTime - endTime

w2_ycv <- all_cv_y_pred2(m5,m6,anim_plot_data2)

e_w2_x <- cv_error_x_p(m5, anim_plot_data2, w2_xcv) #
e_w2_y <- cv_error_y_p(m5, anim_plot_data2, w2_ycv) #

## Coverage Probability

e_w2_df_x <- cv_df_p(m5, anim_plot_data2,w2_xcv)
e_w2_df_y <- cv_df_p(m5, anim_plot_data2,w2_ycv)

sd_mean_w2x <- mean(e_w2_df_x$sd)
sd_mean_w2y <- mean(e_w2_df_y$sd)

e_w2_df_x$low2 <- e_w2_df_x$pred - (2*e_w2_df_x$sd)
e_w2_df_x$up2 <- e_w2_df_x$pred + (2*e_w2_df_x$sd)

e_w2_df_y$low2 <- e_w2_df_y$pred - (2*e_w2_df_y$sd)
e_w2_df_y$up2 <- e_w2_df_y$pred + (2*e_w2_df_y$sd)

e_w2_df_x$inint2 <- 0
e_w2_df_x$inint2[e_w2_df_x$low2 < e_w2_df_x$xmap & e_w2_df_x$xmap < e_w2_df_x$up2] <- 1

e_w2_df_y$inint2 <- 0
e_w2_df_y$inint2[e_w2_df_y$low2 < e_w2_df_y$ymap & e_w2_df_y$ymap < e_w2_df_y$up2] <- 1

sum(e_w2_df_x$inint2)/nrow(e_w2_df_x) #coverage probability - 
sum(e_w2_df_y$inint2)/nrow(e_w2_df_y) #coverage probability - 



###################
## Week 3
###################

w3_xcv <- all_cv_x_pred2(m2,m4,anim_plot_data3)
w3_ycv <- all_cv_y_pred2(m2,m4,anim_plot_data3)


e_w3_x <- cv_error_x_p(m2, anim_plot_data3, w3_xcv) #
e_w3_y <- cv_error_y_p(m2, anim_plot_data3, w3_ycv) #

## Coverage Probability

e_w3_df_x <- cv_df_p(m2, anim_plot_data3,w3_xcv)
e_w3_df_y <- cv_df_p(m2, anim_plot_data3,w3_ycv)

sd_mean_w3x <- mean(e_w3_df_x$sd)
sd_mean_w3y <- mean(e_w3_df_y$sd)

e_w3_df_x$low2 <- e_w3_df_x$pred - (2*e_w3_df_x$sd)
e_w3_df_x$up2 <- e_w3_df_x$pred + (2*e_w3_df_x$sd)

e_w3_df_y$low2 <- e_w3_df_y$pred - (2*e_w3_df_y$sd)
e_w3_df_y$up2 <- e_w3_df_y$pred + (2*e_w3_df_y$sd)

e_w3_df_x$inint2 <- 0
e_w3_df_x$inint2[e_w3_df_x$low2 < e_w3_df_x$xmap & e_w3_df_x$xmap < e_w3_df_x$up2] <- 1

e_w3_df_y$inint2 <- 0
e_w3_df_y$inint2[e_w3_df_y$low2 < e_w3_df_y$ymap & e_w3_df_y$ymap < e_w3_df_y$up2] <- 1

sum(e_w3_df_x$inint2)/nrow(e_w3_df_x) #coverage probability - 
sum(e_w3_df_y$inint2)/nrow(e_w3_df_y) #coverage probability - 
