## Attempt 1 ----------------------------------------------------------------------------------------

#I feel like this method is better at the moment (mathces what doing with gpgp package in paper at the moment)

#Following this example: https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html#sec:geostats
#Spatial Model with Random/Fixed Time Effect


val2 <- data_inla(f_sim2, 4, polyg2, w2, n2, 2) #Get Data set up properly (t=3, int = 1)
val2_known <- filter(val2, !is.na(xmap))

plot(val2_known$xmap, val2_known$ymap)


#Create Mesh - nonconvex (create around known values)
prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(130, 100))

prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(0.1, 5), cutoff = 0.75,
                        offset = c(-0.1, -0.1))

#par(mfrow=c(1,2))
plot(prmesh1, asp = 1, main = "") #see what mesh looks like
#plot(val2_known$xmap, val2_known$ymap)

meuse.spde <- inla.spde2.matern(mesh = prmesh1, alpha = 2)
A.meuse <- inla.spde.make.A(mesh = prmesh1, loc = coordinates(val2_known[,13:14]))
s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = meuse.spde$n.spde)


meuse.stack <- inla.stack(data  = list(x = val2_known$xmap),
                          A = list(A.meuse, 1),
                          effects = list(c(s.index, list(Intercept = 1)),
                                         list(t = val2_known$t)),
                          tag = "meuse.data")

val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = coordinates(val2_pred2[,13:14]))
meuse.stack.pred <- inla.stack(data = list(x = NA),
                               A = list(A.pred, 1),
                               effects = list(c(s.index, list(Intercept = 1)),
                                              list(t = val2_pred2$t)),
                               tag = "meuse.pred")


join.stack <- inla.stack(meuse.stack, meuse.stack.pred)

#Fixed/random time effect - essentially gives same answer
form <- x ~ -1 + Intercept + f(spatial.field, model = spde) + t
  #f(t, model = "rw1")
  

m_new <- inla(form, data = inla.stack.data(join.stack, spde = meuse.spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE))

summary(m_new) #fixed effect mean is 0 (not adding anything to model, not big enough?)

index.pred <- inla.stack.index(join.stack, "meuse.pred")$data

xest <- m_new$summary.fitted.values[index.pred, "mean"]
xsd <- m_new$summary.fitted.values[index.pred, "sd"] #not as small as it was for gpgp


val2_pred2$x_grid ## Little bit different in decimals (so only moving a little bit)
val2_pred2$x_actual



##Attempt 2 -------------------------------------------------------------------------------

#Separable Space-Time Model (Discrete Time) - https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html#discrete-time-domain

val2 <- data_inla(f_sim2, 3, polyg2, w2, n2, 2) #Get Data set up properly (t=3, int = 1)
val2_known <- filter(val2, !is.na(xmap))


#Create Mesh - nonconvex (create around known values)
prdomain <- inla.nonconvex.hull(as.matrix(val2_known[, 5:6]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(120, 100))

prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(4, 3), cutoff = 0.8,
                        offset = c(-0.05, -0.05))

plot(prmesh1, asp = 1, main = "") #see what mesh looks like


spde2 <- inla.spde2.pcmatern(mesh = prmesh1, 
                            prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01


#k <- max(val2_known$t)

k <- 3

iset <- inla.spde.make.index('i', n.spde = spde2$n.spde, n.group = k)


A <- inla.spde.make.A(mesh = prmesh1,
                      loc = cbind(val2_known$xval, val2_known$yval), group = val2_known$t-1) 


sdat <- inla.stack(
  data = list(x = val2_known$xmap), 
  A = list(A), 
  effects = list(c(iset, list(Intercept = 1))),
  tag = 'stdata') 



val2_pred <- filter(val2, is.na(xmap))
val2_pred2 <- val2_pred %>% distinct(x_grid,y_grid,.keep_all= TRUE)

A.pred <- inla.spde.make.A(mesh = prmesh1, loc = cbind(val2_pred2$xval, val2_pred2$yval), group = val2_pred2$t)
#iset2 <- inla.spde.make.index('i', n.spde = spde2$n.spde, n.group = max(val2_pred$t))

sdat.pred <- inla.stack(
  data = list(x = NA), 
  A = list(A.pred), 
  effects = list(c(iset, list(Intercept = 1))),
  tag = 'stdata.pred') 

#Join stack
join.stack <- inla.stack(sdat, sdat.pred)


#Fit model

h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
# Model formula
formulae <- x ~ -1 + Intercept + f(i, model = spde2, group = i.group, 
                      control.group = list(model = 'ar1')) #, hyper = h.spec)) 


# PC prior on the autoreg. param.
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
# Model fitting
res3 <- inla(formulae,  data = inla.stack.data(join.stack), 
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(join.stack)), 
            control.family = list(hyper = list(prec = prec.prior)), 
            control.fixed = list(expand.factor.strategy = 'inla'))

#Summary of results
summary(res3)

ndex.pred <- join.stack$data$index$stdata.pred


spde.mn2 <- res3$summary.fitted.values[ndex.pred, "mean"]
spde.sd2<- res3$summary.fitted.values[ndex.pred, "sd"]

#Gives a very small number for 4th missing data point. Attempt 1 does not do that. 

#Currently using INLA version 22.05.07 (update?)

#After updated Attempt 2 not longer works, not getting good enough initial values (issue with model specification?)

inla.upgrade() #version 22.12.16

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)





###########################################################################


val <- cbind(val2_pred2, pred = xest, sd = xsd)

x_error <- sqrt(sum((val$x_actual-val$pred)^2)/nrow(val)) #1.161 - little bit better than grid

val$ll <- val$pred - (2*val$sd)
val$ul <- val$pred + (2*val$sd)

val$interval <- 0
val$interval[val$ll < val$x_actual & val$x_actual < val$ul] <- 1

sum(val$interval)/nrow(val)
