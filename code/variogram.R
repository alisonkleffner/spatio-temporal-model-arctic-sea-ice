## Create a Variogram

b <- variogram(xmap ~ t , ~xmap+ymap, data= anim_plot_data3)
plot(b)

c <- variogram(xmap ~ t , ~xmap+ymap, data= anim_plot_data2)
plot(c)

a <- variogram(xmap ~ t , ~xmap+ymap, data= anim_plot_data1, width = 10)
plot(a)

a2 <- variogram(xmap ~ t , ~xmap+ymap, data= anim_plot_data1, cutoff = 900)
plot(a2)


#simulated data
d <- variogram(xmap ~ t , ~xmap+ymap, data= f1, cutoff = 50, width = 1)
plot(d)

#Never seems to level out - strong spatial correlation



vv <- variogram(xmap ~ 1, ~ xmap +ymap + t,
                data = anim_plot_data1, 
                width = 10, 
                cutoff = 900,
                alpha = 45)

plot(vv)


rm(data)
rm(XY)
rm(idx.t)
rm(idx.x)
rm(idx.y)
rm(t.joint)
rm(x.joint)
rm(y.joint)

