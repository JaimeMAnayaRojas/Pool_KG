# Make survival plot

surv_link <- function(post, size, G = 1, NK, NG, center){

  z = size - center

  if(G == 1){
    p =  inv_logit(with(post, Intercept_survG + b_NK_survG * NK + b_z_survG * z + b_zNK_survG* NK*z))
  }else{
    p =  inv_logit(with(post,Intercept_survK + b_NG_survK * NG + b_z_survK * z + b_zNG_survK* NG*z))
  }

  return(p)
}



grow_link <-function(post, size, G = 1, NK, NG, center){

  z = size - center

  if(G == 1){
    mu = with(post,Intercept_growG + b_NK_growG * NK + b_z_growG * z + b_zNK_growG * NK * z )
  }else{
    mu = with(post,Intercept_growK + b_NG_growK * NG + b_z_growK * z + b_zNG_growK * NG * z)
  }


  return(mu-size)
}




recr_link <-function(post, size, G = 1, NK, NG, center){

  z = size - center

  if(G == 1){
    lambda = with(post,Intercept_recrG + b_NK_recrG * NK + b_z_recrG * z + b_zNK_recrG * NK * z )
  }else{
    lambda = with(post, Intercept_recrK + b_NG_recrK * NG + b_z_recrK * z + b_zNG_recrK * NG * z )
  }
  return(exp(lambda))
}



jpeg(file = "Guppy_fig.jpeg",width = 7.5, height = 5, units='in', res=800)
op<-par(mfrow=c(2,3), mar = c(4, 4 , 2, 1), oma = c(0.5, 1, 1, 0.5)) # c(bottom, left, top, right)
 
x_axis_min = 5
x_axis_max = 32

z.c_grid <-  seq(from= x_axis_min, to=x_axis_max, length.out = 200)

pred.raw <- sapply(1:length(z.c_grid), function(i) surv_link(post = post, size = z.c_grid[i], G = 1, NK = 0, NG = 0, center = 18))
p.KG= pred.raw
p.KG.mean <- apply(pred.raw, 2, mean)
p.KG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

pred.raw <- sapply(1:length(z.c_grid), function(i) surv_link(post = post, size = z.c_grid[i], G = 1, NK = 1, NG = 0, center = 18))
p.NK= pred.raw
p.NK.mean <- apply(pred.raw, 2, mean)
p.NK.PI <- apply(pred.raw, 2, HPDI,prob = .95)


LOS(c(p.NK)-c(p.KG))

mean(c(p.NK)-c(p.KG))*100
HPDI((c(p.NK)-c(p.KG)), prob = .95)*100

toplot <- Gdata
toplot$SL1_mm
plot( surv ~ SL1_mm, data=toplot,  col=rangi2, ylab="Survival (%)",
      xlab= "", pch="",  xlim=c(5,30), xaxt="n",yaxt="n")
axis(1, labels = c(5,10,15,20,25,30), at = (c(5,10,15,20,25,30)) )
axis(2, labels = c(0, 20, 40, 60, 80, 100), at = (c(0, 0.2, 0.4, 0.6, 0.8, 1)) )


shade((p.KG.PI), z.c_grid, col=col.alpha("black", 0.2))
lines(z.c_grid, (p.KG.PI[1,]), col="black", lwd=1, lty=3)
lines(z.c_grid, (p.KG.PI[2,]), col="black", lwd=1, lty=3)

shade((p.NK.PI), z.c_grid, col=col.alpha("orange", 0.3))
lines(z.c_grid, (p.NK.PI[1,]), col="darkorange", lwd=1, lty=3)
lines(z.c_grid, (p.NK.PI[2,]), col="darkorange", lwd=1, lty=3)


points( surv ~ SL1_mm, bg=col.alpha("black", 0.2), pch=21, cex = 0.5,  data=toplot[toplot$KG == 1,])
points( surv ~ SL1_mm, bg="orange", pch=21, data=toplot[toplot$NK==1,])
lines(z.c_grid, (p.KG.mean), col="black", lwd=2)
lines(z.c_grid, (p.NK.mean), col="darkorange2", lwd=2)
title(main='a)', adj = 0, line = 0.5)
legend(x= 10, y= .3, c('Control (KG)','No killifish (NK)'),lty=c(1,1), pch=c(21,21), col=c('black','orange'), pt.bg = c('gray','orange'),bty='n')


gsize = seq(from= 7, to=14, length.out = 100)
surv10_NK_G <- sapply(1:length(gsize), function(i) surv_link(post = post, size = gsize[i], G = 1, NK = 1, NG = 0, center = 10))
surv10_KG_G <- sapply(1:length(gsize), function(i) surv_link(post = post, size = gsize[i], G = 1, NK = 0, NG = 0, center = 10))
mean(as.numeric(surv10_NK_G - surv10_KG_G))*100
HPDI(as.numeric(surv10_NK_G - surv10_KG_G))*100
100-LOS(as.numeric(surv10_NK_G - surv10_KG_G))


gsize = seq(from= 25, to=30, length.out = 100)
surv25_NK_G <-sapply(1:length(gsize), function(i) surv_link(post = post, size = gsize[i], G = 1, NK = 1, NG = 0, center = 10))
surv25_KG_G <-  sapply(1:length(gsize), function(i) surv_link(post = post, size = gsize[i], G = 1, NK = 0, NG = 0, center = 10))
mean(as.numeric(surv10_NK_G - surv25_NK_G))*100
HPDI(as.numeric(surv10_NK_G - surv25_NK_G))*100
100-LOS(as.numeric(surv10_NK_G - surv25_NK_G))

# Growth


pred.raw <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 1, NK = 0, NG = 0, center = 18))
p.KG= pred.raw
p.KG.mean <- apply(pred.raw, 2, mean)
p.KG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

pred.raw <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 1, NK = 1, NG = 0, center = 18))
p.NK= pred.raw
p.NK.mean <- apply(pred.raw, 2, mean)
p.NK.PI <- apply(pred.raw, 2, HPDI,prob = .95)

LOS(c(p.NK) - c(p.KG))
chainmode(c(p.NK) - c(p.KG))
HPDI(c(p.NK) - c(p.KG), prob = .95)


toplot <- subset(Gdata, surv == 1)
toplot$z1<- toplot$z1 - toplot$SL1_mm

plot( z1 ~ SL1_mm, data=toplot,  col=rangi2, ylab="Somatic growth (mm/28 days)",
      xlab= "", pch="",  xlim=c(5,30), xaxt="n")
axis(1, labels = c(5,10,15,20,25,30), at = (c(5,10,15,20,25,30)) )
#axis(2, labels = c(0, 20, 40, 60, 80, 100), at = (c(0, 0.2, 0.4, 0.6, 0.8, 1)) )

shade((p.KG.PI), z.c_grid, col=col.alpha("black", 0.2))
lines(z.c_grid, (p.KG.PI[1,]), col="black", lwd=1, lty=3)
lines(z.c_grid, (p.KG.PI[2,]), col="black", lwd=1, lty=3)

shade((p.NK.PI), z.c_grid, col=col.alpha("orange", 0.3))
lines(z.c_grid, (p.NK.PI[1,]), col="darkorange", lwd=1, lty=3)
lines(z.c_grid, (p.NK.PI[2,]), col="darkorange", lwd=1, lty=3)

points( z1 ~ SL1_mm, bg=col.alpha("gray", 1), pch=21, cex = 1,  data=toplot[toplot$KG == 1,])
points( z1 ~ SL1_mm, bg=col.alpha("orange", 1),  pch=21, data=toplot[toplot$NK==1,])

lines(z.c_grid, (p.KG.mean), col="black", lwd=2)
lines(z.c_grid, (p.NK.mean), col="darkorange2", lwd=2)


title(main='b)', adj = 0, line = 0.5)






pred.KG <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 1, NK = 0, NG = 0, center = 18))
pred.NK <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 1, NK = 1, NG = 0, center = 18))


mean(as.numeric(pred.NK / pred.KG)) * 100 -100
HPDI(as.numeric(pred.NK / pred.KG))  * 100 -100
LOS(as.numeric(pred.NK - pred.KG))

# Fecundity

z.c_grid <-  seq(from= x_axis_min, to=x_axis_max, length.out = 200)

pred.raw <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 1, NK = 0, NG = 0, center = 18))
p.KG.mean <- apply(pred.raw, 2, mean)
p.KG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

pred.raw <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 1, NK = 1, NG = 0, center = 18))
p.NK.mean <- apply(pred.raw, 2, mean)
p.NK.PI <- apply(pred.raw, 2, HPDI,prob = .95)

toplot <- subset(Gdata, Repr == 1)
toplot$Recr
plot( Recr ~ SL1_mm, data=toplot,  col=rangi2, ylab="Offsping (N)",
      xlab= "", pch="",  xlim=c(5,30), xaxt="n")
axis(1, labels = c(5,10,15,20,25,30), at = (c(5,10,15,20,25,30)) )
#axis(2, labels = c(0, 20, 40, 60, 80, 100), at = (c(0, 0.2, 0.4, 0.6, 0.8, 1)) )

shade((p.KG.PI), z.c_grid, col=col.alpha("black", 0.2))
lines(z.c_grid, (p.KG.PI[1,]), col="black", lwd=1, lty=3)
lines(z.c_grid, (p.KG.PI[2,]), col="black", lwd=1, lty=3)

shade((p.NK.PI), z.c_grid, col=col.alpha("orange", 0.3))
lines(z.c_grid, (p.NK.PI[1,]), col="darkorange", lwd=1, lty=3)
lines(z.c_grid, (p.NK.PI[2,]), col="darkorange", lwd=1, lty=3)

points( Recr ~ SL1_mm, bg=col.alpha("gray", 1), pch=21, cex = 1,  data=toplot[toplot$KG == 1,])
points( Recr ~ SL1_mm, bg=col.alpha("orange", 1),  pch=21, data=toplot[toplot$NK==1,])

lines(z.c_grid, (p.KG.mean), col="black", lwd=2)
lines(z.c_grid, (p.NK.mean), col="darkorange2", lwd=2)

title(main='c)', adj = 0, line = 0.5)


#legend('topleft',c('KG','NK'),lty=c(1,1),pch=c(16,16),col=c('black','orange'),bty='n')
z.c_grid <-  seq(from= 20, to=30, length.out = 50)
pred.KG <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 1, NK = 0, NG = 0, center = 18))
pred.NK <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 1, NK = 1, NG = 0, center = 18))

mean(as.numeric(pred.NK / pred.KG))
HPDI(as.numeric(pred.NK / pred.KG))
LOS(as.numeric(pred.NK - pred.KG))

mean(as.numeric(pred.NK / pred.KG)) * 100 -100
HPDI(as.numeric(pred.NK / pred.KG))  * 100 -100
LOS(as.numeric(pred.NK - pred.KG))

LOS(post$b_zNK_recrG)


# Killifish
z.c_grid <-  seq(from= 5, to=100, length.out = 200)

pred.raw <- sapply(1:length(z.c_grid), function(i) surv_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 0, center = 10))
p.KG.mean <- apply(pred.raw, 2, mean)
p.KG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

pred.raw <- sapply(1:length(z.c_grid), function(i) surv_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 1, center = 10))
p.NK.mean <- apply(pred.raw, 2, mean)
p.NK.PI <- apply(pred.raw, 2, HPDI,prob = .95)

toplot <- Kdata
toplot$SL1_mm
plot( surv ~ SL1_mm, data=toplot,  col=rangi2, ylab="Survival (%)",
      xlab= "Initial size (mm)", pch="",  xlim=c(5,100), xaxt="n",yaxt="n")
axis(1, labels = c(5,20,40, 60, 80, 100), at = c(5,20,40, 60, 80, 100) )
axis(2, labels = c(0, 20, 40, 60, 80, 100), at = (c(0, 0.2, 0.4, 0.6, 0.8, 1)) )


shade((p.KG.PI), z.c_grid, col=col.alpha("black", 0.3))
lines(z.c_grid, (p.KG.PI[1,]), col="black", lwd=1, lty=3)
lines(z.c_grid, (p.KG.PI[2,]), col="black", lwd=1, lty=3)

shade((p.NK.PI), z.c_grid, col=col.alpha("tomato", 0.2))
lines(z.c_grid, (p.NK.PI[1,]), col="tomato", lwd=1, lty=3)
lines(z.c_grid, (p.NK.PI[2,]), col="tomato", lwd=1, lty=3)



points( surv ~ SL1_mm, bg=col.alpha("black", 0.2), pch=21, cex = 0.5,  data=toplot[toplot$KG == 1,])
points( surv ~ SL1_mm, bg="tomato", pch=21, data=toplot[toplot$NG==1,])
lines(z.c_grid, (p.KG.mean), col="black", lwd=2)
lines(z.c_grid, (p.NK.mean), col="tomato", lwd=2)

title(main='d)', adj = 0, line = 0.5)
legend(x= 40, y = 0.3,c('Control (KG)','No guppy (NG)'),lty=c(1,1),pch=c(21,21), col=c(col.alpha("black", 1),'tomato'), pt.bg = c(col.alpha("gray", 1),'tomato'),bty='n')

pred.KG <- sapply(1:length(z.c_grid), function(i) surv_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 0, center = 10))
pred.NG <- sapply(1:length(z.c_grid), function(i) surv_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 1, center = 10))


LOS(pred.NG-pred.KG)

HPDI(pred.NG/pred.KG, prob=.95)


# Growth

pred.raw <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 0, center = 18))
p.KG.mean <- apply(pred.raw, 2, mean)
p.KG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

pred.raw <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 1, center = 18))
p.NG.mean <- apply(pred.raw, 2, mean)
p.NG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

toplot <- subset(Kdata, surv == 1)
toplot$z1 = toplot$z1 - toplot$SL1_mm 
plot( z1 ~ SL1_mm, data=toplot,  col=rangi2, ylab="Final size (mm)",
      xlab= "Initial size (mm)", pch="",  xlim=c(5,100), ylim=c(-6,10), xaxt="n")#, yaxt='n')
axis(1, labels = c(5,20,40, 60, 80, 100), at = c(5,20,40, 60, 80, 100) )
#axis(2, labels = c(5,20,40, 60, 80, 100), at = c(5,20,40, 60, 80, 100) )
#axis(2, labels = c(0, 20, 40, 60, 80, 100), at = (c(0, 0.2, 0.4, 0.6, 0.8, 1)) )


shade((p.KG.PI), z.c_grid, col=col.alpha("black", 0.3))
lines(z.c_grid, (p.KG.PI[1,]), col="black", lwd=1, lty=3)
lines(z.c_grid, (p.KG.PI[2,]), col="black", lwd=1, lty=3)

shade((p.NG.PI), z.c_grid, col=col.alpha("tomato", 0.2))
lines(z.c_grid, (p.NG.PI[1,]), col="tomato", lwd=1, lty=3)
lines(z.c_grid, (p.NG.PI[2,]), col="tomato", lwd=1, lty=3)

points( z1 ~ SL1_mm, bg=col.alpha("black", 0.2), pch=21, cex = 1,  data=toplot[toplot$KG == 1,])
points( z1 ~ SL1_mm, bg="tomato", pch=21, data=toplot[toplot$NG==1,])
lines(z.c_grid, (p.KG.mean), col="black", lwd=2)
lines(z.c_grid, (p.NG.mean), col="tomato", lwd=2)
title(main='e)', adj = 0, line = 0.5)




pred.KG <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 0, center = 10))

pred.NG <- sapply(1:length(z.c_grid), function(i) grow_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 1, center = 10))
mean(pred.NG/pred.KG)

LOS(pred.NG-pred.KG)

HPDI(pred.NG/pred.KG, prob=.95)

# Fecundity


pred.raw <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 0, center = 18))
p.KG.mean <- apply(pred.raw, 2, mean)
p.KG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

pred.raw <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 0, NK =0, NG = 1, center = 18))
p.NG.mean <- apply(pred.raw, 2, mean)
p.NG.PI <- apply(pred.raw, 2, HPDI,prob = .95)

toplot <- subset(Kdata, Repr == 1)
toplot$Recr
plot( Recr ~ SL1_mm, data=toplot,  col=rangi2, ylab="Offsping (N)",
      xlab= "Initial size (mm)", pch="",  xlim=c(5,100), xaxt="n")
axis(1, labels = c(5,20,40, 60, 80, 100), at = c(5,20,40, 60, 80, 100) )

shade((p.KG.PI), z.c_grid, col=col.alpha("black", 0.3))
lines(z.c_grid, (p.KG.PI[1,]), col="black", lwd=1, lty=3)
lines(z.c_grid, (p.KG.PI[2,]), col="black", lwd=1, lty=3)

shade((p.NG.PI), z.c_grid, col=col.alpha("tomato", 0.2))
lines(z.c_grid, (p.NG.PI[1,]), col="tomato", lwd=1, lty=3)
lines(z.c_grid, (p.NG.PI[2,]), col="tomato", lwd=1, lty=3)

points( Recr ~ SL1_mm, bg=col.alpha("black", 0.2), pch=21, cex = 1,  data=toplot[toplot$KG == 1,])
points( Recr ~ SL1_mm, bg="tomato", pch=21, data=toplot[toplot$NG==1,])

lines(z.c_grid, (p.KG.mean), col="black", lwd=2)
lines(z.c_grid, (p.NG.mean), col="tomato", lwd=2)

title(main='f)', adj = 0, line = 0.5)


z.c_grid <-  seq(from= 20, to=100, length.out = 200)
pred.KG <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 0, NK = 0, NG = 0, center = 10))
pred.NG <- sapply(1:length(z.c_grid), function(i) recr_link(post = post, size = z.c_grid[i], G = 0, NK =0, NG = 1, center = 10))


mean(as.numeric(pred.NG / pred.KG))
HPDI(as.numeric(pred.NG / pred.KG))
LOS(as.numeric(pred.NG - pred.KG))

LOS(post$b_zNK_recrG)
graphics.off()


################################################################################################################################################################


jpeg(file = "Fitness.jpeg", width = 8, height = 7, units='in', res=800)
op<-par(mfrow=c(2,2), mar = c(4, 4 , 2, 1), oma = c(0.5, 1, 1, 0.5)) # c(bottom, left, top, right)
toPlot <- rbind(Gsum[c("lamda.GR","lamda.NR"),], Ksum[c("lamda.GR","lamda.NG"),])
toPlot$x <- 1:4
library(latex2exp)
plot(mean ~ x, toPlot, xlim=c(0.5,4.5), xaxt="n", xlab = "", ylim = c(0.8, 2.5), pch="",
     ylab = TeX("Absolute Fitness ($\\lambda$)"))
segments(x0 = c(1,2,3,4), x1 = c(1,2,3,4), y0 = toPlot$`2.5%`, y1 = toPlot$`97.5%`, lwd = 1.5)
segments(x0 = c(1), x1 = c(2), y0 = toPlot$mean[1], y1 = toPlot$mean[2], lwd = 1, lty=3)
segments(x0 = c(3), x1 = c(4), y0 = toPlot$mean[3], y1 = toPlot$mean[4], lwd = 1, lty=3)

points(mean ~ x, toPlot, pch= c(21,21,22,22), bg = c("gray") , cex = 2)
abline(h=1, lwd = 1.25, lty=2, col="gray")
axis(1, at = c(1,2,3,4), labels = c("KG", "NK", "KG", "NG"))
abline(v=2.5, lwd = 1.5, col='black')
legend("topright", c('Guppy', "Killifish"), pch= c(21, 22), pt.bg =  "gray", pt.cex = 1.75, bty='n')
#graphics.off()
title(main='a) Overall fitness effects', adj = 0, line = 0.5)

###-------

toPlot <- rbind(Gsum[c("Sum.lam.effect", "lam.surv","lam.growth", "lam.fec", "lam.rcz"),], Ksum[c("Sum.lam.effect", "lam.surv","lam.growth", "lam.fec", "lam.rcz"),] )
x = c(1,1.2,1.4,1.6,1.8, 2.2,2.4,2.6,2.8,3 )
toPlot$x <- x
library(latex2exp)
plot(mean ~ x, toPlot, xlim=c(0.8,3.1), ylim=c(-0.7,0.6), xaxt="n", xlab = "",  pch="",
     ylab = TeX("Contributions to relative fitness ($\\Delta\\lambda$)"))

library('plotfunctions')
add_bars(x = toPlot$x, y = toPlot$mean, width = 0.18, col= c("gray", 'tomato', 'cyan4', 'purple', 'black') )
segments(x0 = x, x1 = x, y0 = toPlot$`2.5%`, y1 = toPlot$`97.5%`, lwd = 1.5)

abline(h=0, lwd = 1.25, lty=2, col="gray")
axis(1, at = c(1.3,2.3), labels = c("Guppy", "Killifish"))
abline(v=1.95, lwd = 1.5, col='black')
#abline(v=2.6, lwd = 1.5, col='black')
text(x = 1.5, y = 4, labels = "Guppy")
text(x = 3.5, y = 4, labels = "Killifish")
title(main='b) LTRE decompositions', adj = 0, line = 0.5)
legend("bottomright", c('All vital rates','Survival', "Growth", "Fecundity", "Offspring size"), pch=c(22), pt.bg = c("gray", 'tomato', 'cyan4', 'purple', 'black'), pt.cex = 2.5,  bty='n', 
       y.intersp = 0.85)



##-------

# Guppy
nBigMatrix <- 100
# Make the meshpoints
min.size <- (2) #min(IPMdata$z1, na.rm = T) - sd(IPMdata$z1, na.rm = T)*0.5
max.size <- (40) #max(IPMdata$z1, na.rm = T) + sd(IPMdata$z1, na.rm = T)*0.5

U=max.size
L=min.size
m = nBigMatrix
h <- (U - L)/m
meshpts <- z1 <- z <- L + ((1:m) - 1/2) * h
size.cen <- (18)
mean <- apply(Gsurv.mat,2, mean)
ci <- apply(Gsurv.mat,2, HPDI, prob=.95)
meshpts[t(ci)[,1] >0]
meshpts[t(ci)[,2] < 0]

plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.03,0.005), xlim = c(5,30), xlab= " Guppy size (mm)")
abline(h=0, lty=1, col="gray")

# survival
mean <- apply(Gsurv.mat,2, mean)
ci <- apply(Gsurv.mat,2, HPDI, prob=.95)
range(meshpts[which(apply(Gsurv.mat, 2, LOS) > 95)])
lines(mean ~ meshpts, lwd=2, col = "tomato", lty=3)
#shade(ci, meshpts, col = col.alpha('tomato', 0.2))
a = range(meshpts[which(apply(Gsurv.mat, 2, LOS) > 95)])
abline(v=a[1], lty=3, col='red' )
abline(v=a[2], lty=3, col='red' )

# Growth
mean <- apply(Ggrow.mat,2, mean)
ci <- apply(Ggrow.mat,2, HPDI, prob=.95)
range(meshpts[which(apply(Ggrow.mat, 2, LOS) > 95)])
lines(mean ~ meshpts, lwd=2, col = "cyan4", lty= 3)
lines(mean[which(apply(Ggrow.mat, 2, LOS) > 95)] ~ meshpts[which(apply(Ggrow.mat, 2, LOS) > 95)], lwd=3, col = "cyan4", lty= 1)
#shade(ci, meshpts, col = col.alpha('cyan4', 0.2))

# Fecundity
mean <- apply(Gfec.mat,2, mean)
ci <- apply(Gfec.mat,2, HPDI, prob=.95)
range(meshpts[which(apply(Gfec.mat, 2, LOS) > 95)])
lines(mean ~ meshpts, lwd=1.5, col = "purple", lty = 3)
lines(mean[which(apply(Gfec.mat, 2, LOS) > 95)] ~ meshpts[which(apply(Gfec.mat, 2, LOS) > 95)], lwd=3, col = "purple", lty= 1)
#shade(ci, meshpts, col = col.alpha('purple', 0.2))
# a = range(meshpts[which(apply(Gfec.mat, 2, LOS) > 95)])
# abline(v=a[1], lty=3, col='purple' )
# abline(v=a[2], lty=3, col='purple' )

title(main="c) Guppy", adj = 0, line = 0.5)
legend("topleft", c('All vital rates','Survival', "Growth", "Fecundity"), pch=c(22), pt.bg = c("gray", 'tomato', 'cyan4', 'purple', 'black'), pt.cex = 2.5,  bty='n', 
       y.intersp = 0.5, inset = 1)


#legend("bottomright", c('Survival', "Growth", "Fecundity"), lty= 1, col = c('tomato', 'cyan4', 'purple', 'black'), lwd = 2,  bty='n')
#####----------------------------------------------------------
nBigMatrix <- 100

# Make the meshpoints
min.size <- (2) 
max.size <- (110) 

U=max.size
L=min.size
m = nBigMatrix
h <- (U - L)/m
meshpts <- z1 <- z <- L + ((1:m) - 1/2) * h
size.cen <- (18)

mean <- apply(Ksurv.mat,2, mean)
ci <- apply(Ksurv.mat,2, HPDI, prob=.95)
meshpts[t(ci)[,1] >0]
meshpts[t(ci)[,2] < 0]

meshpts <- c(0,meshpts)
# Killifish


plot(mean ~ meshpts, ylab = TeX("Fitness effects ($\\Delta \\lambda$)"), pch="", ylim = c(-0.0001,0.02), xlim = c(5,100), xlab= " Killifish size (mm)")
abline(h=0, lty=1, col="gray")

# survival
mean <- apply(Ksurv.mat,2, mean)
ci <- apply(Ksurv.mat,2, HPDI, prob=.95)
range(meshpts[which(apply(Ksurv.mat, 2, LOS) > 95)])
lines(mean ~ meshpts, lwd=2, col = "tomato", lty=3)
#shade(ci, meshpts, col = col.alpha('tomato', 0.2))
a = range(meshpts[which(apply(Ksurv.mat, 2, LOS) > 95)])
abline(v=a[1], lty=3, col='red' )
abline(v=a[2], lty=3, col='red' )

# Growth
mean <- apply(Kgrow.mat,2, mean)
ci <- apply(Kgrow.mat,2, HPDI, prob=.95)
range(meshpts[which(apply(Kgrow.mat, 2, LOS) > 95)])
lines(mean ~ meshpts, lwd=2, col = "cyan4", lty= 3)
lines(mean[which(apply(Kgrow.mat, 2, LOS) > 95)] ~ meshpts[which(apply(Kgrow.mat, 2, LOS) > 95)], lwd=3, col = "cyan4", lty= 1)
#shade(ci, meshpts, col = col.alpha('cyan4', 0.2))

# Fecundity
mean <- apply(Kfec.mat,2, mean)
ci <- apply(Kfec.mat,2, HPDI, prob=.95)
range(meshpts[which(apply(Kfec.mat, 2, LOS) > 95)])
lines(mean ~ meshpts, lwd=1.5, col = "purple", lty = 3)
lines(mean[which(apply(Kfec.mat, 2, LOS) > 95)] ~ meshpts[which(apply(Kfec.mat, 2, LOS) > 95)], lwd=3, col = "purple", lty= 1)
#shade(ci, meshpts, col = col.alpha('purple', 0.2))
# a = range(meshpts[which(apply(Kfec.mat, 2, LOS) > 95)])
# abline(v=a[1], lty=3, col='purple' )
# abline(v=a[2], lty=3, col='purple' )

title(main="d) Killifish", adj = 0, line = 0.5)


#legend("topright", c('Survival', "Growth", "Fecundity"), lty= 1, col = c('tomato', 'cyan4', 'purple', 'black'), lwd = 2,  bty='n')



graphics.off()

