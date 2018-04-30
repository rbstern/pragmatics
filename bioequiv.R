png("figures/bioequiv1.png",width = 4, height = 4, units = 'in', res = 400)
# Plot 1
plot(c(-1.05,1.05),c(-1.05,1.05),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="n",
     xaxt="n",cex.lab=2)
p <- seq(-1,1,length.out = 1000)
eps_new=0.15
points(p,p,pch=16,cex=0.4,col=2,lwd=3)
axis(1, pos=-0.8, at=-5:5,labels = FALSE, lwd.tick=0,xlab="a")
axis(2, pos=-0.8,labels = FALSE, lwd.tick=0)
text(0.85,-0.95,expression(mu[1]), 
     cex=1)
text(-0.90,0.85,expression(mu[2]), 
     cex=1)
p0=0
points(p0,p0,cex=1.1,pch=16,col=1)
B <- 50
p_circle <- seq(p0-eps_new,p0+eps_new,length.out = B)
for(i in seq_along(p_circle))
{
  for(j in seq_along(p_circle))
  {
    if((p_circle[i])^2+(p_circle[j])^2<eps_new^2)
    {
      points(p_circle[i],p_circle[j],
             cex=0.08,pch=16,col=colors()[colors()=="dodgerblue4"])   
    }
  }
}
points(p0,p0,cex=1.1,pch=16,col=1)
lines(x=c(0,0),y=c(-0.8,0))
lines(x=c(0,-0.8),y=c(0,0))
text(0,-0.9,expression(mu[0]), 
     cex=1)
text(-0.9,0,expression(mu[0]), 
     cex=1)
dev.off()


png("figures/bioequiv2.png",width = 4, height = 4, units = 'in', res = 400)
plot(c(-1.05,1.05),c(-1.05,1.05),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="n",
     xaxt="n",cex.lab=2)
p <- seq(-1,1,length.out = 1000)
eps_new=0.15
points(p,p,pch=16,cex=0.4,col=2,lwd=3)
axis(1, pos=-0.8, at=-5:5,labels = FALSE, lwd.tick=0,xlab="a")
axis(2, pos=-0.8,labels = FALSE, lwd.tick=0)
text(0.85,-0.95,expression(mu[1]), 
     cex=1)
text(-0.90,0.85,expression(mu[2]), 
     cex=1)

points(p,p+eps_new,pch=16,cex=0.2,col=colors()[colors()=="dodgerblue4"],lwd=1)
points(p,p-eps_new,pch=16,cex=0.2,col=colors()[colors()=="dodgerblue4"],lwd=1)
  
B <- 10
for(i in seq_along(p))
{
  points(rep(p[i],B),seq(p[i]-eps_new,p[i]+eps_new,length.out = B),
         cex=0.08,pch=16,col=colors()[colors()=="dodgerblue4"])
}

points(-0.8,eps_new-0.8,pch=16)
text(-0.9,eps_new-0.75,expression(epsilon*sigma), 
     cex=1)
points(-0.8,-eps_new-0.8,pch=16)
text(-0.65,-eps_new-0.85,expression(-epsilon*sigma), 
     cex=1)
points(p,p,pch=16,cex=0.4,col=2,lwd=3)

dev.off()