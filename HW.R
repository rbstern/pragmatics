library(DirichletReg)

posterior <- function(data,prior_par=c(1,1,1))
{
  posterior_par <- data+prior_par
  return(posterior_par)
}


log_posterior <- function(x,posterior_par)
{
  return(ddirichlet(x, posterior_par, log = TRUE))
}

hpd_th <- function(data,prior_par=c(1,1,1),probability=0.95,B=100000)
{
  posterior_par <- posterior(data=data,
                             prior_par=prior_par)
  
  sample_posterior <- rdirichlet(n=B,alpha=posterior_par)
  sample_posterior_densities <- log_posterior(sample_posterior,
                                              posterior_par=posterior_par)
  
  th <- quantile(sample_posterior_densities,probs=1-probability)
  return(th)
}

numerator <- function(theta,p,m=1)
{
  return(m*((theta[1]-p^2)^2+(theta[2]-2*p*(1-p))^2+(theta[3]-(1-p)^2)^2))
}
denominator <- function(p)
{
  return(p^2*(1-p^2)+2*p*(1-p)*(1-2*p*(1-p))+(1-p)^2*(1-(1-p)^2))
}
lambda <- function(theta)
{
  num <- sqrt(theta[1]*theta[3])-theta[2]/2
  den <- sqrt(theta[1]*theta[3])+theta[2]/2
  return(num/den)
}
distance_lambda <- function(theta,p,m=1)
{
  return(abs(lambda(theta)))
}
distance_estimator <- function(theta,p,m=1)
{
  return(numerator(theta,p,m=1)/denominator(p))
}

add_points <- function(p0,B,eps,col=colors()[colors()=="dodgerblue4"],
                       distance=distance_estimator,cex=0.25)
{
  for(x in seq(0.0000001,0.99999999,length.out = B))
  {
    for(y in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(distance(theta=c(x,theta2,y),p=p0,m=1)<eps)
        points(x,y,pch=16,cex=cex,col=col)
    }
  }
}

pdf("figures/hw.pdf")
# Plot 1
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=2)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=2)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)

p0 <- 1/3
eps=0.1/10
B=200
add_points(p0,B,eps)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p0^2,(1-p0)^2,cex=1.5,pch=16,col=1)


# Plot 2
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=2)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=2)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)


eps=0.1/10
B=200
for(p0 in seq(0.0001,0.9999,length.out = 100))
{
  add_points(p0,B,eps)
}
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)

# Plot 3
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=2)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=2)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)


eps=0.1/10
B=200
for(p0 in seq(0.0001,0.9999,length.out = 100))
{
  add_points(p0,B,eps)
}
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)

data <- c(55,83,50)
th <- hpd_th(data)
B=100
for(x in seq(0.0000001,0.99999999,length.out = B))
{
  for(y in seq(0.0000001,0.99999999,length.out = B))
  {
    if(x+y>1)
      next;
    theta2 <- 1-x-y
    if(theta2>1)
      next;
    if(log_posterior(x=t(c(x,theta2,y)),
                     posterior_par = posterior(data))>th)
      points(x,y,pch=16,cex=0.25,col=colors()[colors()=="green4"])
  }
}
post <- posterior(data)-1
text(x=(post/sum(post))[1]+eps*3,
     y=(post/sum(post))[3]+eps*3,label="Case",cex=1.5)

data <- c(24,42,39)
th <- hpd_th(data)
B=1000
for(x in seq(0.0000001,0.99999999,length.out = B))
{
  for(y in seq(0.0000001,0.99999999,length.out = B))
  {
    if(x+y>1)
      next;
    theta2 <- 1-x-y
    if(theta2>1)
      next;
    if(log_posterior(x=t(c(x,theta2,y)),
                     posterior_par = posterior(data))>th)
      points(x,y,pch=16,cex=0.25,col=colors()[colors()=="green4"])
  }
}
post <- posterior(data)
text(x=(post/sum(post))[1]+eps*5,
     y=(post/sum(post))[3]+eps*5,label="Control",cex=1.5)



# Plot 4
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=2)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=2)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)

eps=0.1/10
B=250
for(p0 in seq(0.0001,0.9999,length.out = 100))
{
  add_points(p0,B,eps,cex=0.15)
}
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)

thetaAA <- 0.2
thetaAB <- 0.3
thetaBB <- 1-thetaAA-thetaAB
set.seed(1)
sample_new=sample(1:3,1000,prob = c(thetaAA,thetaAB,thetaBB),replace = TRUE)
sample_new <- table(sample_new)
p <- 0.8
thetaAA <- p^2
thetaAB <- 2*p*(1-p)
thetaBB <- (1-p)^2
sample_new2=sample(1:3,1000,prob = c(thetaAA,thetaAB,thetaBB),replace = TRUE)
sample_new2 <- table(sample_new2)
data <- matrix(c(4,6,57,58,120,206,110,34,
                 18,53,118,97,361,309,148,22,
                 94,74,100,48,194,142,44,12),ncol=3)
data <- rbind(data,sample_new,sample_new2)
for(experiment in 1:nrow(data))
{
  th <- hpd_th(data[experiment,],probability = 0.8)
  B=1200
  for(x in seq(0.0000001,0.99999999,length.out = B))
  {
    previous="small"
    for(y in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(log_posterior(x=t(c(x,theta2,y)),
                       posterior_par = posterior(data[experiment,]))>th)
      {
        points(x,y,pch=16,cex=0.15,col=colors()[colors()=="green2"])
        if(previous=="small")
          points(x,y,pch=16,cex=0.35,col=colors()[colors()=="blue2"])
        previous="large"
      } else if(previous=="large") {
        points(x,y,pch=16,cex=0.35,col=colors()[colors()=="blue2"])
        previous="small"
      }
    }
  }
  post <- posterior(data[experiment,])-1
  text(x=(post/sum(post))[1]-1*eps,
       y=(post/sum(post))[3]+0*eps,label=paste0("",experiment,""),cex=1)
}

for(experiment in 1:nrow(data))
{
  th <- hpd_th(data[experiment,],probability = 0.8)
  B=1000
  for(x in seq(0.0000001,0.99999999,length.out = B))
  {
    previous="small"
    for(y in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(log_posterior(x=t(c(x,theta2,y)),
                       posterior_par = posterior(data[experiment,]))>th)
      {
        #points(x,y,pch=16,cex=0.15,col=colors()[colors()=="green2"])
        if(previous=="small")
          points(x,y,pch=16,cex=0.35,col=colors()[colors()=="blue2"])
        previous="large"
      } else if(previous=="large") {
        points(x,y,pch=16,cex=0.35,col=colors()[colors()=="blue2"])
        previous="small"
      }
    }
  }
  post <- posterior(data[experiment,])-1
  text(x=(post/sum(post))[1]-1*eps,
       y=(post/sum(post))[3]+0*eps,label=paste0("",experiment,""),cex=1)
}
for(experiment in 1:nrow(data))
{
  th <- hpd_th(data[experiment,],probability = 0.8)
  B=1000
  for(y in seq(0.0000001,0.99999999,length.out = B))
  {
    previous="small"
    for(x in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(log_posterior(x=t(c(x,theta2,y)),
                       posterior_par = posterior(data[experiment,]))>th)
      {
        #points(x,y,pch=16,cex=0.15,col=colors()[colors()=="green2"])
        if(previous=="small")
          points(x,y,pch=16,cex=0.35,col=colors()[colors()=="blue2"])
        previous="large"
      } else if(previous=="large") {
        points(x,y,pch=16,cex=0.35,col=colors()[colors()=="blue2"])
        previous="small"
      }
    }
  }
  post <- posterior(data[experiment,])-1
  text(x=(post/sum(post))[1]-1*eps,
       y=(post/sum(post))[3]+0*eps,label=paste0("",experiment,""),cex=1)
}

dev.off()


png("figures/hw1.png",width = 4, height = 4, units = 'in', res = 400)
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=2)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=2)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)

p0 <- 1/3
eps=0.1/10
B=200
add_points(p0,B,eps,cex=0.12)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p0^2,(1-p0)^2,cex=1.5,pch=16,col=1)
dev.off()


png("figures/hw2.png",width = 4, height = 4, units = 'in', res = 400)
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=2)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=2)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)


eps=0.1/10
B=200
for(p0 in seq(0.0001,0.9999,length.out = 100))
{
  add_points(p0,B,eps,cex=0.12)
}
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
dev.off()


png("figures/hw3.png",width = 4, height = 4, units = 'in', res = 400)
plot(c(0,1),c(0,1),pch=NA,xlab="",
     ylab="",  yaxt="n",bty="l",
     xaxt="n",xlim=c(0,1),ylim = c(0,1),
     xaxs="i",yaxs="i",cex.lab=2)
title(ylab=expression(theta[3]), line=0.5, 
      cex.lab=1.5)
title(xlab=expression(theta[1]), line=1.5, 
      cex.lab=1.5)
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)
points(p,(1-p),pch=16,cex=0.4)

eps=0.1/10
B=250
for(p0 in seq(0.0001,0.9999,length.out = 100))
{
  add_points(p0,B,eps,cex=0.1)
}
p <- seq(0,1,length.out = 1000)
points(p^2,(1-p)^2,pch=16,cex=0.4,col=2)

thetaAA <- 0.2
thetaAB <- 0.3
thetaBB <- 1-thetaAA-thetaAB
set.seed(1)
sample_new=sample(1:3,1000,prob = c(thetaAA,thetaAB,thetaBB),replace = TRUE)
sample_new <- table(sample_new)
p <- 0.8
thetaAA <- p^2
thetaAB <- 2*p*(1-p)
thetaBB <- (1-p)^2
sample_new2=sample(1:3,1000,prob = c(thetaAA,thetaAB,thetaBB),replace = TRUE)
sample_new2 <- table(sample_new2)
data <- matrix(c(4,6,57,58,120,206,110,34,
                 18,53,118,97,361,309,148,22,
                 94,74,100,48,194,142,44,12),ncol=3)
data <- rbind(data,sample_new,sample_new2)
for(experiment in 1:nrow(data))
{
  th <- hpd_th(data[experiment,],probability = 0.8)
  B=1200
  for(x in seq(0.0000001,0.99999999,length.out = B))
  {
    previous="small"
    for(y in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(log_posterior(x=t(c(x,theta2,y)),
                       posterior_par = posterior(data[experiment,]))>th)
      {
        points(x,y,pch=16,cex=0.1,col=colors()[colors()=="green2"])
        if(previous=="small")
          points(x,y,pch=16,cex=0.3,col=colors()[colors()=="blue2"])
        previous="large"
      } else if(previous=="large") {
        points(x,y,pch=16,cex=0.3,col=colors()[colors()=="blue2"])
        previous="small"
      }
    }
  }
}

for(experiment in 1:nrow(data))
{
  th <- hpd_th(data[experiment,],probability = 0.8)
  B=1000
  for(x in seq(0.0000001,0.99999999,length.out = B))
  {
    previous="small"
    for(y in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(log_posterior(x=t(c(x,theta2,y)),
                       posterior_par = posterior(data[experiment,]))>th)
      {
        #points(x,y,pch=16,cex=0.15,col=colors()[colors()=="green2"])
        if(previous=="small")
          points(x,y,pch=16,cex=0.3,col=colors()[colors()=="blue2"])
        previous="large"
      } else if(previous=="large") {
        points(x,y,pch=16,cex=0.3,col=colors()[colors()=="blue2"])
        previous="small"
      }
    }
  }
}
for(experiment in 1:nrow(data))
{
  th <- hpd_th(data[experiment,],probability = 0.8)
  B=1000
  for(y in seq(0.0000001,0.99999999,length.out = B))
  {
    previous="small"
    for(x in seq(0.0000001,0.99999999,length.out = B))
    {
      if(x+y>1)
        next;
      theta2 <- 1-x-y
      if(theta2>1)
        next;
      if(log_posterior(x=t(c(x,theta2,y)),
                       posterior_par = posterior(data[experiment,]))>th)
      {
        #points(x,y,pch=16,cex=0.15,col=colors()[colors()=="green2"])
        if(previous=="small")
          points(x,y,pch=16,cex=0.3,col=colors()[colors()=="blue2"])
        previous="large"
      } else if(previous=="large") {
        points(x,y,pch=16,cex=0.3,col=colors()[colors()=="blue2"])
        previous="small"
      }
    }
  }
  post <- posterior(data[experiment,])-1
  text(x=(post/sum(post))[1]-1*eps,
       y=(post/sum(post))[3]+1*eps,label=paste0("",experiment,""),cex=0.7)
}
dev.off()
