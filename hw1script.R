####
####  R Script to simulate data and fit normal-normal-inverse-gamma model 
####

set.seed(2021)
#mu.true=3
#s2.true=2^2
#y=rnorm(10,mu.true,sqrt(s2.true))


df = read.csv("/Users/sergi/Documents/Bayesian Statistical Methods/hw1/NFL_pos_spending_2022.csv")
differences = df[["QB"]]-df[["WR"]]
y=differences

plot(df[["QB"]]-df[["WR"]])

mu.0=0
s2.0=100
source("norm.invchisq.mcmc.R")
n.mcmc=100000
n.burn=round(.1*n.mcmc)
out=norm.invchisq.mcmc(y,mu.0,s2.0,1,mean(y),n.mcmc)

####
####  Trace Plots
####

layout(matrix(1:2,2,1))
plot(out$mu.save,type="l",main="",ylab=bquote(mu),xlab="MCMC Iteration")
abline(h=mu.true,lwd=2,col=rgb(0,1,0,.5))
plot(out$s2.save,type="l",main="",ylab=bquote(sigma^2),xlab="MCMC Iteration")
abline(h=s2.true,lwd=2,col=rgb(0,1,0,.5))

####
####  Posterior Means and CIs 
####

mean(out$mu.save[-(1:n.burn)])
quantile(out$mu.save[-(1:n.burn)],c(0.025,0.975))

mean(out$s2.save[-(1:n.burn)])
quantile(out$s2.save[-(1:n.burn)],c(0.025,0.975))

####
####  Posterior Histograms 
####

dinvchisq <- function(x,nu){
  (2^(-nu/2) / gamma(nu/2)) * x^(-nu/2-1) * exp(-1/(2*x))
}

d = density(out$mu.save[n.burn:n.mcmc])
xx = d$x 
dx = xx[2L] - xx[1L] 
yy = d$y 
result = sum(yy[xx >= 0]) * dx
result

layout(matrix(1:1,1,1))
plot(density(out$mu.save[n.burn:n.mcmc]),lwd=2,main="Posterior and Prior: mu",add=TRUE)
curve(dnorm(x,mu.0,sqrt(s2.0)),col=2,lwd=2,add=TRUE)
abline(v=mean(y),col=4)
abline(v=quantile(out$mu.save[-(1:n.burn)],c(0.025,0.975)),col=rgb(1,0,0,.5),lwd=3,lty=2)
legend("topright",lty=c(1,1,1,2),col=c(1,2,4,2),lwd=2,legend=c("Posterior","Prior","Sample","95% ETI"))
plot(density(out$s2.save[n.burn:n.mcmc]),lwd=2,main="Posterior and Prior: s2")
curve(dinvchisq(x,out$nu),col=2,lwd=2,add=TRUE)
abline(v=var(y),col=4)
var(y)
abline(v=quantile(out$s2.save[-(1:n.burn)],c(0.025,0.975)),col=rgb(1,0,0,.5),lwd=3,lty=2)
legend("topright",lty=c(1,1,1,2),col=c(1,2,4,2),lwd=2,legend=c("Posterior","Prior","Sample","95% ETI"))

