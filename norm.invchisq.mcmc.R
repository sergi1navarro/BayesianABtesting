norm.invchisq.mcmc <- function(y,mu.0,s2.0,nu,mu.strt,n.mcmc){
  
#
#  Gibbs Sampler for Gaussian Data Model, with Gaussian prior on mu and invchisq on s2.
#

####
####  Set up variables 
####
  
n.burn=round(n.mcmc/10)
n=length(y)
mu.save=rep(0,n.mcmc)
s2.save=rep(0,n.mcmc)
mu=mu.strt
  
#nu = 1
  
####
####  Begin Gibbs Loop 
####
  
for(k in 1:n.mcmc){
  if((k%%1000)==0) cat(k," ") # this prints k every 1000 iterations
    
  ####
  ####  Sample s2 
  ####
    
  tmp.nu=n+nu

  s2=1/rchisq(1, tmp.nu, 1)

  ####
  ####  Sample mu 
  ####
    
  tmp.mn=(s2*mu.0+s2.0*sum(y))/(s2+n*s2.0)
  tmp.var=s2*s2.0/(s2+n*s2.0)

  mu=rnorm(1,tmp.mn,sqrt(tmp.var))
    
  ####
  ####  Save Samples 
  ####
    
  mu.save[k]=mu
  s2.save[k]=s2
    
  }
  cat("\n")
  
  ####
  ####  Write Output 
  ####
  
  list(mu.save=mu.save,s2.save=s2.save,mu.0=mu.0,s2.0=s2.0,nu=nu,n.mcmc=n.mcmc,n.burn=n.burn)
  
}