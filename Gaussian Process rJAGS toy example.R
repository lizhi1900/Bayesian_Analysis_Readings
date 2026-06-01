#https://www.niamhcahill.com/post/gp_tutorial/
#https://hendersontrent.github.io/posts/2024/05/gaussian-process-time-series/
rm(list=ls(all=TRUE))
library(tidyverse)
library(tidybayes)
require(MASS)
require(rjags)


set.seed(123)
y <- 3 * sin(2 * seq(0, 4 * pi, length.out = 100)) + runif(100) * 2
trend <- 0.08 * seq(from = 1, to = 100, by = 1)
y <- y + trend
tmp <- data.frame(timepoint = 1:length(y), y = y)

ggplot(data = tmp) +
  geom_point(aes(x = timepoint, y = y), colour = "black") +
  geom_line(aes(x = timepoint, y = y), colour = "black") +
  labs(x = "Timepoint",
       y = "Values") +
  theme_bw()


###########################################################
gp_model <- '
model{

  gp ~ dmnorm(mu,Sigma.inv)
  Sigma.inv <- inverse(Sigma)

  #off diagonals
  for(i in 1:(n_obs-1))
  {
    for(j in (i+1):n_obs) {
     cov_exp_quad[i,j] <- sigma_1^2*exp(-(1/l_1^2)*(d[i,j]^2))
     cov_periodic[i,j] <- sigma_2^2*exp(-(2/l_2^2)*sin(pi*d[i,j]/p)^2)
     cov_exp_quad[j,i] <- cov_exp_quad[i,j]
     cov_periodic[j,i] <- cov_periodic[i,j]
     Sigma[i,j] <- cov_exp_quad[i,j]+cov_periodic[i,j]
     Sigma[j,i] <- Sigma[i,j]
    }
  }
  #diagonals
   for(i in 1:n_obs){
    mu[i] <- 0
    cov_exp_quad[i,i] <- sigma_1^2 + 0.00001
    cov_periodic[i,i] <- sigma_2^2 
    Sigma[i,i] <- cov_exp_quad[i,i]+cov_periodic[i,i]
    y[i]      ~ dnorm(gp[i],sigma_y^-2)
    y_pred[i] ~ dnorm(gp[i],sigma_y^-2)
   } 
  p ~ dt(0,10^-2,1)T(0,)
  sigma_1 ~ dt(0,10^-2,1)T(0,)
  sigma_2 ~ dt(0,10^-2,1)T(0,)
  l_1~ dt(0,4^-2,1)T(0,)
  l_2 ~ dt(0,4^-2,1)T(0,)
  sigma_y ~ dt(0,10^-2,1)T(0,)
}'

### get data and estimation years
x <- tmp%>% pull(timepoint)
y <- tmp%>% pull(y)
n_obs <- length(x)
d <- abs(outer(x,x,"-"))
 

data_list<-list(y = y,pi=pi,
                n_obs=n_obs,
                d = d)  

##parameters to save
jags_pars <- c("y_pred") 
 
model <- jags.model(textConnection(gp_model), data=data_list, n.chains = 5)
    
update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                    variable.names = jags_pars ,                   
                    n.iter = 400, 
                    progress.bar = "none")

out=summary(posterior_sample)
attributes(out)
lower=out$quantiles[,1]
upper=out$quantiles[,5]
dat <- data.frame(timepoint = 1:length(y), y = y,lower=lower,upper=upper)

ggplot(data = dat,aes(x, y)) +
  geom_point(aes(x = timepoint, y = y), colour = "black") +
  geom_line(aes(x = timepoint, y = y), colour = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3,fill = "blue") +
  labs(x = "Timepoint",
       y = "Values") +
  theme_bw()

####################################################


