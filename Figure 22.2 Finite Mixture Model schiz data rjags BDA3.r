rm(list=ls(all=TRUE))

require(BayesDA)
library(ggplot2)
library(tidyverse)
library(tidybayes)
require(rjags)


data(schiz)
y<-log(schiz)
S=c(rep(0,11),rep(1,6))

mix_model <- '
model{
for(j in 1:17){
       alpha[j]~dnorm(mu+beta*S[j],1/sigma_a^2)
  for(i in 1:30){
       z[i,j]~dbern(lambda*S[j])
       y[i,j]~dnorm(alpha[j]+tau*z[i,j],1/sigma_y^2)
       y_rep[i,j]~dnorm(alpha[j]+tau*z[i,j],1/sigma_y^2)
     }
}
  lambda~dunif(0,1)
  beta~dnorm(0,1)
  tau~dnorm(0,1)T(0,)#identify the model
  mu~dnorm(0,1)
  sigma_a ~  dt(0,10^-2,1)T(0,)#error term
  sigma_y ~  dt(0,10^-2,1)T(0,)#error term
}'


data_list<-list(y = y, S = S  )
                  

##parameters to save
jags_pars <- c("lambda","beta","tau","y_rep") #
 
model <- jags.model(textConnection(mix_model), data=data_list, n.chains = 5)
    
update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                    variable.names = jags_pars ,                   
                    n.iter = 200, 
                    progress.bar = "none")

par_dat <- posterior_sample %>% spread_draws(lambda,beta,tau)#tidybayes

#Table 22.1
summary(par_dat)

###########################################################################3
#Figure 22.2 BDA3
ppdata=posterior_sample %>%  spread_draws(y_rep[i,v] |i)%>%filter(v>11)
#print(ppdata) 

sd_matrix=matrix(apply(ppdata[,5:34],1,sd),nrow=1000)

T_min=apply(sd_matrix,1,min)
T_max=apply(sd_matrix,1,max)
T=apply(y[,12:17],2,sd)

plot(T_min,T_max,type="p",main="Figure 22.2 BDA3",pch = 20,cex=1,xlim=c(.1,.4),ylim=c(.25,.6))
points(min(T),max(T),pch = 4,cex=2)


