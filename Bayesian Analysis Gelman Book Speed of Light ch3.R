

#Reference: Bayesian Data Analysis By Andrew Gelman, 2003
#zli3  [at] live.com


rm(list=ls(all=TRUE))

setwd("C:\\tuDrive2025\\bayesian\\speed_light\\")
mydata=read.csv("light.asc",skip=4,header=0,sep=" ")


mydata=as.vector(unlist(t(mydata)))[1:66]

y61=sort(mydata)[61]
y6=sort(mydata)[6]
n=66
require(rjags)
model_string <- "model{

  # Likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mu, 1 / sigma ^ 2)
   }

  # Prior for mu
  mu ~ dnorm(0, 1 /100 ^ 2)
 

  # Prior for sigma
  sigma ~ dexp(1 / 1.5)

}"

model <- jags.model(textConnection(model_string), 
                    data=list(n = n, y = mydata ),
                    n.chains = 5)

update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = c("mu", "sigma"), 
                                 n.iter = 1000, progress.bar = "none")

summary(posterior_sample)

params=as.matrix(posterior_sample)
head(params)


####posterior predictive plot
mu=tail(params,200)[,1]

Ty=abs(y61-mu)-abs(y6-mu)

summary(Ty)
rep=apply(tail(params,200),1,function(x){rnorm(n,x[1],x[2])})

Trep=1:200
for(i in 1:200){
 a=sort(rep[,i])
 Trep[i]=abs(a[61]-mu[i])-abs(a[6]-mu[i])
}

plot(Ty,Trep,xlim=c(-15,15),ylim=c(-15,15))

#pvalue figure 6.4
sum(ifelse(Trep>Ty,1,0))/200



