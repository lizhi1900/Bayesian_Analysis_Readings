
#Reference: Bayesian Data Analysis By Andrew Gelman, 2003
#zli3  at   live.com
#Dec/22/2025
 
rm(list=ls(all=TRUE))
setwd("C:\\tuDrive2025\\bayesian\\presidential\\")
mydata=read.table("C:\\tuDrive2025\\bayesian\\datasets\\presidential.asc",
                  skip=28,header=1)

mydata=na.omit(mydata)
mydata=mydata[which(mydata$year!=1992),]
n=nrow(mydata)#511 check
head(mydata)
year=mydata$year

mydata=mydata[,c(1,5:24)]
mydata$Dvote=mydata$Dvote*100

ordinary=lm(Dvote ~.-1 ,data=mydata)
summary(ordinary)

names(ordinary)

ordinary$residuals

df1 <- data.frame(ordinary$residuals,year)
sqrt(mean((aggregate(df1, list(df1$year), FUN=mean)$ordinary.residuals)^2))

M=ncol(mydata)-1

require(rjags)
model_string <- "model{

  # Likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mu[i], 1 / sigma ^ 2)
    mu[i] <- beta[]%*% x[i,]
  }

  # Prior for beta
  for (j in 1:M) {
    beta[j] ~ dnorm(2, 1 / 2 ^ 2)
  }

  # Prior for sigma
  sigma ~ dexp(1 / 1.5)

}"

model <- jags.model(textConnection(model_string), 
                    data=list(n = n, y = mydata$Dvote, x=mydata[,2:ncol(mydata)],M=M), 
                    n.chains = 5)

update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = c("beta", "sigma"), 
                                 n.iter = 100, progress.bar = "none")

summary(posterior_sample)

params = as.matrix(posterior_sample)
nrow(params )

Beta=tail(params,200)[,1:(ncol(mydata)-1)]
Sigma=tail(params,200)[,ncol(mydata)]
X=as.matrix(mydata[,2:ncol(mydata)])


###posterior predictive plot

y_hat=X%*%t(Beta)
dim(y_hat)
res=y_hat-mydata$Dvote

df <- data.frame(res,year)
T_y=rep(0,200)
               
for(i in 1:200){
  T_y[i]=sqrt(mean((aggregate(df[,i], list(df$year), FUN=mean)$x)^2))
}



#T_rep
T_rep=rep(0,200)
               
for(i in 1:200){
  musigma=cbind(X%*%Beta[i,],Sigma[i])
  y_rep=apply(musigma,1,function(x){rnorm(1,x[1],x[2])})
  res=y_rep-X%*%Beta[i,]
  df <- data.frame(res,year)
  
  T_rep[i]=sqrt(mean((aggregate(df$res, list(df$year), FUN=mean)$x)^2))
}
summary(T_y)
summary(T_rep)
plot(T_y,T_rep,xlim=c(0,1.5),ylim=c(0,1.5))

#
