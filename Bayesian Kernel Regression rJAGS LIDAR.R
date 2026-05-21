rm(list=ls(all=TRUE))
#zli3 $at$ live com
require(rjags)
library(dplyr)
library(ggplot2)

mydata=read.table("http://www.stat.cmu.edu/~larry/all-of-nonpar/=data/lidar.dat",header=TRUE)
mydata[1:6,]
x=seq(min(mydata$range),max(mydata$range),by=5)
N <- nrow(mydata)
ggplot(mydata, aes(range, logratio)) +
    geom_point(aes(x = range, y = logratio),alpha=0.5,color ="red")

######################################################
model_code <- "
model {
  #likelihood
  for (i in 1:N) {
    for (j in 1:N) {
     kernel[i,j]=exp(-(X[i]-X[j])^2/2/a^2)
    }
    Y[i]~ dnorm(mu[i],sigma^-2)
    mu[i]=kernel[i,]%*%w/sum(kernel[i,])
  }
 #prediction
    for (i in 1:length(x)) {
     for (j in 1:N) {
      k[i,j]=exp(-(x[i]-X[j])^2/2/a^2)
     }
     y_pred[i]~ dnorm(m[i],sigma^-2)
     m[i]=k[i,]%*%w/sum(k[i,]) 
    }  

 #priors
 for(j in 1:N){
    w[j] ~ dnorm(0,1)
 }
 a~dnorm(0,1)T(0,) #bandwidth 
 sigma~dnorm(0,1)T(0,)
}
"
data_list <- list(
  N = N,
  Y=mydata$logratio,
  X=mydata$range,
  x=x
  
)
model <- jags.model(textConnection(model_code), 
                    data=data_list, 
                    n.chains = 5);
update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = c("a","mu","sigma","w","y_pred"),
                                 n.iter = 100, progress.bar = "none")

#summary(posterior_sample)

####################################################

weight_chains <- data.frame(posterior_sample[[1]])

w= weight_chains|> select(matches("w"))
a=weight_chains|> select(matches("^a$"))
sigma=weight_chains|> select(matches("sigma"))
y_pred=weight_chains|> select(matches("y_pred"))

interval95=apply(y_pred,2,function(x) quantile(x,c(.025,.975)))

###################################################
#plot predicted mean and 95% confidence interval

df<-data.frame(x,colMeans(y_pred),t(interval95))
 colnames(df)<-c("x","pred_mean","low","high")


ggplot(df, aes(x, pred_mean)) +
    geom_point(mydata, mapping =aes(x = range, y = logratio),alpha=0.5,color ="red") +
    geom_point(aes(x = x, y = pred_mean),alpha=0.5,color ="blue") +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "blue") +
    #geom_abline(intercept = mean(unlist(posterior_sample[,1,])), slope = mean(unlist(posterior_sample[,177,])), color = "orange", size = 1.2) +
    #geom_smooth(aes(x = x, y = pred_mean),method = "loess", color = "red", se = FALSE)  # Add regression line without confidence interval
      geom_line(aes(x = x, y = pred_mean))
 





