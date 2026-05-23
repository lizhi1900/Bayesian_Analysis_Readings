rm(list=ls(all=TRUE))
require(FNN)
require(rjags)
library(dplyr)
library(ggplot2)

#https://www.itl.nist.gov/div898/handbook/pmd/section1/dep/dep144.htm
#http://127.0.0.1:25624/library/stats/html/lowess.html
#Weighted linear regression, BDA13 page 372
#https://www.statology.org/weighted-least-squares-in-r/
#zli3 ~a t~ live com

train=read.table("https://www.itl.nist.gov/div898/handbook/datasets/LOESSCMP.DAT",skip=25)
#train=read.table("C:\\tuDrive2026\\bayesian\\loess\\LOESSCMP.DAT",skip=25)

x=train[,1]
y=train[,2]

xi=seq(min(x),max(x),by=.1)#for prediction
n=length(xi)

ID=rep(1:n,each=7)#k=7
indices<-as.vector(t(knnx.index(x, xi, k=7)))
mydata=cbind(ID,x=rep(xi,each=7),train[indices,])
mydata[1:7,]
N=nrow(mydata)

wt =NULL
for(i in 1:n){
  d=mydata[which(mydata$ID==i),]
  model<-lm(d$V2~d$V1)
  tem=1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
 wt <- c(wt,tem)
}

#######################################################
model_string <- "model{

  # Likelihood
  for(i in 1:N){
     y[i] ~ dnorm(a[ID[i]] + b[ID[i]]*x[i], w[i]/sigma[ID[i]]^2);
     
  }

   # local regression specific parameters
   for(j in 1:n){
      a[j] ~ dnorm(mu_a, 1/sigma_a^2);
      b[j]  ~ dnorm(mu_b, 1/sigma_b^2);
      sigma[j] ~ dexp(1/100);
   }
    
   # Priors
      mu_a    ~ dnorm(0, 1/500^2);
      sigma_a ~ dnorm(0,1/100^2)T(0,);
      mu_b     ~ dnorm(0, 1/500^2);
      sigma_b  ~ dnorm(0,1/100^2)T(0,);
   #prediction
    for(j in 1:n){
       y_pred[j]~dnorm(mu[j],1/sigma[j]^2)
       mu[j]=a[j]+b[j]*xi[j]
    }

}"

data_list<-list(N=N,n = n,ID=ID, 
                y = mydata$V2, x=mydata$V1, 
                w = wt,xi=xi) 
 
model <- jags.model(textConnection(model_string), data=data_list, n.chains = 5)
    
update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                    variable.names = c( "mu","y_pred"),                   
                    n.iter = 100, 
                    progress.bar = "none")


summary(posterior_sample)

weight_chains <- data.frame(posterior_sample[[1]])
y_pred=weight_chains|> select(matches("y_pred"))
y_mu=weight_chains|> select(matches("mu"))

interval95=apply(y_pred,2,function(x) quantile(x,c(.025,.975)))

###################################################
#plot predicted mean and 95% confidence interval

df<-data.frame(xi,colMeans(y_mu),t(interval95))
 colnames(df)<-c("x","pred_mean","low","high")


ggplot(df, aes(x, pred_mean)) +
    geom_point(train, mapping =aes(x = x, y = y),alpha=0.5,color ="red") +
    geom_point(aes(x = x, y = pred_mean),alpha=0.5,color ="blue") +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "blue") +
    #geom_abline(intercept = mean(unlist(posterior_sample[,1,])), slope = mean(unlist(posterior_sample[,177,])), color = "orange", size = 1.2) +
    #geom_smooth(aes(x = x, y = pred_mean),method = "loess", color = "red", se = FALSE)  # Add regression line without confidence interval
      geom_line(aes(x = x, y = pred_mean))
 









