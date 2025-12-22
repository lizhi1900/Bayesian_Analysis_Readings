
#Reference: Bayesian Data Analysis By Andrew Gelman, 2003
#zli3  [at] live.com


rm(list=ls(all=TRUE))
setwd("C:\\tuDrive2025\\bayesian\\presidential\\")
mydata=read.table("C:\\tuDrive2025\\bayesian\\datasets\\presidential.asc",
                  skip=28,header=1)

mydata=na.omit(mydata)
mydata=mydata[which(mydata$year!=1992),]
n=nrow(mydata)#######511 check
head(mydata)

#######year 11 indicators #########
year=mydata$year
head(year)
y=unique(year)
y_ind=sapply(y,function(x){ifelse(year==x,1,0)})
colnames(y_ind)<-paste0("y",y)
head(y_ind)

######region x year, 44 indicators #####
Northeast=c(  7,  8, 19, 20, 21, 29, 30, 32, 38, 39, 45, 48)
South=c(  1,  4,  9, 10, 17, 18, 24, 33, 36, 40, 42, 43, 46)
Midwest=c(  13, 14, 15, 16, 22, 23, 25, 27, 34, 35, 41, 49)
West=c(  2,  3,  5,  6, 11, 12, 26, 28, 31, 37, 44, 47, 50)
ne=ifelse(mydata$state %in% Northeast,1,0)
sth=ifelse(mydata$state %in% South,1,0)
mdw=ifelse(mydata$state %in% Midwest,1,0)
wst=ifelse(mydata$state %in% West,1,0)

ne_y=ne*y_ind
colnames(ne_y)<-paste0("ne",y)
sth_y=sth*y_ind
colnames(sth_y)<-paste0("sth",y)
mdw_y=mdw*y_ind
colnames(mdw_y)<-paste0("mdw",y)
wst_y=wst*y_ind
colnames(wst_y)<-paste0("wst",y)
ry_ind=cbind(ne_y,mdw_y,wst_y,sth_y)
head(ry_ind)

########Table 15.1 variables, exclude regional subregional varibles ######
mydata=mydata[,c(1,5:18)]
mydata$Dvote=mydata$Dvote*100
head(mydata)

###### full set 1+14+44+11=70 #######
mydata=cbind(mydata,ry_ind,y_ind)
ncol(mydata)
head(mydata)






require(rjags)
model_string <- "model{

  # Likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mu[i], 1 / sigma ^ 2)
    mu[i] <- beta[]%*% x[i,]
  }

  # Prior for beta
  for (j in 1:14) { #vague Prior 
    beta[j] ~ dnorm(0, 1 / 100 ^ 2)
  }
  for (j in 15:47) { #non south Prior 
    beta[j] ~ dnorm(0, 1 / tau[1] ^ 2)
  }
  for (j in 48:58) { # south Prior 
    beta[j] ~ dnorm(0, 1 / tau[2] ^ 2)
  }
  for (j in 59:69) { #year Prior 
    beta[j] ~ dnorm(0, 1 / tau[3] ^ 2)
  }

  # Prior for sigma
  sigma ~ dexp(1 / 1.5)

  # Prior for tau
  for (j in 1:3) { 
    tau[j] ~ dexp(1 / 1.5)
  }

}"

model <- jags.model(textConnection(model_string), 
                    data=list(n = n, y = mydata$Dvote, x=mydata[,2:ncol(mydata)]), 
                    n.chains = 5)

update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = c("beta","tau", "sigma"), 
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
plot(T_y,T_rep,xlim=c(0,3),ylim=c(0,3))

#
