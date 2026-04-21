rm(list=ls(all=TRUE))

library(fastDummies)
library(tidyverse)
require(rjags)

mydata=read.table("C:\\tuDrive2026\\bayesian\\datasets\\presidential.asc",
                  skip=28,header=1)

mydata=na.omit(mydata)
mydata=mydata[which(mydata$year!=1992),]
mydata$Dvote=mydata$Dvote*100
nrow(mydata)#511 check
Y=mydata$Dvote
X=as.matrix(mydata[,5:24])

n=11
N=511
#head(mydata)
#The south states. 
the_south_states=unique(mydata[which(mydata$r1==1),]$state)
max(mydata$state)
############################generate IDs

mydata$yearID=abs(mydata$year-1988)/4+1

#4 hypothetical regions 
r1=1:10
r2=11:20
r3=21:30
r4=31:50
region=rep(0,511)
region[which(mydata$state %in% r1)]<-"N"
region[which(mydata$state %in% r2)]<-"S"
region[which(mydata$state %in% r3)]<-"E"
region[which(mydata$state %in% r4)]<-"W"

mydata=cbind(mydata,region)

mydata$yearregion=paste0(mydata$year,mydata$region)

unique_id=unique(mydata$yearregion)
yearregionID<-rep(0,511)
for(i in unique_id){
 id=which(i==unique_id)
 yearregionID[ which(mydata$yearregion==i)]<-id
}
max(yearregionID)
mydata=cbind(mydata,yearregionID)

unique_id=unique(mydata$region)
regionID<-rep(0,511)
for(i in unique_id){
 id=which(i==unique_id)
 regionID[ which(mydata$region==i)]<-id
}

mydata=cbind(mydata,regionID)

########################generate ID maps


# Create mapping from year to region code
map_yearregion_to_region <- mydata %>%
  select(yearregionID, regionID) %>%
  distinct() %>%
  arrange(yearregionID) %>%
  column_to_rownames("yearregionID") %>%
  pull(regionID) %>%
  as.integer()
map_yearregion_to_region

toregion=ifelse(map_yearregion_to_region==2,2,1) #S is 2

##################################rJAGS
model_string <- "model{

# Likelihood
  for(i in 1:511){
     y[i] ~ dnorm(mu[i], 1/sigma^2);
     mu[i] <- b[]%*%x[i,] +gamma[yearregionID[i]]+delta[yearID[i]]+delta[yearID[i]];
  }
  for(j in 1:20){          
     b[j]~ dnorm(0, 1/100^2);   
  }
# year-region specific parameters
   for(j in 1:44){          
    gamma[j] ~ dnorm(0,1/tau_gamma[toregion[j]]^2);
   }
# region specific parameters
  for(i in 1:2) {
    tau_gamma[i] ~ dexp(1/100);
  }
# delta parameters
    for(j in 1:11){
     delta[j]~ dnorm(0, 1/tau_delta^2);
    }

   # Priors
      tau_delta ~ dexp(1/100);#need to be half normal
      sigma     ~ dexp(1/100); #need to be half normal

}"

model <- jags.model(textConnection(model_string), 
                    data=list(yearID=mydata$yearID,
                              yearregionID=mydata$yearregionID,
                              toregion=toregion,y = Y, x=X), 
                    n.chains = 5);

update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = c("b","gamma","delta","sigma","tau_gamma","tau_delta" ),
                                 n.iter = 100, progress.bar = "none")

summary(posterior_sample)


###posterior predictive plot
b=posterior_sample[[5]]
#b[1:3,]
dim(b)

# Create dummy variables
df <- data.frame(mydata$yearID,mydata$yearregionID)
df<- dummy_cols(df, select_columns = "mydata.yearID")
df <- dummy_cols(df, select_columns = "mydata.yearregionID")
#head(df)

mydata=cbind(mydata,df)
dim(mydata)
x=cbind(X,mydata[,32:86])
#head(x)
dim(b)
dim(x)

yhat=b[,1:75]%*%t(x)
dim(yhat)
sigma=b[,76]
#generate from each sigma 511 times.
yrep=yhat+t(sapply(sigma,function(x){rnorm(511,0,x)}))

T_y=rep(0,100)
            
for(i in 1:100){
 res=yhat[i,]-mydata$Dvote
 df <- data.frame(res,mydata$year)
  T_y[i]=sqrt(mean(aggregate(df[,1], list(df$mydata.year), FUN=mean)$x^2))
}
summary(T_y)

T_rep=rep(0,100)
            
for(i in 1:100){
 res=yhat[i,]-yrep[i,]
 df <- data.frame(res,mydata$year)
  T_rep[i]=sqrt(mean(aggregate(df[,1], list(df$mydata.year), FUN=mean)$x^2))
}
summary(T_rep)

plot(T_y,T_rep,xlim=c(0,5),ylim=c(0,5))
abline(a=0,b=1,col="red")

#Note: I have nowhere to find the correct region labels, 
#if anyone have more data of the BDA book please contact zli3 %at% live com.     
#https://sites.stat.columbia.edu/gelman/book/data/
