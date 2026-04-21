rm(list=ls(all=TRUE))
require(LearnBayes)

mydata=read.table("C:\\tuDrive2026\\bayesian\\datasets\\presidential.asc",
                  skip=28,header=1)
mydata=na.omit(mydata)
mydata=mydata[which(mydata$year!=1992),]
nrow(mydata)#511 check
mydata$Dvote=mydata$Dvote*100
#head(mydata,20)
Y=mydata$Dvote
X=as.matrix(mydata[,5:24])
year=mydata$year

####################################################

n=200
theta.sample <- blinreg(Y ,X ,n)
y.hat=blinregexpected(X,theta.sample)
y.rep = blinregpred(X,theta.sample)
y0=theta.sample$beta%*%t(X)
#head(y0-y.hat)
dim(y.hat)
dim(y.rep)

#####################################################
###posterior predictive scatter plot

T_y=rep(0,n)
              
for(i in 1:n){
 res=(y.hat[i,]-Y)
 df <- data.frame(res,year)
 T_y[i]=sqrt(mean(aggregate(df[,1], list(df$year), FUN=mean)$x^2))
}

T_rep=rep(0,n)

for(i in 1:n){
 res=(y.hat[i,]-y.rep[i,])
 df <- data.frame(res,year)
 T_rep[i]=sqrt(mean(aggregate(df[,1], list(df$year), FUN=mean)$x^2))
}

plot(T_y,T_rep,xlim=c(0,1.5),ylim=c(0,1.5))
abline(a=0,b=1,col="red")

#https://sites.stat.columbia.edu/gelman/book/data/