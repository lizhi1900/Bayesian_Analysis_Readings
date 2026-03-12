#table 14.1
rm(list=ls(all=TRUE))
library(LearnBayes)
years=1988
#https://sites.stat.columbia.edu/gelman/book/data/
#zli3 a_t live  com, looking for more data of the book. 
for(i in years){
 this=read.table(paste0("C:\\tuDrive2026\\bayesian\\datafolder\\",i,".asc"))
 last=read.table(paste0("C:\\tuDrive2026\\bayesian\\datafolder\\",i-2,".asc"))

 last$id=last$V1+last$V2*.01 #id
 this$id=this$V1+this$V2*.01
 last=last[apply(last!=-9, 1, all),] #Drop -9
 this=this[apply(this!=-9, 1, all),]
 this=this[which(!(this$V4==0|this$V5==0)),] #Drop 0 votes
 last=last[which(!(last$V4==0|last$V5==0)),] 

 mydata=merge(last,this,by="id")

#incumbent party, winner in last, democrat 1, repub -1
 winner=ifelse(mydata$V4.x>mydata$V5.x,1,-1)
 winnervote=ifelse(winner==1,
                   mydata$V4.x/(mydata$V4.x+mydata$V5.x),
                   mydata$V5.x/(mydata$V4.x+mydata$V5.x))
 R=ifelse(mydata$V3.y==0,0,1)
 y=ifelse(winner==1,
          mydata$V4.y/(mydata$V4.y+mydata$V5.y),
          mydata$V5.y/(mydata$V4.y+mydata$V5.y))

 Incumbentparty=winner 
 Incumbentpastvote=winnervote

 fit <- lm(y ~ R + Incumbentpastvote + Incumbentparty, x=TRUE, y=TRUE )
 theta.sample <- blinreg(fit$y, fit$x, 5000)
}
#posterior distributions
par(mfrow=c(2,2))
hist(theta.sample$beta[,2], main="Incumbency",
  xlab=expression(beta[1]))
hist(theta.sample$beta[,3], main="Vote porportion in 1986",
  xlab=expression(beta[2]))
hist(theta.sample$beta[,4], main="Incumbent party",
  xlab=expression(beta[3]))
hist(theta.sample$sigma, main="ERROR SD",
  xlab=expression(sigma))

#Table 14.1 posterior quantiles 
apply(theta.sample$beta, 2, quantile, c(.05, .5, .95))
quantile(theta.sample$sigma, c(.05, .5, .95))

