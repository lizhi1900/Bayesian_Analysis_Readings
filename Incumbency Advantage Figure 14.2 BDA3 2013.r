#Figure 14.2

rm(list=ls(all=TRUE))
library(LearnBayes)
years=seq(1900,1992,by=2)
bars=NULL

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
 R=abs(mydata$V3.y)
 y=ifelse(winner==1,
          mydata$V4.y/(mydata$V4.y+mydata$V5.y),
          mydata$V5.y/(mydata$V4.y+mydata$V5.y))

 Incumbentparty=winner 
 Incumbentpastvote=winnervote


 fit <- lm(y ~ R + Incumbentpastvote + Incumbentparty, x=TRUE, y=TRUE )
 theta.sample <- blinreg(fit$y, fit$x, 5000)
 out <- quantile(theta.sample$beta[,2], c(.05, .5, .95))
 bars <-cbind(bars,out)
}

#posterior median and 95% interval
matplot(rbind(years, years), bars[c(1,3),], 
        type="l", lty=1, col=1,
        xlab="Year", ylab="Estimated incumbency advantage")
points(years, bars[2,], pch=19)


