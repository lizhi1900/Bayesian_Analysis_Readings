#Table14.2

rm(list=ls(all=TRUE))
library(LearnBayes)
#https://sites.stat.columbia.edu/gelman/book/data/
#zli3 a_t live  com, looking for more data of the book. 

years=seq(1900,1992,by=2)
residuals=NULL 
incumbent_running=NULL 
n=1000
R0outliertotal=rep(0,n)
R1outliertotal=rep(0,n)
R0total=rep(0,n)
R1total=rep(0,n)

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
 print(i)

 #####################################################
 model <- lm(y ~ R+Incumbentpastvote+Incumbentparty, x=TRUE, y=TRUE )
 residuals=c(residuals, model$residuals)
 incumbent_running=c(incumbent_running,R)

 ##############################################
 #Bayesian predictive replications
 theta.sample <- blinreg(model$y, model$x, n)
 y.rep = blinregpred(model$x,theta.sample)
 model.rep.residuals <-apply(y.rep,1,function(x){lm(x ~ model$x)$residuals})
 outlier.indcator <-   model.rep.residuals>.2
 
 R0= outlier.indcator[which(R==0),]
 R0outliertotal=R0outliertotal+colSums(R0)
 R0total=R0total+rep(nrow(R0),n)

 R1=outlier.indcator[which(R==1),]
 R1outliertotal=R1outliertotal+colSums(R1)
 R1total=R1total+rep(nrow(R1),n)
}


##########################################################3
#Table 14.2 Observed proportion of outliers
library(gmodels)

#perform cross-tabulation of incumbent_running and outliers preference
CrossTable(x=incumbent_running, y=abs(residuals)>.2, digits = 5,
           prop.c=FALSE,prop.t=FALSE,prop.chisq=FALSE)

#####################################################################
#Table 14.2 Posterior predictive distribution of proportion of outliers
propo.open.outliers=R0outliertotal/R0total
propo.incumbent.running.outliers=R1outliertotal/R1total

quantile(propo.open.outliers,names = TRUE,c(.025,.5,.975),digits = 5)
quantile(propo.incumbent.running.outliers,names = TRUE,c(.025,.5,.975),digits = 5)




