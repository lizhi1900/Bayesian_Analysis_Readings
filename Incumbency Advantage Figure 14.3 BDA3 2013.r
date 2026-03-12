#Figure 14.3
rm(list=ls(all=TRUE))
library(ggplot2)
#https://sites.stat.columbia.edu/gelman/book/data/
#zli3 a_t live  com, looking for more data of the book. 

demvote1980s=NULL
res=NULL
incumbent=NULL
year=c(1980,1982,1984,1986,1988)

for(i in year){
 this=read.table(paste0("C:\\tuDrive2026\\bayesian\\datafolder\\",i,".asc"))
 last=read.table(paste0("C:\\tuDrive2026\\bayesian\\datafolder\\",(i-2),".asc"))

 last$id=last$V1+last$V2*.01 #id
 this$id=this$V1+this$V2*.01
 last=last[apply(last!=-9, 1, all),] #Drop -9
 this=this[apply(this!=-9, 1, all),]
 this=this[which(!(this$V4==0|this$V5==0)),] #Drop 0 votes
 last=last[which(!(last$V4==0|last$V5==0)),] 
 mydata=merge(last,this,by="id")

#incumbent party, winner in past, democrat 1, repub -1
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
 demvote=mydata$V4.x/(mydata$V4.x+mydata$V5.x)

 model <- lm(y ~ R + Incumbentpastvote + Incumbentparty)

 res <- c(res, rstandard(model))
 demvote1980s=c(demvote1980s,demvote)
 incumbent <- c(incumbent, ifelse(R==0,"Open","Not Open"))
}

########################################################
df_plot <- data.frame(
  fitted = demvote1980s,
  std_resid = res,
  label = incumbent
)
dim(df_plot)
ggplot(df_plot, aes(x = fitted, y = std_resid)) +
  geom_point(color='black', shape=21, size=1, aes(fill=factor(label ))) + 
  scale_fill_manual(values=c('black','white' ))+
  labs(
    x = "Democratic vote in the previous election",
    y = "Standardized Residuals",
    title = "Figure 14.3"
  )
