#Table14.2

rm(list=ls(all=TRUE))
############################################
library(ggplot2)
years=1988

for(i in years){
 y1988=read.table(paste0("C:\\tuDrive2026\\bayesian\\datafolder\\",i,".asc"))
 y1986=read.table(paste0("C:\\tuDrive2026\\bayesian\\datafolder\\",i-2,".asc"))

 y1986$id=y1986$V1+y1986$V2*.01 #id
 y1988$id=y1988$V1+y1988$V2*.01
 y1986=y1986[apply(y1986!=-9, 1, all),] #Drop -9
 y1988=y1988[apply(y1988!=-9, 1, all),]
 y1988=y1988[which(!(y1988$V4==0|y1988$V5==0)),] #Drop 0 votes
 y1986=y1986[which(!(y1986$V4==0|y1986$V5==0)),] 

 mydata=merge(y1986,y1988,by="id")
}
head(mydata)
demvote1986=(mydata$V4.x/(mydata$V4.x+mydata$V5.x))
demvote1988=(mydata$V4.y/(mydata$V4.y+mydata$V5.y))
R=abs(mydata$V3.y)

df_plot <- data.frame(demvote1986,demvote1988,R)

ggplot(df_plot, aes(x = demvote1986, y = demvote1988)) +
  geom_point(color='black', shape=21, size=1, aes(fill=factor(R ))) + 
  scale_fill_manual(values=c('white', 'black'))+
  labs(
    x = "Democratic vote in 1986",
    y = "Democratic vote in 1988",
    title = "Figure 14.1"
  )

