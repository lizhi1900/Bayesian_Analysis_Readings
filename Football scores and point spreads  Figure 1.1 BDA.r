rm(list=ls(all=TRUE))
football=read.table("https://sites.stat.columbia.edu/gelman/book/data/football.asc",
                   skip=8,header=FALSE)
#zli3 a_t live  com, looking for more data of the book. 


names(football) <- c("home", "favorite", 
"underdog", "spread", "favorite.name",
"underdog.name", "week")

head(football)
nrow(football)


##################################################################

mydata=football[1:672,]

outcome=mydata$favorite-mydata$underdog
#figure 1.1
plot(mydata$spread,outcome)

d=outcome-mydata$spread
#figure 1.2 a
plot(mydata$spread,d)
#figure 1.2 b  
hist(d,prob=TRUE,breaks=40,main="outcome - point spread")

x_vals<-seq(min(d),max(d),length=40)
y_vals<-dnorm(x_vals,mean(d),sd(d))
lines(x_vals,y_vals,col="red",lwd=2)

