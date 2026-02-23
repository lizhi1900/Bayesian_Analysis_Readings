rm(list=ls(all=TRUE))

kidney=read.table("C:\\tuDrive2026\\bayesian\\datasets\\KidneyCancerClean.csv",
                   skip=4,header=TRUE,sep=",")
#zli3 a_t live  com, looking for more data of the book. 

mydata=kidney[,c('state', 'Location', 'fips', 'dc', 'pop', 'dc.2', 'pop.2')]

mydata$dct = mydata$dc + mydata$dc.2 #death count total
mydata$popm = (mydata$pop + mydata$pop.2) / 2
head(mydata)

mydata$naiveproportion = mydata$dct / mydata$popm
summary(mydata$naiveproportion)
hist(mydata$naiveproportion, 300)

threshold = mydata$naiveproportion[order(mydata$naiveproportion, decreasing=TRUE)[100]]
print(threshold)
mydata$cancerhigh = mydata$naiveproportion >= threshold

###########################################################
library(ggplot2)
library(stringr)

state <- map_data("state")
counties <- map_data("county")
counties$Location=paste0(str_to_title(counties$subregion)," County, ",str_to_title(counties$region)) 

naive <- merge(counties,mydata,by="Location")
naive <- subset(naive, select = -c(Location,state))
naive <- naive[order(naive$order), ]
naive <- subset(naive, cancerhigh==TRUE)

#https://jtr13.github.io/cc19/different-ways-of-plotting-u-s-map-in-r.html
ca_map <- ggplot(data=state, mapping=aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black", fill="gray") + 
  geom_polygon(data=naive, fill=NA, color="white") + 
  geom_polygon(color="black", fill=NA) + 
  ggtitle('Map Kidney Cancer Naive Highest 100 Counties') + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ca_map

####################################################

largepops <- subset(mydata,popm>=300000)
cancerhigh <- subset(mydata,cancerhigh==TRUE)

nrow(largepops)
nrow(cancerhigh)

summary(largepops$naiveproportion)
summary(cancerhigh$popm)
summary(largepops$popm)


##################################
m=mean(largepops$naiveproportion)
V=var(largepops$naiveproportion)
c(m,V)
a = ((m*m*(1-m))/V) - m
b = (((1-m)*(1-m)*m)/V) - (1-m)
c(a,b)

df<-data.frame(x= seq(0,8e-4,length=40))
eq<-function(x){dbeta(x,a,b)}
ggplot(data=df,aes(x=x))+stat_function(fun=eq)

################################################
mydata$bayes1=(a + mydata$dct) / (a + b + mydata$popm)
threshold = mydata$bayes1[order(mydata$bayes1, decreasing=TRUE)[100]]
print(threshold)
mydata$bayeshigh1 = mydata$bayes1 >= threshold

table(mydata$cancerhigh,mydata$bayeshigh1) ##check

bayes1 <- merge(counties,mydata,by="Location")
bayes1 <- subset(bayes1, select = -c(Location,state))
bayes1 <- bayes1[order(bayes1$order), ]
bayes1 <- subset(bayes1, bayeshigh1==TRUE)

ca_map <- ggplot(data=state, mapping=aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black", fill="gray") + 
  geom_polygon(data=bayes1, fill=NA, color="white") + 
  geom_polygon(color="black", fill=NA) + 
  ggtitle('Map Kidney Cancer Bayesian Estimation 1. Highest 100 Counties') + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ca_map

hist(mydata$naiveproportion, col='blue',main="Histograms of the Naive and Bayesian estimates", xlab='Estimates',breaks=100)
hist(mydata$bayes1, col='red', add=TRUE,breaks=100)
#add legend
legend('topright', c('Naive Estimate', 'Bayes Estimate'), fill=c('blue', 'red'))

#############################################
#predict the proportions for the second time period. 

mydata$naiveproportion1 = mydata$dc / mydata$pop
mydata$proportion2 = mydata$dc.2 / mydata$pop.2

proportions1_largecounties = mydata$naiveproportion1[mydata$pop >= 300000]
m1 = mean(proportions1_largecounties)
V1 = var(proportions1_largecounties)
a1 = ((m1*m1*(1-m1))/V1) - m1
b1 = (((1-m1)*(1-m1)*m1)/V1) - (1-m1)

mydata$bayesestimate1 = (mydata$dc + a1)/(mydata$pop + a1 + b1)

diff_naiveproportion1 = mydata$naiveproportion1 - mydata$proportion2
sum_abs_diff_naiveproportion1 = sum(abs(diff_naiveproportion1))
sum_squared_diff_naiveproportion1 = sum((diff_naiveproportion1**2))
kl_loss_naiveproportion1 = sum(mydata$proportion2 * log(mydata$proportion2 / (mydata$naiveproportion1 + 1e-10) + 1e-10)) + sum((1 - mydata$proportion2) * log((1 - mydata$proportion2) / (1 - mydata$naiveproportion1 + 1e-10) + 1e-10))
negloglik_naiveproportion1 = -sum(mydata$dc.2 * log(mydata$naiveproportion1 + 1e-10)) - sum((mydata$pop.2 - mydata$dc.2) * log(1 - mydata$naiveproportion1 + 1e-10))

diff_bayesestimate1 = mydata$bayesestimate1 - mydata$proportion2
sum_abs_diff_bayesestimate1 = sum(abs(diff_bayesestimate1))
sum_squared_diff_bayesestimate1 = sum(diff_bayesestimate1**2)
kl_loss_bayesestimate1 = sum(mydata$proportion2 * log(mydata$proportion2 / (mydata$bayesestimate1 + 1e-10) + 1e-10)) + sum((1 - mydata$proportion2) * log((1 - mydata$proportion2) / (1 - mydata$bayesestimate1 + 1e-10) + 1e-10))
negloglik_bayesestimate1 = -sum(mydata$dc.2 * log(mydata$bayesestimate1 + 1e-10)) - sum((mydata$pop.2 - mydata$dc.2) * log(1 - mydata$bayesestimate1 + 1e-10)) 

bothresults_absdiff = c(sum_abs_diff_naiveproportion1, sum_abs_diff_bayesestimate1)
bothresults_squareddiff = c(sum_squared_diff_naiveproportion1, sum_squared_diff_bayesestimate1)
bothresults_klloss = c(kl_loss_naiveproportion1, kl_loss_bayesestimate1)
bothresults_negloglik = c(negloglik_naiveproportion1, negloglik_bayesestimate1)

print(bothresults_absdiff)
print(bothresults_squareddiff)
print(bothresults_klloss)
print(bothresults_negloglik)

###########################################################################
#Empirical Bayes

require(VGAM)

loglik<-function(lambda){
 terms=dbetabinom.ab(mydata$dct,mydata$popm, lambda[1],lambda[2],log=TRUE)
 return(sum(terms))
}

out <- optim(par = c(8, 100000), fn = loglik, control = list(fnscale = -1),
            hessian = TRUE)
cat('Posterior mode:\n')# Empirical Bayes Estimates: 15.70324 151539.38273
print((post_mode <- out$par))

#Reference
#https://github.com/Alexchen666/kidney-cancer/blob/main/KidneyCancerClean.csv
#https://stat238.berkeley.edu/spring-2026/codelecturenine238spring2026/
https://sites.stat.columbia.edu/gelman/book/data/