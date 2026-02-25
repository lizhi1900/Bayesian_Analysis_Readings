rm(list=ls(all=TRUE))

library(haven)
library(tidyverse)
require(VGAM)

data <- read_dta('C:\\tuDrive2026\\bayesian\\datasets\\pew_research_center_june_elect_wknd_data.dta')
obama<- read.csv("C:\\tuDrive2026\\bayesian\\datasets\\2008ElectionResult.csv",header=TRUE)
#zli3 a_t live  com, looking for more data of the book. 


#Count very very liberal in states
state_count<-data %>% filter(ideo == 5)%>% group_by(state) %>%  count()
  print(state_count,n=49)

#Missing is not a respond, exclude 15 hawaii
state_respondent<-data %>% filter(ideo != 0 & state != 15) %>% 
  group_by(state) %>% count
  print(state_respondent,n=49)

propotion = state_count / state_respondent
obamapct=obama$vote_Obama_pct[which(obama$state!="Alaska" & obama$state!="Hawaii")]/100

###########################################################################
#Empirical Bayes

require(VGAM)

loglik<-function(lambda){
 terms=dbetabinom.ab(state_count$n,state_respondent$n, lambda[1],lambda[2],log=TRUE)
 return(sum(terms))
}

out <- optim(par = c(8, 40), fn = loglik, control = list(fnscale = -1),
            hessian = TRUE)
cat('Posterior mode:\n')
print((post_mode <- out$par))


############################################

#hyper prior p(log(a),log(b))=1, jacobian is 1/a/b
as.numeric(state_count$n)

logpost<-function(lambda){
 b=exp(lambda[2])
 a=exp(lambda[1])
 terms=dbetabinom.ab(state_count$n,state_respondent$n, a,b,log=TRUE)
 return(sum(terms)-log(a)-log(b))
}

out <- optim(par = c(2, 8), fn = logpost, control = list(fnscale = -1),
            hessian = TRUE)
cat('Posterior mode:\n')
print((post_mode <- out$par))
cat('\n')
cat('Covariance matrix: \n')
print((post_cov <- -solve(out$hessian)))

###################### ## ############################### 
lambda_1 <-seq(-2.5,8,by=.1)
lambda_2 <-seq(1.5 , 10,by=.1)
z <- matrix(0, nrow = length(lambda_1), ncol = length(lambda_2))

for (i in 1:length(lambda_1)){
  for (j in 1:length(lambda_2)){
    lambda <- c(lambda_1[i], lambda_2[j])
    z[i,j] <- logpost(lambda)
  }
}


contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)

##############################333
library(LearnBayes)
library(coda)

set.seed(11)

iters <- 10^4
proposal <- list(var = post_cov, scale = 2)

# random walk metropolis
fit1 <- rwmetrop(logpost, proposal, start = post_mode, iters )

# overlaying last 5000 draws on contour plot of the log posterior
contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)
points(x = fit1$par[5001:iters,1], y = fit1$par[5001:iters,2], col = "red")

###################################################
a=exp(fit1$par[5001:iters,1])
b=exp(fit1$par[5001:iters,2])
theta=1:49
for(j in 1:49){
theta[j] = mean((a+state_count$n[j])/(a+b+state_respondent$n[j]))
}

df <- data.frame(theta,propotion$n,state_respondent$n,obamapct )

ggplot(df, aes(x = state_respondent.n)) + 
  geom_line(aes(y = theta, color = 'bayes mean')) + 
  geom_line(aes(y = obamapct, color = 'vote share'))+
  geom_line(aes(y = propotion.n, color = 'very liberal'))

Reference
#https://stat238.berkeley.edu/spring-2026/codelectureten238spring2026/
#https://sites.stat.columbia.edu/gelman/book/data/
#https://jytan.net/blog/2023/eight-schools/

