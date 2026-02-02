rm(list=ls(all=TRUE))
rats=read.csv("https://sites.stat.columbia.edu/gelman/book/data/rats.asc",sep=" ", skip=3,header=1)
#zli3 @ live  com

require(VGAM)

post<-function(lambda,mydata){
 b=exp(lambda[2])/(exp(lambda[1])+1)
 a=b*exp(lambda[1])
 terms=dbetabinom.ab(mydata[,1],mydata[,2], a,b)
 return(a*b*(a+b)^(-5/2)*prod(terms))
}

logpost<-function(lambda,mydata){
 b=exp(lambda[2])/(exp(lambda[1])+1)
 a=b*exp(lambda[1])
 terms=dbetabinom.ab(mydata[,1],mydata[,2], a,b,log=TRUE)
 return(log(a*b*(a+b)^(-5/2))+sum(terms))
}

###################### 5.2 ############################### 
lambda_1 <-seq(-2.5,-1,by=.1)
lambda_2 <-seq(1.5 , 3,by=.1)
z <- matrix(0, nrow = length(lambda_1), ncol = length(lambda_2))

for (i in 1:length(lambda_1)){
  for (j in 1:length(lambda_2)){
    lambda <- c(lambda_1[i], lambda_2[j])
    z[i,j] <- post(lambda, rats)
  }
}


contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)

##################### 5.3a ############################
lambda_1 <-seq(-2.3,-1.3,by=.1)
lambda_2 <-seq(1 , 5,by=.1)
z <- matrix(0, nrow = length(lambda_1), ncol = length(lambda_2))

for (i in 1:length(lambda_1)){
  for (j in 1:length(lambda_2)){
    lambda <- c(lambda_1[i], lambda_2[j])
    z[i,j] <- post(lambda, rats)
  }
}



contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)


##################################### 5.3b ############################

out <- optim(par = c(8, 2), fn = logpost, control = list(fnscale = -1),
            hessian = TRUE, mydata=rats)
cat('Posterior mode:\n')
print((post_mode <- out$par))
cat('\n')
cat('Covariance matrix: \n')
print((post_cov <- -solve(out$hessian)))

library(LearnBayes)
library(coda)

set.seed(11)

iters <- 10^4
proposal <- list(var = post_cov, scale = 2)

# random walk metropolis
fit1 <- rwmetrop(logpost, proposal, start = post_mode, iters, rats)

# overlaying last 5000 draws on contour plot of the log posterior
contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)
points(x = fit1$par[5001:iters,1], y = fit1$par[5001:iters,2], col = "red")

#References
#Learning Bayesian Hierarchical Modeling from 8 Schools, Jun Yu Tan 2023
#Bayesian Statistics class, ST4234 in NUS, taught by Prof Li Cheng
#Bayesian Data Analysis, Third Edition A. Gelman, J.B. Carlin, H.S. Stern, 
#and 3 more authors 2013  