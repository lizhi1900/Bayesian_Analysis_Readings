rm(list=ls(all=TRUE))
# zli3 at  live com

# given data
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

# defining the log posterior for lambda
logpost <- function(lambda, sigma, y){
  a1=lambda[1]
  a2=exp(lambda[2])^2+sigma^2
  sum(dnorm(y,a1,sqrt(a2),log =1))+lambda[2]
}

# grids
lambda_1 <- seq(from = -18, to = 37, by = 0.1)
lambda_2 <- seq(from = -6, to = 4.1, by = 0.1)
z <- matrix(0, nrow = length(lambda_1), ncol = length(lambda_2))

for (i in 1:length(lambda_1)){
  for (j in 1:length(lambda_2)){
    lambda <- c(lambda_1[i], lambda_2[j])
    z[i,j] <- logpost(lambda, sigma, y)
  }
}

contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)


out <- optim(par = c(8, 2), fn = logpost, control = list(fnscale = -1),
            hessian = TRUE, sigma = sigma, y = y)
cat('Posterior mode:\n')
print((post_mode <- out$par))#7.926685 1.841525

cat('\n')
cat('Covariance matrix: \n')
print((post_cov <- -solve(out$hessian)))


library(LearnBayes)
library(coda)

set.seed(11)

iters <- 10^4
proposal <- list(var = post_cov, scale = 2)

# random walk metropolis
fit1 <- rwmetrop(logpost, proposal, start = post_mode, iters, sigma, y)

# overlaying last 5000 draws on contour plot of the log posterior
contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)
points(x = fit1$par[5001:iters,1], y = fit1$par[5001:iters,2], col = "red")


cat('Acceptance rate: \n')
print(fit1$accept)

par(mfrow=c(2,1))
plot(density(fit1$par[5001:iters,1]), main = "", xlab = expression(lambda[1]))
plot(density(fit1$par[5001:iters,2]), main = "", xlab = expression(lambda[2]))

#######################################################################

#thetaj|mu,tau,y,sigma ~ N(thetajhat, Vhat)
tau=1:40
mu=7.926685

tm=sapply(tau,function(x)(y/sigma^2+mu/x^2)/(1/sigma^2+1/x^2))
V=sapply(tau,function(x)(1/(1/sigma^2+1/x^2)))

matplot(t(tm), type = "b",pch=1,col = 1:8)
legend("topright", legend = 1:8, col=1:8, pch=1)

matplot(t(sqrt(V)), type = "b",pch=1,col = 1:8)
legend("topright", legend = 1:8, col=1:8, pch=1)

#References
#Learning Bayesian Hierarchical Modeling from 8 Schools, Jun Yu Tan 2023
#Bayesian Statistics class, ST4234 in NUS, taught by Prof Li Cheng
#Bayesian Data Analysis, Third Edition A. Gelman, J.B. Carlin, H.S. Stern, 
#and 3 more authors 2013  

