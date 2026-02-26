library(BayesDA)

data(meta)
names(meta)
# Calculating empirical log-odds and its sampling variances:
y <- apply(meta, 1, function(x) log( (x[4]/(x[5]-x[4]))/(x[2]/(x[3]-x[2])) ) )
s2 <- apply(meta, 1, function(x) 1/(x[5]-x[4]) + 1/x[4] +1/(x[3]-x[2]) + 1/x[2] )

sigma<-sqrt(s2)  

logpost <- function(lambda, sigma, y){
  sum(-0.5*log(exp(2*lambda[2])+sigma^2) - 
        ((lambda[1]-y)^2)/(2*(sigma^2+exp(2*lambda[2])))-log(sqrt(pi*2))) +
        lambda[2]
}

# defining the log posterior for lambda
logpost2 <- function(lambda, sigma, y){
  a1=lambda[1]
  a2=exp(lambda[2])^2+sigma^2
  sum(dnorm(y,a1,sqrt(a2),log =1))+lambda[2]
}
t=c(-.2,-4)
logpost(t,sigma,y)#check
logpost2(t,sigma,y) 

#####################################################3
# grids
lambda_1 <- seq(from = min(y), to = max(y), by = 0.1)
lambda_2 <- seq(from = -6, to = 4.1, by = 0.1)
z <- matrix(0, nrow = length(lambda_1), ncol = length(lambda_2))

for (i in 1:length(lambda_1)){
  for (j in 1:length(lambda_2)){
    lambda <- c(lambda_1[i], lambda_2[j])
    z[i,j] <- logpost2(lambda, sigma, y)
  }
}

contour(x = lambda_1, y = lambda_2, z = z, col = "blue", nlevels = 40,
        xlab = expression(lambda[1]), ylab = expression(lambda[2]),
        cex.axis = 1.1, cex.lab = 1.3)

out <- optim(par = c(-.2, -2), fn = logpost2, control = list(fnscale = -1),
            hessian = TRUE, sigma = sigma, y = y)
cat('Posterior mode:\n')
print((post_mode <- out$par))
cat('\n')
cat('Covariance matrix: \n')
print((post_cov <- -solve(out$hessian)))

#######################################################
library(LearnBayes)
library(coda)

set.seed(11)

iters <- 10^4
proposal <- list(var = post_cov, scale = 2)

# random walk metropolis
fit1 <- rwmetrop(logpost2, proposal, start = post_mode, iters, sigma, y)

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

#######################################################
mcmcobj1 <- mcmc(fit1$par[5001:iters,])
colnames(mcmcobj1) <- c("lambda_1", "lambda_2")
par(mfrow=c(2,1))
traceplot(mcmcobj1)

par(mfrow=c(2,1))
autocorr.plot(mcmcobj1, auto.layout = FALSE)

# the last 5000 MCMC samples (lambda_1, lambda_2)
lambda_samples <- fit1$par[5001:iters,]

# function to compute mean
theta_hat <- function(lambda, y_j, sigma_j){
    ((y_j/sigma_j^2)+(lambda[,1]/exp(2*lambda[,2]))) /
    ((1/sigma_j^2)+(1/exp(2*lambda[,2])))
}

# function to compute variance
V <- function(lambda, y_j, sigma_j){
    1 / (1/sigma_j^2 + 1/exp(2*lambda[,2]))
}

# drawing 5000 samples of theta_j
theta_samples <- function(lambda, y_j, sigma_j){
    rnorm(5000, mean = theta_hat(lambda, y_j, sigma_j),
          sd = sqrt(V(lambda, y_j, sigma_j)))
}

theta_mean <- rep(0, 22)
theta_sd <- rep(0,22)

# the joint posterior density of (theta_1,...,theta_j)
theta_all <- matrix(0, nrow = 5000, 22)
for (j in 1:22){
        thetas <- theta_samples(lambda_samples, y[j], sigma[j])
        theta_all[,j] <- thetas
        theta_mean[j] <- mean(thetas)
        theta_sd[j] <- sd(thetas)
}

print(theta_dist <- cbind(theta_mean, theta_sd))

#table 5.5
quantile(fit1$par[5001:iters,1], probs = c(.025, .25, .5, .75, .975))
quantile(exp(fit1$par[5001:iters,2]), probs = c(.025, .25, .5, .75, .975))
quantile(theta_all, probs = c(.025, .25, .5, .75, .975))
 

#Reference
#https://jytan.net/blog/2023/eight-schools/
#https://sites.stat.columbia.edu/gelman/book/data/
