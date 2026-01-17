rm(list=ls(all=TRUE))
require(rjags)
library("ggdist")
#logit <- qlogis
invlogit <- plogis
set.seed(1)

set.seed(1234)
N <- 50
x <- rnorm(N)
a <- 0
b <- .8
p <- invlogit(a + b*x)
y <- rbinom(N, 1, p)
df <- 4
data_1 <- list(N=N, x=x, y=y, df=df)

#################################################################
#zli3@live.com

# Jags code to fit the robust probit model to the simulated data

robit_string <- "
model
{
	for (i in 1:n)
	{
		y[i] ~ dbern(p[i])
            p[i] <- pt(z[i], 0,1,1)
		z[i] <- a  + b*x[i]
	}
	
	a ~ dnorm(0, 0.0001)
	b ~ dnorm(0, 0.0001)
}
"

# Set up the data
model_data <- list(n = N, y = y, x = x, df=df)

# Choose the parameters to watch
model_parameters <- c("a", "b" )


model <- jags.model(textConnection(robit_string), 
                    data=model_data, 
                    n.chains = 5)

update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = model_parameters, 
                                 n.iter = 100, progress.bar = "none")

summary(posterior_sample)