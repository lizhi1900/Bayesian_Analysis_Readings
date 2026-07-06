rm(list=ls(all=TRUE))
#https://www.pymc.io/projects/examples/en/latest/mixture_models/dp_mix.html
#zli3 **at** live &&dot&& com
library(ggplot2)
library(dplyr)
library(rjags)
library(tidybayes)
library(rmp)

data(ethylene)
mydata=unique(ethylene[,-c(3,4,6,7)])
y=mydata[which(mydata$dose==0),]$impl

N <- length(y )
K <- 30


###########################################################
model_code <- "
model {
  for (i in 1:N) {
    y[i] ~ dpois(lambda[zeta[i]])
    zeta[i] ~ dcat(pi[])
  }
  for (h in 1:H) {
    lambda[h] ~ dgamma(25, 2)
  }
  
  # Stick breaking prior
  for (h in 1:(H-1)) {
    V[h] ~ dbeta(1,a)
  }
  V[H]=1
  pi[1] <- V[1]
  for (h in 2:H) {
    pi[h] <- V[h] * (1-V[h-1]) * pi[h-1] / V[h-1]
  }
  a=2#~dgamma(1,1)
}
"

data_list <- list(
  N = N,
  H = K,
  y = y 
)
model <- jags.model(textConnection(model_code), 
                    data=data_list, 
                    n.chains = 5);

update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                                 variable.names = c("pi","lambda"),
                                 n.iter = 100, progress.bar = "none")

#summary(posterior_sample)

pi_dat<-posterior_sample %>% spread_draws(pi[i]|i)#tidybayes
lambda_dat<-posterior_sample %>% spread_draws(lambda[i]|i)#tidybayes

# Plot 1: Bar plot of posterior expected mixture weight
plot_w <- 1:K
w_mean<-colMeans(pi_dat[4:(3+K)])
df_bar <- data.frame(Component = plot_w, Weight = w_mean)

ggplot(df_bar, aes(x = Component, y = Weight)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  scale_x_continuous(limits = c(0.5, K)) +
  labs(x = "Component", y = "Posterior expected mixture weight") +
  theme_minimal()


theta <- as.matrix(lambda_dat[,4:(3+K)])
f <- function(x, theta) dpois(x, theta)

# Assuming x_plot is defined somewhere, for example:
x_plot <- 0:50

# Compute dpm_pdf_components: array of dimension N x K x length(x_plot)
dpm_pdf_components <- array(NA, dim = c(nrow(theta ), ncol (theta ), length(x_plot)))
for (i in 1:nrow(theta )) {
  for (j in 1:ncol (theta )) {
    dpm_pdf_components[i, j, ] <- f(x_plot, theta[i, j])
  }
}

# Compute weighted sum over K components for each N and x_plot
w=as.matrix(pi_dat[,4:(3+K)])
dpm_pdfs <- matrix(0, nrow = nrow(theta ), ncol = length(x_plot))
for (i in 1:nrow(theta )) {
  for (k in 1:length(x_plot)) {
    dpm_pdfs[i, k] <- sum(w[i, ] * dpm_pdf_components[i, , k])
  }
}


ggplot() +
  # Histogram of standardized waiting times
  geom_histogram(
    aes(x = y ,y =after_stat(density)),
    bins = 20,
    color = "NA",
    fill = "steelblue",
    alpha = 0.5,
    position = "identity"
  ) +

  # Plot posterior sample densities
  geom_line(
    aes(x = rep(x_plot, times = nrow(theta )), y = as.vector(t(dpm_pdfs)),group = rep(1:nrow(theta ), each = length(x_plot)),linetype = "Posterior sample densities" ),
    color = "gray",
    alpha = 0.5
  ) +
 
  # Plot posterior expected density
  geom_line(
    aes(x = x_plot, y =colMeans(dpm_pdfs),linetype = "Posterior expected density" ),
    color = "black"
  ) + labs(linetype ="")+ #"Reference Lines"

  # Labels and theme
  labs(
    x = "Control Group Implantation of Mice",
    y = "Density"
  )+theme(legend.position = "top")






