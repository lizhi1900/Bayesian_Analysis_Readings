library(ggplot2)
library(tidyverse)
library(tidybayes)
require(MASS)
require(rjags)
rm(list=ls(all=TRUE))

set.seed(1234)
#https://juanitorduz.github.io/hsgp_intro/
#zli3 _a_t_ live com 
generate_synthetic_data <- function(start, stop, num, scale) {
  x <- seq(from = start, to = stop, length.out = num)
  y <- sin(4 * pi * x) + sin(7 * pi * x)
  y_obs <- y + rnorm(num, mean = 0, sd = scale)
  list(x = x, y = y, y_obs = y_obs)
}

n_train <- 20
n_test <- 30
scale <- 0.3

train_data <- generate_synthetic_data(start = 0, stop = 1, num = n_train, scale = scale)
test_data <- generate_synthetic_data(start = -0.2, stop = 1.2, num = n_test, scale = scale)

df_train <- data.frame(x = train_data$x, y = train_data$y, y_obs = train_data$y_obs, type = "train")
df_test <- data.frame(x = test_data$x, y = test_data$y, y_obs = test_data$y_obs, type = "test")

df_obs <- rbind(
  data.frame(x = df_train$x, y_obs = df_train$y_obs, type = "observed (train)"),
  data.frame(x = df_test$x, y_obs = df_test$y_obs, type = "observed (test)")
)

ggplot() +
  geom_point(data = df_obs, aes(x = x, y = y_obs, color = type)) +
  geom_line(data = df_train, aes(x = x, y = y), color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "blue", alpha = 0.3, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 1, color = "blue", alpha = 0.3, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("observed (train)" = "#1f77b4", "observed (test)" = "#ff7f0e")) +
  labs(
    x = "x",
    y = "y",
    color = NULL,
    title = "Synthetic Data"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal"
)
### get data and estimation years
x <- df_train %>% pull(x)
y <- df_train %>% pull(y_obs)
n_obs <- length(x)


 m=20 #very important to be correct
 L=1.5 #very important to be correct

 a=1.28
 l=.077

S=a^2*sqrt(2*pi)*l*exp(-.5*l^2*((1:m)*pi/2/L)^2)
phi=sqrt(1/L)*sin(outer(x+L,pi*1:m,"*")/2/L)
K_approx=phi%*%diag(sqrt(S))%*%t(phi)


plot_covariance <- function(matrix){
  
  mat_df <- as.data.frame(matrix)
  mat_df$x <- seq_len(nrow(mat_df))
  
  mat_df <- reshape(mat_df, 
                    direction = "long",
                    v.names = "values",
                    varying = 1:(ncol(mat_df) - 1),
                    times = 1:nrow(mat_df),
                    idvar = "x")
  
  p <- ggplot(data = mat_df) +
    geom_tile(aes(x = x, y = time, fill = values)) +
    labs(title = "Heatmap of covariance matrix",
         x = "x",
         y = "x",
         fill = "k(x,x')") +
    scale_fill_viridis_c(limits = c(0, 1),
                         breaks = seq(from = 0, to = 1, by = 0.2)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}

plot_covariance(matrix = K_approx)


##########################################################
gp_model <- '
model{

  for(i in 1:m){
    S[i]=a^2*sqrt(2*pi)*l^1*exp(-.5*l^2*((i)*pi/2/L)^2)
    beta[i] ~ dnorm(0,1);
  }

  gp <- phi %*% (sqrt(S)*beta );

  for(i in 1:n){
    y[i]~dnorm(gp[i],sigma_y^-2);
    y_pred[i]~dnorm(gp[i],sigma_y^-2);
  }
  a ~ dgamma(a1,b1)#amplitude 
  l ~ dgamma(a2,b2)#length_scale
  sigma_y ~  dt(0,10^-2,1)T(0,)#error term 

}'

data_list<-list(y = y, pi = pi, m = m, L=L,
                a1=8,b1=5,a2=1,b2=1,
                phi=phi,n=n_obs)
                  

##parameters to save
jags_pars <- c("a","l","sigma_y","y_pred") #
 
model <- jags.model(textConnection(gp_model), data=data_list, n.chains = 5)
    
update(model, 1000, progress.bar = "none"); # Burnin for 1000 samples

posterior_sample <- coda.samples(model, 
                    variable.names = jags_pars ,                   
                    n.iter = 400, 
                    progress.bar = "none")

par_dat <- posterior_sample %>% spread_draws(a,l,sigma_y)#tidybayes
summary(par_dat[])
#y_rep_mat=matrix(par_dat$y_pred,nrow=2000)
####################################################

out=summary(posterior_sample)

lower=out$quantiles[4:23,1]
upper=out$quantiles[4:23,5]
dat <- data.frame(df_train,lower=lower,upper=upper)

ggplot(data = dat,aes(x, y)) +
  geom_point(aes(x = x, y = y), colour = "black") +
  geom_line(aes(x = x, y = y), colour = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3,fill = "blue") +
  labs(x = "Timepoint",
       y = "Values") +
  theme_bw()

