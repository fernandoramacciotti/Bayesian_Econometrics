# Monte Carlo Methods with R
# Robert & Casella
# Cap 6


# Exercise 6.1 ------------------------------------------------------------

# Markov Chain defined as X_t+1 = rho * X_t + e_t
# e ~ N(0, 1)
# X_0 ~ N(0, 1)
# Plot the histogram for t <= 1e4 and rho = 0.9
# Check the potential fit to the stationary dist N(0, 1/(1 - rho^2))

library(ggplot2)

t <- 1e4
rho <- 0.9
x <- rep(0, t)

set.seed(10)
x[1] <- rnorm(1, mean = 0, sd = 1)


for(i in 2:t){
  e <- rnorm(1, mean = 0, sd = 1)
  x[i] <- rho * x[i-1] + e
}

dens <- dnorm(x, mean = 0, sd = sqrt(1 / (1 - rho^2)))
df <- data.frame(x = x, dens = dens)

ggplot(data = df, aes(x = x)) +
  geom_histogram(aes(y = ..density..)) + 
  geom_line(data = df, aes(x = x, y = dens))




# Example 6.1 -------------------------------------------------------------

# target density f is Beta(2.7, 6.3)
# proposed q is U[0,1]

rm(list = ls())

n <- 5e3
alpha <- 2.7
beta <- 6.3
f <- function(x, alpha, beta) dbeta(x, shape1 = alpha, shape2 = beta) # target distribution beta
q <- function(k) runif(k) # proposed distribution uniform
x <- rep(NA, n) # initializa chain
accepted <- rep(0, n) # vector to calculate acceptance ratio
x[1] <- runif(1)

for(i in 2:n){
  x_prop <- q(1)
  log_alpha <- log(f(x_prop, alpha, beta)) - log(f(x[i-1], alpha, beta)) # calculating log_alpha due to computational efficiency
  rho <- min(1, exp(log_alpha))
  
  u <- runif(1)
  
  if(rho > u){
    x[i] <- x_prop
    accepted[i] <- 1
  } 
  else {
    x[i] = x[i-1]
  }
}

acc_ratio <- mean(accepted) # acceptance ratio
dens <- f(x, alpha, beta) # calculating the target distribution for the generated chain

df <- data.frame(x = x, dens = dens)

# plotting

ggplot(data = df, aes(x = x)) + 
  geom_histogram(aes(y = ..density..)) +
  geom_line(aes(y = dens), colour = 'red')

# Kolmogorov-Smirnov test for sample equality

ks.test(jitter(x), rbeta(n, alpha, beta))

# Moments MC and analytical
E_beta <- alpha / (alpha + beta)
var_beta <- alpha * beta / (alpha + beta)^2 / (alpha + beta + 1)

E_MCMC <- mean(x)
var_MCMC <- var(x)

tab <- matrix(c(E_beta, var_beta, E_MCMC, var_MCMC), ncol = 2, nrow = 2)
colnames(tab) <- c("Real", "MCMC")
rownames(tab) <- c("E", "Var")
tab