# Rejection Sampling

rm(list = ls())

n <- 1e5

M <- 10
alpha <- 4
beta <- 6
x <- runif(n) # generating random sampling
f <- dbeta(x, shape1 = alpha, shape2 = beta) # target density
g <- dnorm(x) # sampling density

c <- f / (M * g)
result <- x[runif(n) < c]
acceptance_ratio <- length(result) / n

E_beta <- alpha / (alpha + beta)
var_beta <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))

E_MC <- mean(result)
var_MC <- var(result)

tab <- matrix(c(E_beta, var_beta, E_MC, var_MC), nrow = 2, ncol = 2)
rownames(tab) <- c('Expected Value', 'Var')
colnames(tab) <- c('Real', 'MC')
tab


df = data.frame(x, f, g)
ggplot(data = data.frame(result), aes(result)) + 
  geom_histogram(aes(y = ..density..), bins = 25) + 
  geom_line(data = df, aes(x = x, y = f, colour = 'Target Density')) + 
  geom_line(data = df, aes(x = x, y = M * g, colour = 'Proposal (M*g)')) + 
  ggtitle(paste('Acceptance Ratio: ', acceptance_ratio)) + 
  xlab('x')
