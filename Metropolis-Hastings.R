# Example 5.5.2
# Metropolis-Hastings
# Sample from Rayleigh distribution
# f(x) = x/sigma^2 * exp(-x^2 / 2*sigma^2)
# The proposal distribution is
# chi-square(x[i-1])
library(ggplot2)

n <- 1e5
x.star <- numeric(n)
f <- function(x){
  x/sigma^2 * exp(-x^2 / 2*sigma^2)
}
sigma <- 4

# initialization
set.seed(3214)
x.star[1] <- rchisq(1, df = 1)

# Iteration
accept <- 0
for(i in 2:n){
  x <- x.star[i-1]
  y <- rchisq(1, df = x)
  ratio <- (f(y) * dchisq(x, df = y)) / (f(x) * dchisq(y, df = x))
  u <- runif(1)
  if(u <= ratio){
    x.star[i] <- y
    accept <- accept + 1
  }
  else{
    x.star[i] <- x
  }
}

# Accept rate
accept / (n-1)

# Vis
index <- 10000:15000
ggplot(data = data.frame(x = index, y = x.star[index]), aes(x = x, y = y)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  xlab("Iteration") +
  ylab("Sample")

# Final
x.final <- sample(x.star[(n/2):n], size = 1000)
ggplot(data = data.frame(x = x.final), aes(x = x)) +
  geom_histogram() + 
  theme_classic()
