# Example 5.5.3
# Metropolis-Hastings (Random walk)
# Sample from T distribution
# f(x) = (1 + x^2/v)^(-(v+1)/2)
# The proposal distribution is
# N(x,sigma^2)
library(ggplot2)

rw.Metropolis <- function(n, sigma, x0, v){
  # initialization
  x.star <- numeric(n)
  x.star[1] <- x0
  
  # Iteration
  accept <- 0
  for(i in 2:n){
    x <- x.star[i-1]
    y <- rnorm(1, x, sigma)
    ratio <- dt(y, v) / dt(x, v)
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
  Accept.rate <- accept / (n-1)
  
  return(list(x.star, Accept.rate))
}

n <- 1e5
sigma <- c(0.5, 1, 2.5, 10)
x0 <- 25
v <- 4
res <- lapply(sigma, function(i){
  rw.Metropolis(n, i, x0, v)
})

sapply(1:4, function(i) res[[i]][[2]])
## [1] 0.8557586 0.7262973 0.4711347 0.1543615


# Final for sigma = 2.5
x.star <- res[[3]][[1]]
x.final <- sample(x.star[(n/2):n], size = 1000)
ggplot(data = data.frame(x = x.final), aes(x = x)) +
  geom_histogram() + 
  theme_classic()
