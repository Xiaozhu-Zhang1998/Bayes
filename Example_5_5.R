# Independent sampling
# Basic Parameters 
m <- 5000         # length of chain
xt <- numeric(m)
a <- 1            # parameter of Beta(a,b)
b <- 1            # parameter of Beta(a,b)
p <- 0.2          # mixing parameter
n <- 30           # sample size
mu <- c(0, 5)     # parameters of the normal densities
sigma <- c(1, 1)

# Generate observations ----
i <- sample(1:2, size = n, replace = TRUE, prob = c(p, 1-p))
x <- rnorm(n, mu[i], sigma[i])

# Metropolis-Hastings Chain ----
u <- runif(m)
y <- rbeta(m, a, b)
xt[1] <- 0.5
for(i in 2:m){
  fy <- y[i] * dnorm(x, mu[1], sigma[1]) + 
    (1 - y[i]) * dnorm(x, mu[2], sigma[2])
  fx <- xt[i-1] * dnorm(x, mu[1], sigma[1]) + 
    (1 - xt[i-1]) * dnorm(x, mu[2], sigma[2])
  r <- prod(fy/fx) * (xt[i-1] ** (a-1) * (1 - xt[i-1]) ** (b-1)) /
    (y[i] ** (a-1) * (1 - y[i]) ** (b-1))
  if(u[i] <= r) xt[i] <- y[i]
  else xt[i] <- xt[i-1]
}

# Plot ----
plot(xt, type = 'l', ylab = 'p')
hist(xt[101:m], main = '', xlab = 'p', prob = TRUE)
print(mean(xt[101:m]))