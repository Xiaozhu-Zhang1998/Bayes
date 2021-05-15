# Gelman.Rubin
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i, 1:j])
  # for chan in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, var)
  W <- mean(psi.w)
  v.hat <- W * (n-1) / n + (B / (n*k))
  r.hat <- v.hat / W
  return(r.hat)
}

normal.chain <- function(sigma, N, X1){
  # generates a Metropolis chain for N(0,1)
  # with Normal(X[t], sigma) proposal distribution
  # and starting value X1
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for(i in 2:N){
    xt <- x[i-1]
    y <- rnorm(1, xt, sigma) # candidate point
    r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
    r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
    r <- r1 / r2
    if(u[i] <= r) x[i] <- y
    else x[i] <- xt
  }
  return(x)
}

# GO! ====
sigma <- 0.2  # parameter of proposal distribution
k <- 4        # number of chains to generate
n <- 15000    # length of chains
b <- 1000     # burn-in length

# choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

# generate the chains
X <- matrix(0, nrow = k, ncol = n)
for(i in 1:k)
  X[i,] <- normal.chain(sigma, n, x0[i])

# compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for(i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

# plot psi for the four chains
par(mfrow = c(2,2))
for(i in 1:k)
  plot(psi[i, (b+1):n], type = 'l', xlab = i, ylab = bquote(psi))
par(mfrow = c(1,1)) # restore default

# plot the sequence of R-hat statitics
rhat <- rep(0, n)
for(j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[, 1:j])
plot(rhat[(b+1):n], type = 'l', xlab = '', ylab = "R")
abline(h = 1.1, lty = 2)