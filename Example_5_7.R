# initialize constants and parameters
N <- 5000               # length of chain
burn <- 1000            # burn-in length
X <- matrix(0, N, 2)    # the chain, a bivariate sample
rho <- 0.75             # correlation
mu1 <- 0; mu2 <- 2; sigma1 <- 1; sigma2 <- 0.5
s1 <- sqrt(1 - rho ** 2) * sigma1
s2 <- sqrt(1 - rho ** 2) * sigma2

# generate chain
x[1,] <- c(mu1, mu2)    # initialize
for(i in 2:N){
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1 / sigma2
  X[i,1] <- rnorm(1, m1, s1)
  x1 <- X[i,1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2 / sigma1
  X[i,2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N,]

# compare sample statistics to parameters
colMeans(x)
cov(x)
plot(x, main = "", cex = 0.5, xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(x[,2]))
