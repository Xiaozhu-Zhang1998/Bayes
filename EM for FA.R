# EM algorithm for Factor Analysis
library(ggplot2)
library(mvtnorm)
## The true model ----
z <- rnorm(n = 1000, mean = 0, sd = 1)
lambda <- c(10, 15)
mu <- c(3,5)
phi <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
epsilon <- rmvnorm(n = 1000, mean = c(0,0), sigma = phi)
dat <- t(sapply(1:1000, function(i){
  mu + lambda * z[i] + epsilon[i,]
}))

ggplot(data.frame(x = dat[,1], y = dat[,2]), aes(x = x, y = y)) + 
  geom_point() + 
  theme_bw()

## E-step
Estep <- function(dat, mu, lambda, phi){
  Ez <- apply(dat, 1, function(x){
    t(lambda) %*% solve(lambda %*% t(lambda) + phi) %*% (x - mu)
  })
  Ezz <- sapply(Ez, function(x){
    1 - t(lambda) %*% solve(lambda %*% t(lambda) + phi) %*%
      lambda + x %*% t(x)  
  })
  return(list(Ez = Ez, Ezz = Ezz))
}

## M-step
Mstep <- function(dat, Ez, Ezz){
  # obtain lambda
  temp1 <- t(sapply(1:N, function(i){
    (dat[i,] - mu) * Ez[i]
  }))
  temp2 <- apply(temp1, 2, sum)
  temp3 <- 1 / sum(Ezz)
  lambda <- temp2 * temp3
  # obtain phi
  temp4 <- t(sapply(1:N, function(i){
    lambda %*% matrix(Ezz[i]) %*% t(lambda)
  }))
  temp5 <- matrix(apply(temp4, 2, sum)/N, 2)
  phi <- cov(dat) - temp5
  return(list(lambda = lambda, phi = phi))
}

# start to go!
## Initialization ----
lambda <- c(12,18)
phi <- diag(2)

## GO
N <- dim(dat)[1]
mu <- apply(dat, 2, sum) / N
for (i in 1:50){
  # Estep
  par1 <- Estep(dat, mu, lambda, phi)
  Ez <- par1$Ez
  Ezz <- par1$Ezz
  # Mstep
  par2 <- Mstep(dat, Ez, Ezz)
  lambda <- par2$lambda
  phi <- par2$phi
}

## log-likelihood
A <- -0.5 * N * log(det(lambda %*% t(lambda) + phi))
B <- -0.5 * sum(sapply(1:N, function(i){
  t(dat[i,] - mu) %*% solve(lambda %*% t(lambda) + phi) %*% (dat[i,] - mu)
}))
loglikelihood <- A + B


## Vis
slope <- lambda[2]/lambda[1]
intercept <- mu[2] - slope * mu[1]

ggplot(data.frame(x = dat[,1], y = dat[,2]), aes(x = x, y = y)) + 
  geom_point() + 
  geom_abline(slope = slope, intercept = intercept, 
              color = 'red', size = 1) +
  theme_bw()

