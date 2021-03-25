# EM algorithm for GMM
library(ggplot2)
## The true model ----
c1 <- rnorm(n = 1000, mean = 10, sd = 3)
c2 <- rnorm(n = 1000, mean = 20, sd = 3.5)
c3 <- rnorm(n = 1000, mean = 30, sd = 4)
dat <- c(c1, c2, c3)
dat <- sample(dat)
ggplot(data.frame(x = dat), aes(x = x)) + 
  geom_density(size = 1) + 
  theme_bw()

## E-step
Estep <- function(dat, mu1, mu2, mu3, 
                  sigma1, sigma2, sigma3,
                  pi1, pi2, pi3){
  w1 <- sapply(dat, 
               function(x){dnorm(x, mean = mu1, sd = sigma1) * pi1})
  w2 <- sapply(dat, 
               function(x){dnorm(x, mean = mu2, sd = sigma2) * pi2})
  w3 <- sapply(dat, 
               function(x){dnorm(x, mean = mu3, sd = sigma3) * pi3})
  w <- cbind(w1, w2, w3)
  deno <- apply(w, 1, sum)
  w <- t(sapply(1:dim(w)[1], function(i){w[i,]/deno[i]}))
  return(w)
}

## M-step
Mstep <- function(dat, w){
  mu1 <- sum(w[,1] * dat) / sum(w[,1])
  mu2 <- sum(w[,2] * dat) / sum(w[,2])
  mu3 <- sum(w[,3] * dat) / sum(w[,3])
  sigma1 <- sqrt(sum(w[,1] * (dat-mu1)^2) / sum(w[,1]))
  sigma2 <- sqrt(sum(w[,2] * (dat-mu2)^2) / sum(w[,2]))
  sigma3 <- sqrt(sum(w[,3] * (dat-mu3)^2) / sum(w[,3]))
  pi1 <- sum(w[,1])/dim(w)[1]
  pi2 <- sum(w[,2])/dim(w)[1]
  pi3 <- sum(w[,3])/dim(w)[1]
  return(list(mu1 = mu1, mu2 = mu2, mu3 = mu3,
              sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3,
              pi1 = pi1, pi2 = pi2, pi3 = pi3))
}

# start to go!
## Initialization ----
pi1 = 0.2; pi2 = 0.3; pi3 = 0.5
mu1 = 10; mu2 = 22; mu3 = 29
sigma1 = 2; sigma2 = 2; sigma3 = 2

## GO
for (i in 1:100){
  w <- Estep(dat, mu1, mu2, mu3, sigma1, sigma2, sigma3,
             pi1, pi2, pi3)
  par <- Mstep(dat, w)
  # repar
  mu1 <- par$mu1; mu2 <- par$mu2; mu3 <- par$mu3
  sigma1 <- par$sigma1; sigma2 <- par$sigma2; sigma3 <- par$sigma3
  pi1 <- par$pi1; pi2 <- par$pi2; pi3 <- par$pi3
}

## likelihood
mu <- c(mu1, mu2, mu3)
sigma <- c(sigma1, sigma2, sigma3)
pi <- c(pi1, pi2, pi3)

likelihood <-sum(sapply(1:N, function(i){
  sum(sapply(1:3, function(j){
    dnorm(dat[i], mean = mu[label[j]], sd = sigma[label[j]]) * pi[[label[j]]]
  }))
}))
likelihood

## Vis
N <- length(dat)
label <- apply(w, 1, function(x){which(x == max(x))})
dat.n <- data.frame(x = dat, label = as.factor(label))
ggplot(data.frame(x = dat), aes(x = x)) + 
  geom_histogram(data = dat.n, 
                 aes(x = x, y =..density.., fill = label)) +
  theme_bw()

