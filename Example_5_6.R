# Exercise 5.5.6
# Solution 1: random walk
wais <- read.table("wais.txt", header = T)
x <- wais[,1]
y <- wais[,2] 
m <- 55000                 # length of chain
mu.beta <- c(0, 0)         # prior dist
sigma.beta <- c(100, 100)  # prior dist
prop.s <- c(0.1, 0.1)      # proposal distribution standard variance
beta <- matrix(nrow = m, ncol = 2)
acc.prob <- 0
current.beta <- c(0,0)     # init beta
for (t in 1:m){
  prop.beta <- rnorm(2, current.beta, prop.s)       # candidate point
  cur.eta <- current.beta[1] + current.beta[2] * x  
  prop.eta <- prop.beta[1] + prop.beta[2] * x
  loga <- -(sum(y * prop.eta - log(1 + exp(prop.eta)))
            - sum(y * cur.eta - log(1 + exp(cur.eta)))
            + sum(dnorm(prop.beta, mu.beta, sigma.beta, log = T))
            - sum(dnorm(current.beta, mu.beta, sigma.beta, log = T)))
  if(is.na(loga)) loga <- 0
  u <- runif(1)
  u <- log(u)
  if(u < loga){
    current.beta <- prop.beta
    acc.prob <- acc.prob + 1
  }
  beta[t, ] <- current.beta
}
acc.prob <- acc.prob / m
print(acc.prob)


# convergence diagnostics plot
erg.mean <- function(x){
  n <- length(x)
  result <- cumsum(x) / cumsum(rep(1, n))
}

burnin <- 15000
idx <- seq(1, m, 50)
idx2 <- seq(burnin + 1, m)
par(mfrow = c(2,2))
plot(idx, beta[idx, 1], type = 'l', xlab = 'Iterations', ylab = 'Values of beta0')
plot(idx, beta[idx, 2], type = 'l', xlab = 'Iterations', ylab = 'Values of beta1')

ergbeta0 <- erg.mean(beta[, 1])
ergbeta02 <- erg.mean(beta[idx2, 1])
ylims0 <- range(c(ergbeta0, ergbeta02))

ergbeta1 <- erg.mean(beta[, 2])
ergbeta12 <- erg.mean(beta[idx2, 2])
ylims1 <- range(c(ergbeta1, ergbeta12))

plot(idx, ergbeta0[idx], type = 'l', xlab = 'Iterations', ylab = 'Values of beta0',
     main = '(c) Ergodic Mean Plot of beta0', ylim = ylims0)
lines(idx2, ergbeta02[idx2 - burnin], col = 2, lty = 2)

plot(idx, ergbeta1[idx], type = 'l', xlab = 'Iterations', ylab = 'Values of beta1',
     main = '(d) Ergodic Mean Plot of beta1', ylim = ylims1)
lines(idx2, ergbeta12[idx2 - burnin], col = 2, lty = 2)
apply(beta[(burnin + 1):m, ], 2, mean)
apply(beta[(burnin + 1):m, ], 2, sd)