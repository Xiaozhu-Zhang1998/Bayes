# Exercise 5.5.6
# Solution 3: Coordinate M-H
wais <- read.table("wais.txt", header = T)

# Basic parameters
y <- wais[,2]
x <- wais[,1]
m <- 10000
beta0 <- c(0, 0)        # initial value
mu.beta <- c(0, 0)      # prior
s.beta <- c(100, 100)   # prior
prop.s <- c(1.75, 0.2)  # sd of proposal normal
beta <- matrix(nrow = m, ncol = 2)
acc.prob <- c(0, 0)
current.beta <- beta0

# Start
for(t in 1:m){
  for(j in 1:2){
    prop.beta <- current.beta
    prop.beta[j] <- rnorm(1, current.beta[j], prop.s[j])
    cur.eta <- current.beta[1] + current.beta[2] * x
    prop.eta <- prop.beta[1] + prop.beta[2] * x
    loga <- -(sum(y * prop.eta - log(1 + exp(prop.eta)))
              - sum(y * cur.eta - log(1 + exp(cur.eta)))
              + sum(dnorm(prop.beta, mu.beta, s.beta, log = T))
              - sum(dnorm(current.beta, mu.beta, s.beta, log = T)))
    if(is.na(loga)) loga <- 0
    u <- runif(1)
    u <- log(u)
    if(u < loga){
      current.beta <- prop.beta
      acc.prob[j] <- acc.prob[j] + 1
    }
    beta[t,] <- current.beta
  }
}
print(acc.prob / m)
  