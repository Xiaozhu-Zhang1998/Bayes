# Example 5.5.4: A investment model based on Bayesian Inference

# Generate the observations ----
b <- 0.2          # actual value of the beta
w <- 0.25         # width of the uniform support set
m <- 5000         # length of the chain
burn <- 1000      # burn-in time
days <- 250       
x <- numeric(m)   # the chain

# generate the observed frequencies of winners
i <- sample(1:5, size = days, replace = TRUE,
            prob = c(1, 1 - b, 1 - 2 * b, 2 * b, b))
win <- tabulate(i)
print(win)


# Metropolis ----
# computes (without the constant) the target density
prob <- function(y, win){
  if(y < 0 || y >= 0.5)
    return(0)
  return((1/3) ** win[1] * ((1-y)/3) ** win[2] * ((1-2*y)/3) ** win[3] *
           ((2*y)/3) ** win[4] * (y/3) ** win[5])
}

# random walk metropolis algorithm
u <- runif(m)         # for accept / reject step
v <- runif(m, -w, w)  # proposal distribution
x[1] <- 0.25
for(i in 2:m){
  y <- x[i-1] + v[i]
  if(u[i] <= prob(y, win) / prob(x[i-1], win))
    x[i] <- y
  else
    x[i] <- x[i-1]
}


# Plot ----
par(mfrow = c(1,2))
plot(x, type = 'l')
abline(h = b, v = 501, lty = 3)
xb <- x[-(1:501)]
hist(xb, probability = TRUE, xlab = bquote(beta), ylab = 'X', main = '')
z <- seq(min(xb), max(xb), length = 100)
lines(z, dnorm(z, mean(xb), sd(xb)))

# Estimate ----
print(win)
print(round(win / days, 3))
xb <- x[(burn + 1) : m]
print(mean(xb))
print(sd(xb))
