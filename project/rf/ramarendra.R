ramarendra <- function(n, theta) {
  denom <- theta^3 + theta^2 + 2 * theta + 6
  p1 <- theta^3 / denom     
  p2 <- theta^2 / denom     
  p3 <- 2 * theta / denom   
  p4 <- 6 / denom           
  
  u <- runif(n)
  mixture_samples <- numeric(n)
  
  for (i in 1:n) {
    if (u[i] < p1) {
      mixture_samples[i] <- rexp(1, rate = theta)
    } else if (u[i] < (p1 + p2)) {
      mixture_samples[i] <- rgamma(1, shape = 2, rate = theta)
    } else if (u[i] < (p1 + p2 + p3)) {
      mixture_samples[i] <- rgamma(1, shape = 3, rate = theta)
    } else {
      mixture_samples[i] <- rgamma(1, shape = 4, rate = theta)
    }
  }
  return(mixture_samples)
}
