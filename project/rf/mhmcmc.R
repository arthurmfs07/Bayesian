# -------------------------------------------------------------------------------------------------------------------#
# Amarendra likelihood and prior gamma(alpha, beta)

library(coda)

# posterior kernel
m <- function(theta, n, x, alpha, beta) {
  
  term1 <- theta^(4*n+alpha-1)
  term2 <- (theta^3 + theta^2 + 2 * theta + 6)^(-n)
  term3 <- exp(-theta*(beta+sum(x)) )
  
  return(term1 * term2 * term3)
}

AMHtheta <- function(n, x, alpha, beta, niter, start) {
  if (missing(niter) || !is.numeric(niter) || niter <= 0) {
    stop("Error: 'niter' must be a positive numeric value.")
  }
  
  theta = numeric(niter)
  
  # Initialize
  theta[1] = start
  taxa = 0
  
  # Gamma proposal distribution 
  shape_prop = 14
  rate_prop = 14
  
  for (i in 2:niter) {
    
    # Step 2: Propose a new value
    thetac = rgamma(1, shape = shape_prop, rate = rate_prop)
    # Step 3: Compute acceptance ratio
    
    m_thetac = m(thetac, n, x, alpha, beta)
    m_theta_prev = m(theta[i - 1], n, x, alpha, beta)
    
    # Check for zero or invalid posterior values
    if (m_thetac == 0 || m_theta_prev == 0 || is.na(m_thetac) || is.na(m_theta_prev)) {
      alfa = 0 # Reject the proposed value
    } else {
      posterior_ratio = m_thetac / m_theta_prev
      
      # Proposal density ratio 
      proposal_ratio = dgamma(theta[i - 1], shape = shape_prop, rate = rate_prop) / 
        dgamma(thetac, shape = shape_prop, rate = rate_prop)
      
      # Full acceptance ratio
      alfa = min(1, posterior_ratio * proposal_ratio)
    }
    
    # Step 4: Accept or reject
    u = runif(1)
    if (u < alfa) { 
      theta[i] = thetac
      taxa = taxa + 1
    } else {
      theta[i] = theta[i - 1]
    }
  }
  
  taxa = taxa / niter 
  return(list(theta = theta, taxa = taxa))
}


plot_trace <- function(result1, result2){
  plot(result2, type = "l", main = "Trace of Theta", 
       ylab = "Theta", col = "red")
  lines(result1, col = "darkblue")
}


# summary(chain)
# plot(chain)
# autocorr.plot(chain)
# effectiveSize(chain)
# gelman.diag(mcmc.list(chain, chain2, chain3), autoburnin = FALSE)



# -------------------------------------------------------------------------------------------------------------------#
# Amarendra likelihood and jeffreys prior



m1 <- function(theta, n, x) {
  
  term1 <- theta^(4*n-1)
  term2 <- (theta^3 + theta^2 + 2 * theta + 6)^(-(n+1))
  term3 <- exp(-theta*sum(x))
  term4 <- sqrt(theta^6+4*theta^5+18*theta^4+96*theta^3+72*theta^2+96*theta+144)
  
  return(term1 * term2 * term3 * term4)
}

AMHjeffreys <- function(n, x, niter, start) {
  if (missing(niter) || !is.numeric(niter) || niter <= 0) {
    stop("Error: 'niter' must be a positive numeric value.")
  }
  
  theta = numeric(niter)
  
  # Initialize
  theta[1] = start
  taxa = 0
  
  # Gamma proposal distribution 
  shape_prop = 200
  rate_prop = 200
  
  for (i in 2:niter) {
    
    # Step 2: Propose a new value
    thetac = rgamma(1, shape = shape_prop, rate = rate_prop)
    
    # Step 3: Compute acceptance ratio
    
    m_thetac = m1(thetac, n, x)
    m_theta_prev = m1(theta[i - 1], n, x)
    
    # Check for zero or invalid posterior values
    if (m_thetac == 0 || m_theta_prev == 0 || is.na(m_thetac) || is.na(m_theta_prev)) {
      alfa = 0 # Reject the proposed value
    } else {
      posterior_ratio = m_thetac / m_theta_prev
      
      # Proposal density ratio 
      proposal_ratio = dgamma(theta[i - 1], shape = shape_prop, rate = rate_prop) / 
        dgamma(thetac, shape = shape_prop, rate = rate_prop)
      
      # Full acceptance ratio
      alfa = min(1, posterior_ratio * proposal_ratio)
    }
    
    # Step 4: Accept or reject
    u = runif(1)
    if (u < alfa) { 
      theta[i] = thetac
      taxa = taxa + 1
    } else {
      theta[i] = theta[i - 1]
    }
  }
  
  taxa = taxa / niter 
  return(list(theta = theta, taxa = taxa))
}


# res_test <- AMHjeffreys(n = n, x = x, niter = niter, start = 1)
# cat("Acceptance rate:", res_test$taxa, "\n")
# chain_test = mcmc(res_test$theta[(burnin + 1):niter])
# chain_test_thin <- window(chain_test, thin = 10)


# checking if integral is proper
# plot(seq(0.1, 5, length.out = 100), sapply(seq(0.1, 5, length.out = 100), function(theta) m1(theta, n, x)), type = "l", col = "blue", lwd = 2,
#      xlab = expression(theta), ylab = "m1(theta, n, x)",
#      main = "Plot of m1(theta, n, x)")
# integrate_function <- function(n, x) {
#   result <- integrate(m1, lower = 0, upper = Inf, n = n, x = x)
#   return(result)
# }
# integration_result <- integrate_function(n, x)
# cat("Numerical Integration Result:", integration_result$value, "\n")




