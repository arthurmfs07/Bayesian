# ctrl + shift + c comments

# library(cmdstanr)
# 
# stan_file <- "C:/Users/arthu/OneDrive/Área de Trabalho/project/rf/mcmcstan.stan"
# 
# # Load and compile the Stan model
# mod <- cmdstan_model(stan_file)
# 
# set.seed(42)
# n = 9
# x = ramarendra(n, 1)
# 
# 
# data_list <- list(
#   N = n,
#   y = x
# )
# 
# 
# fit <- mod$sample(
#   data = data_list,       # Data for the model
#   seed = 42,             # Seed for reproducibility
#   chains = 4,             # Number of MCMC chains
#   parallel_chains = 4,    # Number of chains to run in parallel
#   iter_warmup = 500,      # Number of warm-up iterations
#   iter_sampling = 1000    # Number of sampling iterations
# )
# 
# # Print the results
# print(fit)

 
# fit_summary <- fit$summary()
# print(fit_summary)
# 
# # Check for divergences
# fit$cmdstan_diagnose()
# 
# 
# posterior_samples <- fit$draws(format = "df")
# 
# # Extract samples for a specific parameter (e.g., theta)
# theta_samples <- posterior_samples$theta
# 
# # Extract all chains separately
# all_chains <- fit$draws(format = "array")
# 
# library(coda)
# mcmc_samples <- as.mcmc(theta_samples)
# autocorr.plot(mcmc_samples, main = "Autocorrelation for theta")
# 
# 
# acf(theta_samples, main = "Autocorrelation Function for theta")
# 
# 
# library(bayesplot)
# mcmc_trace(fit$draws("theta"))
# mcmc_acf(fit$draws("theta"))
# 
# theta_summary <- fit$summary(variables = "theta")
# print(theta_summary)
# 

