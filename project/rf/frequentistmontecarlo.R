run_mle <- function(x, param) {
  nloglike_mc <- function(l, x) {
    n <- length(x)
    -1*(n*log(l^4) - n*log(l^3 + l^2 + 2*l + 6) - l*sum(x) + sum(log(1 + x + x^2 + x^3)))
  }
  result <- optim(1, nloglike_mc, x = x, method = "L-BFGS-B", lower = 0.01, upper = 100, hessian = TRUE)
  se <- sqrt(diag(solve(result$hessian)))
  ci <- result$par + qnorm(c(0.025, 0.975)) * se
  list(mle = result$par, hessian = result$hessian, ci = ci)
}

run_simulation <- function(n_mc_values, parameters, n_sims) {
  results_all <- list()
  
  for (n_mc in n_mc_values) {
    results <- lapply(parameters, function(param) {
      
      sims_mc <- replicate(n_sims, ramarendra(n_mc, theta = param))
      
      
      mle_results <- lapply(1:n_sims, function(i) run_mle(sims_mc[, i], param))
      
      mle_estimates <- sapply(mle_results, `[[`, "mle")
      ci_estimates <- do.call(rbind, lapply(mle_results, `[[`, "ci"))
      
      bias <- mean(mle_estimates) - param
      mse <- mean((mle_estimates - param)^2)
      proportion_ci <- mean(apply(ci_estimates, 1, function(ci) param >= ci[1] && param <= ci[2]))
      
      list(
        sims_mc = sims_mc,
        mle_estimates = mle_estimates,
        ci_estimates = ci_estimates,
        bias = bias,
        mse = mse,
        proportion_ci = proportion_ci
      )
    })
    
    names(results) <- as.character(parameters)
    results_all[[as.character(n_mc)]] <- results
  }
  
  return(results_all)
}



plot_mc_tables <- function() {
  for (param in parameters) {
    table_data_mc <- data.frame(
      Parameter = param,  
      SampleSize = n_mc_values,
      MLE_mean = sapply(results_mc_all, function(x) mean(x[[as.character(param)]]$mle_estimates)),
      Bias = sapply(results_mc_all, function(x) x[[as.character(param)]]$bias),
      MSE = sapply(results_mc_all, function(x) x[[as.character(param)]]$mse),
      PC = sapply(results_mc_all, function(x) x[[as.character(param)]]$proportion_ci)
    )
    
    print(table_data_mc)
  }
}


 # -------------------------- simulation for fixed param plot evolution in funcamarendra.R/ -------------------------- #


run_simulation_fixed_bias <- function(sample_sizes, fixed_param, n_sims) {
  results <- sapply(sample_sizes, function(n_mc) {
    
    sims_mc <- replicate(n_sims, ramarendra(n_mc, theta = fixed_param))
    
    mle_results <- lapply(1:n_sims, function(i) run_mle(sims_mc[, i], fixed_param))
    
    mle_estimates <- sapply(mle_results, `[[`, "mle")
    
    bias <- mean(mle_estimates) - fixed_param
    return(bias)
  })
  
  data.frame(SampleSize = sample_sizes, Bias = results)
}


run_simulation_fixed_mse <- function(sample_sizes, fixed_param, n_sims) {
  results <- sapply(sample_sizes, function(n_mc) {
    
    sims_mc <- replicate(n_sims, ramarendra(n_mc, theta = fixed_param))
    
    mle_results <- lapply(1:n_sims, function(i) run_mle(sims_mc[, i], fixed_param))
    
    mle_estimates <- sapply(mle_results, `[[`, "mle")
    
    mse <- mean((mle_estimates- fixed_param)^2)
    return(mse)
  })
  
  data.frame(SampleSize = sample_sizes, MSE = results)
}




library(ggplot2)

# plot_confidence_intervals(results_mc_all, n_mc = 50, theta = 0.5, n_sims = 50)

plot_confidence_intervals <- function(results_all, n_mc, theta, n_sims) {
  results <- results_all[[as.character(n_mc)]][[as.character(theta)]]
  ci_estimates <- results$ci_estimates
  
  plot_data <- data.frame(
    Simulation = 1:n_sims,  
    Lower = ci_estimates[, 1],
    Upper = ci_estimates[, 2]
  )
  
  plot_data$CoversTheta <- ifelse(
    theta >= plot_data$Lower & theta <= plot_data$Upper,
    "Covers",
    "Doesn't Cover"
  )
  
  plot_data <- plot_data[!duplicated(plot_data$Simulation), ]
  
  if (anyDuplicated(plot_data$Simulation) > 0) {
    stop("Duplicate simulation numbers found after deduplication.")
  }
  

  ggplot(plot_data, aes(x = Simulation)) +
    geom_segment(
      aes(
        xend = Simulation,
        y = Lower,
        yend = Upper,
        color = CoversTheta
      ),
      size = 1
    ) +
    geom_hline(
      yintercept = theta,
      linetype = "dashed",
      color = "red",
      size = 1
    ) +
    scale_color_manual(
      values = c("Covers" = "blue", "Doesn't Cover" = "orange"),
      labels = c("Covers Theta", "Doesn't Cover Theta")
    ) +
    labs(
      title = paste("Confidence Intervals for Theta =", theta, "with Sample Size =", n_mc),
      x = "Simulation Number",
      y = "Confidence Interval",
      color = "Confidence Interval Status"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      text = element_text(size = 12)
    )
}
