# --------------------- functions --------------------- #

pdf <- function(x, par) {
  theta <- par[1]
  fx = theta^4 / (theta^3 + theta^2 + 2 * theta + 6) * (1 + x + x^2 + x^3) * exp(-theta * x)
  return(fx)
}


cdf <- function(x, par) {
  theta <- par[1]
  Fx = 1-(1+((x^3*theta^3+theta^2*(theta+3)*x^2+theta*(theta^2+2*theta+6)*x)/(theta^3+theta^2+2*theta+6)))*exp(-theta*x)
  return(Fx)
}

hrf <- function(x, par){
  theta <- par[1]
  numerator <- theta^4 * (1 + x + x^2 + x^3)
  denominator <- theta^3 * (1 + x + x^2 + x^3) + theta^2 * (3 * x^2 + 2 * x + 1) + 2 * theta * (3 * x + 1) + 6
  hrf <- numerator / denominator
  return(hrf)
}


# --------------------- plots --------------------- #

library(ggplot2)


plot_pdf <- function(space, theta_values) {
  data_pdf <- data.frame(
    x = rep(space, length(theta_values)),
    theta = factor(rep(theta_values, each = length(space))),
    y = unlist(lapply(theta_values, function(theta) pdf(space, theta)))
  )
  
  ggplot(data_pdf, aes(x = x, y = y, color = theta)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + 
    geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +  
    scale_color_brewer(palette = "Set1") +
    labs(color = expression(theta), title = "PDF Amarendra Distribution", x = "Values", y = "Density") +
    theme_minimal()
}


plot_cdf <- function(space, theta_values) {
  data_cdf <- data.frame(
    x = rep(space, length(theta_values)),
    theta = factor(rep(theta_values, each = length(space))),
    y = unlist(lapply(theta_values, function(theta) cdf(space, theta)))
  )
  
  ggplot(data_cdf, aes(x = x, y = y, color = theta)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + 
    geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +  
    scale_color_brewer(palette = "Set1") +
    labs(color = expression(theta), title = "CDF Amarendra Distribution", x = "Values", y = "Density") +
    theme_minimal()
}


plot_hrf <- function(space, theta_values) {
  data_hrf <- data.frame(
    x = rep(space, length(theta_values)),
    theta = factor(rep(theta_values, each = length(space))),
    y = unlist(lapply(theta_values, function(theta) hrf(space, theta)))
  )
  
  ggplot(data_hrf, aes(x = x, y = y, color = theta)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + 
    geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +  
    scale_color_brewer(palette = "Set1") +
    labs(color = expression(theta), title = "hrf Amarendra Distribution", x = "Values", y = "Density") +
    theme_minimal()
}


plot_mix <- function(n, space, theta_values) {
  data_dist <- data.frame(theta = factor(rep(theta_values, each = length(space))),
                          xvalues = unlist(lapply(theta_values, function(theta) ramarendra(n, theta))))
  
  
  ggplot(data_dist, aes(x = xvalues, fill = theta, color = theta)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2) +
    geom_density(aes(color = theta), alpha = 0, size = 1) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    labs(color = "Theta", fill = "Theta", title = "Empirical Density and Histograms", x = "Values", y = "Density") +
    theme_light()
  
}


# --------------------- plots for fixed parameter simulation --------------------- #


plot_fixed_bias <- function(){
  ggplot(bias_results, aes(x = SampleSize, y = Bias)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) + 
    labs(
      title = paste("Bias Evolution for Fixed Parameter =", fixed_param),
      x = "Sample Size",
      y = "Bias"
    ) +
    theme_minimal()
}



plot_fixed_mse <- function(){
  ggplot(mse_results, aes(x = SampleSize, y = MSE)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) + 
    labs(
      title = paste("MSE Evolution for Fixed Parameter =", fixed_param),
      x = "Sample Size",
      y = "MSE"
    ) +
    theme_minimal()
}

