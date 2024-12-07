library(R2jags)

#--------------------- Modelo Priori Jeffreys não informativa  ---------------------#



MCMC_gamma <- function(N_values, parameters, sims){

  cat("model{
	for(i in 1:N){


		# Likelihood

		L[i] <- ((theta^4)/(theta^3 + theta^2 + 2*theta + 6))*(1+y[i]+y[i]^2+y[i]^3)*exp(-theta*y[i])

		pro[i] <- L[i]/C
		ones[i] ~ dbern(pro[i])
	}
	
    theta ~ dgamma(0.01, 0.01)  

    C <- 1000000
}",
      file="gammanoinformative.jags")
  
  
    
  initialization <- function(){
    list("theta"=3)
  }
  
  df <- data.frame(theta = numeric(), N = numeric(), media = numeric(), 
                   vies = numeric(), eqm = numeric(), pc = numeric(),
                   stringsAsFactors = FALSE)
  
  S <- sims
  for (theta0 in parameters) {
    for (N in N_values) {
      theta_c=numeric()
      
      theta_inf=numeric()
      theta_sup=numeric()
      theta_vies=numeric()
      theta_eqm=numeric()
      
      for(j in 1:S){
        y <- ramarendra(N, theta0)
        
        ones <- rep(1, N)
        a.dat <- list("N", "y", "ones")
        
        parameters <- c("theta")
        
        model.fit <- jags(data = a.dat,
                 inits = initialization,
                 parameters.to.save = parameters,
                 n.chains = 1,
                 n.iter = 20000,
                 n.burnin = 5000,
                 n.thin = 10,
                 model.file ='gammanoinformative.jags')
        
        model.fit.mcmc <- as.mcmc(model.fit)
        thetac=as.matrix(model.fit.mcmc[,2])   
        
        
        theta_c[j]=mean(thetac)
        theta_inf[j]=quantile(thetac, 0.025)
        theta_sup[j]=quantile(thetac, 0.975)
        theta_vies[j]=theta_c[j]-theta0
        theta_eqm[j]=(theta_c[j]-theta0)^2
      }
      
      pc_theta=numeric()
      for(i in 1:S){
        if((theta_inf[i]<theta0) & (theta0<theta_sup[i])){
          pc_theta[i]=1
        } else {
          pc_theta[i]=0
        }
      }
      
      pc=mean(pc_theta)
      
      media=mean(theta_c)
      eqm=mean(theta_eqm)
      vies=mean(theta_vies)
      
      new_row <- c(theta0, N, media, vies, eqm, pc)
      df[nrow(df) + 1,] = new_row
    }
  }
  
  return(df)
  
  
}


table_mcmc <- function(df) {
  
  unique_thetas <- unique(df$theta)
  
  for (theta_value in unique_thetas) {
    df_filtered <- df %>%
      filter(theta == theta_value) %>%
      arrange(N)
    
    print(df_filtered, row.names = FALSE, align = "c")
  }
}

 
# denplot(model.fit.mcmc)
# 
# traplot(model.fit.mcmc)
# 
# mcmcplot(model.fit.mcmc)


#--------------------- Modelo Priori Empírica ---------------------#

log_likelihood_empirical <- function(params, x_data) {
  alpha <- params[1]  
  beta <- params[2]   
  
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  log_likelihood_value <- 0
  
  for (x in x_data) {
    integral_result <- tryCatch({
      integrate(
        function(theta) {
          theta^(alpha+3)*exp(-theta*(x+beta))/(theta^3+theta^2+2*theta+6)
        },
        lower = 0,
        upper = Inf
      )$value
    }, error = function(e) {
      return(0)  
    })
    
    log_likelihood_value <- log_likelihood_value+log(beta^alpha)-lgamma(alpha)+log(1+x+x^2+x^3)+log(integral_result)
  }
  
  return(log_likelihood_value)
}


library(nloptr)

# x0 = c(/,/) # lb = c(/,/) # ub = c(/,/)  
estimates_empirical_prior <- function(x_data, x0, lb, ub){
  
  objective_fn <- function(params) {
    -log_likelihood_empirical(params, x_data) 
  }
  
  result <- nloptr(
    x0 = x0,  
    eval_f = objective_fn,
    lb = lb,  
    ub = ub,
    opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 1000)
  )
  
  return(result)
  
}


library(plotly)


plot_empirical_loglikelihood <- function(x_data, alpha_values, beta_values){
  
  alpha_values <- alpha_values
  beta_values <- beta_values  
  log_likelihood_grid <- outer(
    alpha_values, beta_values,
    Vectorize(function(a, b) log_likelihood_empirical(c(a, b), x_data))
  )
  
  plot_ly(
    x = ~alpha_values, 
    y = ~beta_values, 
    z = ~log_likelihood_grid, 
    type = "surface"
  ) %>%
    layout(
      title = "Log-Likelihood Surface",
      scene = list(
        xaxis = list(title = "Alpha"),
        yaxis = list(title = "Beta"),
        zaxis = list(title = "Log-Likelihood")
      )
    )
  
}

# UTILIZAR O MODELO. ESTIMAR OS HIPER-PARAMETROS ALFA E BETA, E USÁ-LOS NA DISTRIBUIÇÃO GAMA DO MODELO JAGS
# set.seed(42)
# estimates_empirical <- estimates_empirical_prior(ramarendra(50,4), c(2,2), c(0,0), c(1000,1000))
# cat("Alpha:", estimates_empirical$solution[1], "\nBeta:", estimates_empirical$solution[2], "\n")



MCMC_empirical <- function(N_values, parameters, sims){
  
  cat("model{
	for(i in 1:N){


		# Likelihood

		L[i] <- ((theta^4)/(theta^3 + theta^2 + 2*theta + 6))*(1+y[i]+y[i]^2+y[i]^3)*exp(-theta*y[i])

		pro[i] <- L[i]/C
		ones[i] ~ dbern(pro[i])
	}
	
    theta ~ dgamma(2.12227, 0.3923691)  

    C <- 1000000
}",
      file="gammaempirical.jags")
  
  
  
  initialization <- function(){
    list("theta"=3)
  }
  
  df <- data.frame(theta = numeric(), N = numeric(), media = numeric(), 
                   vies = numeric(), eqm = numeric(), pc = numeric(),
                   stringsAsFactors = FALSE)
  
  S <- sims
  for (theta0 in parameters) {
    for (N in N_values) {
      theta_c=numeric()
      
      theta_inf=numeric()
      theta_sup=numeric()
      theta_vies=numeric()
      theta_eqm=numeric()
      
      for(j in 1:S){
        y <- ramarendra(N, theta0)
        
        ones <- rep(1, N)
        a.dat <- list("N", "y", "ones")
        
        parameters <- c("theta")
        
        model.fit <- jags(data = a.dat,
                          inits = initialization,
                          parameters.to.save = parameters,
                          n.chains = 1,
                          n.iter = 20000,
                          n.burnin = 5000,
                          n.thin = 10,
                          model.file ='gammaempirical.jags')
        
        model.fit.mcmc <- as.mcmc(model.fit)
        thetac=as.matrix(model.fit.mcmc[,2])   
        
        
        theta_c[j]=mean(thetac)
        theta_inf[j]=quantile(thetac, 0.025)
        theta_sup[j]=quantile(thetac, 0.975)
        theta_vies[j]=theta_c[j]-theta0
        theta_eqm[j]=(theta_c[j]-theta0)^2
      }
      
      pc_theta=numeric()
      for(i in 1:S){
        if((theta_inf[i]<theta0) & (theta0<theta_sup[i])){
          pc_theta[i]=1
        } else {
          pc_theta[i]=0
        }
      }
      
      pc=mean(pc_theta)
      
      media=mean(theta_c)
      eqm=mean(theta_eqm)
      vies=mean(theta_vies)
      
      new_row <- c(theta0, N, media, vies, eqm, pc)
      df[nrow(df) + 1,] = new_row
    }
  }
  
  return(df)
  
  
}
