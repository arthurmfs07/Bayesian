library(AdequacyModel)


# Best model with lowest, W, A, KS, -loglike
 
# 1. Amarendra (AMAR) - "famarendra"
 
# 2. ExponentiatedHL(EHL) - "fghl"
 
# 3. GeneralizedHL(GHL) - "fehl"
 
# 4. Lindley(Li) - "flind"
 
# 5. InverseLindley(ILi) - "filind"
 
# 6. Exponential(Exp) - "fexp"




#########Amarendra#####

pdfamarendra <- function(par, x) 
  {
  theta <- par[1]
  fx = theta^4 / (theta^3 + theta^2 + 2 * theta + 6) * (1 + x + x^2 + x^3) * exp(-theta * x)
  return(fx)
}


cdfamarendra <- function(par, x) 
  {
  theta <- par[1]
  Fx = 1-(1+((x^3*theta^3+theta^2*(theta+3)*x^2+theta*(theta^2+2*theta+6)*x)/(theta^3+theta^2+2*theta+6)))*exp(-theta*x)
  return(Fx)
}


########EHL#####
cdfehl=function(par,x)
{
  lambda=par[1]
  g=2*exp(-x)/(1+exp(-x))^2
  G=(1-exp(-x))/(1+exp(-x))
  f=G^lambda
  return(f)
}

pdfehl=function(par,x)
{
  lambda=par[1]
  g=2*exp(-x)/(1+exp(-x))^2
  G=(1-exp(-x))/(1+exp(-x))
  f=lambda*g*G^(lambda-1)
  return(f)
}


#############GHL############
cdfghl=function(par,x)
{
  lambda=par[1]
  f=1-((2*exp(-x))/(1+exp(-x)))^lambda
  return(f)
}

pdfghl=function(par,x)
{
  lambda=par[1]
  f=2^lambda*lambda*exp(x)*(1/(exp(x)+1))^(lambda+1)
  return(f)
}


#############Lindley############
cdflind=function(par,x)
{
  lambda=par[1]
  f=1-(1+lambda*x/(1+lambda))*exp(-lambda*x)
  return(f)
}

pdflind=function(par,x)
{
  lambda=par[1]
  
  f=lambda^2/(lambda+1)*(1+x)*exp(-lambda*x)
  return(f)
}


#############Inverselindley############
cdfilind=function(par,x)
{
  theta=par[1]
  #if (is.infinite(x)) return(1)
  f=(1+theta/((1+theta)*x))*exp(-theta/x)
  return(f)
}

pdfilind=function(par,x)
{
  theta=par[1]

  f=theta^2/(1+theta)*((1+x)/x^3)*exp(-theta/x)
  return(f)
}


#############Exp############
cdfexp=function(par,x)
{
  lambda=par[1]
  f=1-exp(-lambda*x)
  return(f)
}

pdfexp=function(par,x)
{
  lambda=par[1]
  
  f=lambda*exp(-lambda*x)
  return(f)
}



#----------------------------Initialize adequacy--------------------------- #

initialize_fitting <- function(data_test, starts){

  famarendra=goodness.fit(pdf=pdfamarendra, cdf=cdfamarendra,
                          starts = c(starts[1]), data = data_test,
                          method="N", domain=c(0.01,1e3))

  fghl=goodness.fit(pdf=pdfghl, cdf=cdfghl,
                    starts = c(starts[2]), data = data_test,
                    method="N", domain=c(0.01,1e3))

  # fehl=goodness.fit(pdf=pdfehl, cdf=cdfehl,
  #                    starts = c(starts[3]), data = data_test,
  #                    method="N", domain=c(0.01,1e3))


  flind=goodness.fit(pdf=pdflind, cdf=cdflind,
                     starts = c(starts[4]), data = data_test,
                     method="N", domain=c(0.01,1e3))

  filind=goodness.fit(pdf=pdfilind, cdf=cdfilind,
                      starts = c(starts[5]), data = data_test,
                      method="N", domain=c(0.01,1e3))

  fexp=goodness.fit(pdf=pdfexp, cdf=cdfexp,
                    starts = c(starts[6]), data = data_test,
                    method="N", domain=c(0.01,1e3))

  list(famarendra = famarendra, fghl = fghl, flind = flind, 
       filind = filind, fexp = fexp)
    
  
#   list(famarendra = famarendra, fghl = fghl, fehl = fehl,
#         flind = flind, filind = filind, fexp = fexp)
}

#----------------------------Plot Adequacy--------------------------- #

plots_adequacy <- function(data_test, fits, ylim) {

  max_val <- max(data_test)+10
  min_val <- 0
  hist(data_test, prob = TRUE, main="Data Histogram", xlab="x",
       ylab="f(x)", ylim = ylim, lwd=1, col="white")

  curve(pdfamarendra(fits$famarendra$mle, x), min_val, max_val, add = TRUE, col = "1", lwd = 2, lty = 1)
  curve(pdfghl(fits$fghl$mle, x), min_val, max_val, add = TRUE, col = "2", lwd = 2, lty = 1)
  # curve(pdfehl(fits$fehl$mle, x), min_val, max_val, add = TRUE, col = "3", lwd = 2, lty = 1)
  curve(pdflind(fits$flind$mle, x), min_val, max_val, add = TRUE, col = "1", lwd = 2, lty = 2)
  curve(pdfilind(fits$filind$mle, x), min_val, max_val, add = TRUE, col = "2", lwd = 2, lty = 2)
  curve(pdfexp(fits$fexp$mle, x), min_val, max_val, add = TRUE, col = "3", lwd = 2, lty = 2)

 legend(x = "topright", y = "topright", legend = c("AMAR", "GHL", "Li", "Invli", "Exp"),
        cex = 1, lty = c(1, 1, 2, 2, 2), lwd = 2, col = c(1, 2, 1, 2, 3))
 
 # legend(x = "topright", y = "topright", legend = c("AMAR", "GHL", "EHL", "Li", "Invli", "Exp"),
 #        cex = 1, lty = c(1, 1, 1, 2, 2, 2), lwd = 2, col = c(1, 2, 3, 1, 2, 3))
 
 }


#----------------------------Plot table -------------------------- #

table_adequacy <- function(){
  results_df <- data.frame(
    Distributions = c("Amarendra", "GHL", "Li", "InvLi", "Exp"),
    `MLE` = c(
      fit$famarendra$mle,
      fit$fghl$mle,
      fit$flind$mle,
      fit$filind$mle,
      fit$fexp$mle
    ),
    SEs = c(
      fit$famarendra$Erro,
      fit$fghl$Erro,
      fit$flind$Erro,
      fit$filind$Erro,
      fit$fexp$Erro
    ),
    `neg_loglike` = c(
      fit$famarendra$Value,
      fit$fghl$Value,
      fit$flind$Value,
      fit$filind$Value,
      fit$fexp$Value
    ),
    `A*` = c(
      fit$famarendra$A,
      fit$fghl$A,
      fit$flind$A,
      fit$filind$A,
      fit$fexp$A
    ),
    `W*` = c(
      fit$famarendra$W,
      fit$fghl$W,
      fit$flind$W,
      fit$filind$W,
      fit$fexp$W
    ),
    KS = c(
      fit$famarendra$KS$statistic,
      fit$fghl$KS$statistic,
      fit$flind$KS$statistic,
      fit$filind$KS$statistic,
      fit$fexp$KS$statistic
    ),
    `p-value` = c(
      fit$famarendra$KS$p.value,
      fit$fghl$KS$p.value,
      fit$flind$KS$p.value,
      fit$filind$KS$p.value,
      fit$fexp$KS$p.value
    )
  )
  
  print(results_df)
}



table_info_criteria <- function(){
  results_df <- data.frame(
    Distributions = c("Amarendra", "GHL", "Li", "InvLi", "Exp"),
    `MLE` = c(
      fit$famarendra$mle,
      fit$fghl$mle,
      fit$flind$mle,
      fit$filind$mle,
      fit$fexp$mle
    ),
    `AIC` = c(
      fit$famarendra$AIC,
      fit$fghl$AIC,
      fit$flind$AIC,
      fit$filind$AIC,
      fit$fexp$AIC
    ),
    `BIC` = c(
      fit$famarendra$BIC,
      fit$fghl$BIC,
      fit$flind$BIC,
      fit$filind$BIC,
      fit$fexp$BIC
    ),
    `CAIC` = c(
      fit$famarendra$CAIC,
      fit$fghl$CAIC,
      fit$flind$CAIC,
      fit$filind$CAIC,
      fit$fexp$CAIC
    ),
    `HQIC` = c(
      fit$famarendra$HQIC,
      fit$fghl$HQIC,
      fit$flind$HQIC,
      fit$filind$HQIC,
      fit$fexp$HQIC
    )
    
  )
  
  print(results_df)
}


