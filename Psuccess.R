library(MASS)
library(stats)
library(simsurv)
library(dplyr)

# Set parameters for the inverse gamma distributions
alpha_trt <- 36
beta_trt <- 547
alpha_c <- 44
beta_c <- 438
follow <- 15    # define followup time

#######################################################
# Simulate samples from inverse gamma distributions
num_samples <- 100000
theta_trt_samples <- nimble::rinvgamma(num_samples, shape = alpha_trt, scale = beta_trt)
theta_c_samples <- nimble::rinvgamma(num_samples, shape = alpha_c, scale = beta_c)

# Compute differences
delta_samples <- theta_trt_samples - theta_c_samples

# KDE for PDF of delta
kde <- density(delta_samples)


# pi function (density function)
pi_delta <- approxfun(kde$x, kde$y, rule = 2)
########################################################



# Define the clinically significant difference
d0 <- 0.5  # Example value
z_alpha <- qnorm(0.975)  # 95% confidence level


# function parameters
r <- 0.5
n <- 800

# simulate data
n1 <- round(n*r)
n2 <- n-n1

######################### Functions define outside P_success  ########################
# function for computing single power
compute_power <- function(d, sigma2_trt, sigma2_c, d0, z_alpha) {
  1 - pnorm((d0 - d) / sqrt(sigma2_trt + sigma2_c) + z_alpha)
  
}


# define the integrand
integrand <- function(d) {
  
  # simulate data based on mean difference
  theta_trt = nimble::rinvgamma(1, shape = alpha_trt, scale = beta_trt)
  theta_c = theta_trt - d
  
  # length(theta_c) != 1 ||
  # Prevent negative or zero rate parameters
  if ( theta_c < 0) {
    return(0)
  }
  
  else if (theta_c >= 0) {
    
    
    cov_trt = data.frame(id = 1:n1, trt = rep(1,length(1:n1)))
    t_trt = simsurv(dist = c( "exponential"), lambdas =  1/theta_trt, x = cov_trt, maxt = follow)
    
    cov_c = data.frame(id = (n1+1):(n1+n2), trt = rep(0,length((n1+1):(n1+n2))))
    t_c = simsurv(dist = c( "exponential"), lambdas =  1/theta_c, x = cov_c, maxt = follow) %>%
      mutate(id=c((n1+1):(n1+n2)))
    
    # Define the variances
    sigma2_theta_trt <- (beta_trt+sum(t_trt$eventtime))^2/((alpha_trt+sum(t_trt$status)-1)^2*(alpha_trt+sum(t_trt$status)-2))
    
    sigma2_theta_c <- (beta_c+sum(t_c$eventtime))^2/((alpha_c+sum(t_c$status)-1)^2*(alpha_c+sum(t_c$status)-2))
    
    
    
    power <- compute_power(d, sigma2_theta_trt, sigma2_theta_c, d0, z_alpha)
    
    
    power * pi_delta(d)
  }
}
##############################################################################################

#P_success <- stats::integrate(integrand, lower = min(kde$x), upper = max(kde$x))$value



# approximation using summation
step_size <- 0.1
diff <- seq(min(kde$x),max(kde$x),by=step_size)
out <- rep(NA,length(diff))
for (i in 1:length(diff)){
  cat(i,"\n")
  d <- diff[i]
  out[i] <- integrand(d)
}
P_success <- sum(out)*step_size

#############################################################################

P_success <- function(n,r,step_size,
                      alpha_trt,
                      beta_trt,
                      alpha_c,
                      beta_c,
                      follow, 
                      num_samples,
                      d0,
                      z_alpha){
  #######################################################
  # Simulate samples from inverse gamma distributions
  theta_trt_samples <- nimble::rinvgamma(num_samples, shape = alpha_trt, scale = beta_trt)
  theta_c_samples <- nimble::rinvgamma(num_samples, shape = alpha_c, scale = beta_c)
  
  # Compute differences
  delta_samples <- theta_trt_samples - theta_c_samples
  
  # KDE for PDF of delta
  kde <- density(delta_samples)
  
  
  # pi function (density function)
  pi_delta <- approxfun(kde$x, kde$y, rule = 2)
  ###########################################################
  
  ######################### Functions define for P_success  ########################
  # function for computing single power
  compute_power <- function(d, sigma2_trt, sigma2_c, d0, z_alpha) {
    1 - pnorm((d0 - d) / sqrt(sigma2_trt + sigma2_c) + z_alpha)
    
  }
  
  
  # define the integrand
  integrand <- function(d) {
    
    # simulate data based on mean difference
    theta_trt = nimble::rinvgamma(1, shape = alpha_trt, scale = beta_trt)
    theta_c = theta_trt - d
    
    # length(theta_c) != 1 ||
    # Prevent negative or zero rate parameters
    if ( theta_c < 0) {
      return(0)
    }
    
    else if (theta_c >= 0) {
      
      
      cov_trt = data.frame(id = 1:n1, trt = rep(1,length(1:n1)))
      t_trt = simsurv(dist = c( "exponential"), lambdas =  1/theta_trt, x = cov_trt, maxt = follow)
      t_trt <- t_trt %>% mutate_at(c("eventtime"),~(scale(.) %>% as.vector))
      
      cov_c = data.frame(id = (n1+1):(n1+n2), trt = rep(0,length((n1+1):(n1+n2))))
      t_c = simsurv(dist = c( "exponential"), lambdas =  1/theta_c, x = cov_c, maxt = follow) %>%
        mutate(id=c((n1+1):(n1+n2)))
      t_c <- t_c %>% mutate_at(c("eventtime"),~(scale(.) %>% as.vector))
      
      
      # Define the variances
      sigma2_theta_trt <- (beta_trt+sum(t_trt$eventtime))^2/((alpha_trt+sum(t_trt$status)-1)^2*(alpha_trt+sum(t_trt$status)-2))
      
      sigma2_theta_c <- (beta_c+sum(t_c$eventtime))^2/((alpha_c+sum(t_c$status)-1)^2*(alpha_c+sum(t_c$status)-2))
      
      
      
      power <- compute_power(d, sigma2_theta_trt, sigma2_theta_c, d0, z_alpha)
      
      
      power * pi_delta(d)
    }
  }
  ##############################################################################################
  
  
  ########################################################
  # simulate data
  n1 <- round(n*r)
  n2 <- n-n1
  
  # approximation using summation
  diff <- seq(min(kde$x),max(kde$x),by=step_size)
  out <- rep(NA,length(diff))
  for (i in 1:length(diff)){
    set.seed(66)
    d <- diff[i]  
    out[i] <- integrand(d)
  }
  P_success <- sum(out)*step_size
  
  return(P_success)
}