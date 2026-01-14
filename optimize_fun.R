rm(list = ls())

library(MASS)
library(stats)
library(simsurv)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(xtable)


## The P_success function calculates the probability of success based on the historical priors
## n: sample size
## r: randomization ratio
## alpha_trt, beta_trt, alpha_c, beta_c: IG distribution parameters for prior distributions
## follow: follow up time
## num_samples: number of samples to compute the density of d, difference between treatment and control
## d0: clinically significant difference in hypothesis
## z_alpha: z value of the alpha level


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

## The profit function calculates the expected profit if the trial succeeds based on the revenue and cost 
profit <- function(r,n,step_size,d0,num_samples,z_alpha,rv,budget,c1,c2,c3,alpha_trt,
                   alpha_c,beta_trt,beta_c,follow,p0){
  
  cat("r:",r,"\n")
  
  n1 <- round(n*r)
  n2 <- n-n1
  
  follow_cost <- c3*follow*n
  
  # define the profit if trial success
  v1 <- rv - (c1*n1+c2*n2+follow_cost)
  
  # define the cost if trial fails
  v2 <- c1*n1 + c2*n2 + follow_cost
  
  cat("cost:",v2,"\n")
  
  # round the ratio
  #digits <- length(as.numeric(strsplit(as.character(n),"")[[1]])) + 1
  #r <- round(r, digits = 3)
  
  
  # compute the power
  prob_success <- P_success (n,r,step_size,
                             alpha_trt,
                             beta_trt,
                             alpha_c,
                             beta_c,
                             follow, 
                             num_samples,
                             d0,
                             z_alpha)
  
  cat("prob_success:",prob_success,"\n")
  
  
  ###################
  
  if (prob_success >= p0 && v2 <= budget){  #  only continue if the probability of success is greater than 0.6
    
    expected_profit <- v1 * prob_success - v2 * (1-prob_success)}
  
  else{
    expected_profit <- -10000000000}
  
  #############
  
  cat("profit:", expected_profit,"\n")
  
  return(expected_profit)
  
}



# trial parameters
alpha_trt <- 36   # distribution parameter of theta, generated from phase 2
alpha_c <- 44    # distribution parameter of theta, generated from phase 2
beta_trt <- 547  # distribution parameter of theta, generated from phase 2
beta_c <- 438    # distribution parameter of theta, generated from phase 2
follow <- 15  # follow up time

# specify the total revenue
rv <- 20000

# specify the costs
c1 <- 2  # cost per person in trt arm
c2 <- 4  # cost per person in control arm
c3 <- 0.1 # cost per person per followup time related to study duration

# parameters needed for P_success
num_samples <- 100000
d0 <- 0.5  # Example value
z_alpha <- qnorm(0.975)  # 95% confidence level
step_size <- 0.5


# specify expected probability of success
p0 <- 0.7

# specify sample size list
n_list <- seq(200, 400, by=50)

# specify budget list
budget <- 200000

# specify minimum power?

# pre-define another list of r
r1_list <- rep(NA,length(n_list))


#################################

# specify start value for randomization ratio r
start_values <- 0.6

# optimizing the profit based on 

for (i in 1:length(n_list)){
  n <- n_list[i]
  cat("n:",n,"\n")
  out <- optim(par=start_values ,
               fn= profit, 
               method = "Brent", 
               lower = 0.1 , 
               upper= 0.9,
               control=list(fnscale=-1),
               step_size= step_size,
               d0=d0,
               num_samples=num_samples,
               z_alpha=z_alpha,
               rv= rv,
               budget=budget,
               c1= c1,
               c2=c2,
               c3=c3,
               n=n,
               alpha_trt=alpha_trt,
               alpha_c=alpha_c,
               beta_trt=beta_trt,
               beta_c=beta_c,
               follow=follow,
               p0=p0)
  r1_list[i] <- out$par
}



# summarise in table
table <- data.frame(Sample_size = n_list, Ratio = r1_list)
new_table <- round(table,2)
print(xtable(new_table, type = "latex"), file = "sample_size_fixed_budget.tex")

#define scatterplot
my_plot <- ggplot(table, aes(x=Sample_size, y=Ratio)) +
  geom_point()+ylim(0.4, 0.9)

jpeg('sample_size_fixed_budget.jpg')
my_plot
dev.off()
#define table
my_table <- tableGrob(table)

#create scatterplot and add table underneath it
grid.arrange(my_plot,my_table)

save(my_plot, table, my_table, file = "sample_size_fixed_budget.rdata")
