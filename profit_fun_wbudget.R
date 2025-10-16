# trial parameters
alpha_trt <- 36   # distribution parameter of theta, generated from phase 2
alpha_c <- 44    # distribution parameter of theta, generated from phase 2
beta_trt <- 547  # distribution parameter of theta, generated from phase 2
beta_c <- 438    # distribution parameter of theta, generated from phase 2
follow <- 15  # follow up time

# specify the total revenue
rv <- 1000

# specify the costs
c1 <- 2  # cost per person in trt arm
c2 <- 2  # cost per person in control arm
c3 <- 0.1  # cost per person per followup time related to study duration

# specify number of iteration for power
n_iter <- 100

# specify sample size
n <- 200

# define profit function, with budget
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