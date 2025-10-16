rm(list = ls())

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(xtable)

source("Psuccess.R")
source("profit_fun_wbudget.R")

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
n_list <- seq(200, 2000, by=50)

# specify budget list
budget <- 200000

# specify minimum power?

# pre-define another list of r
r1_list <- rep(NA,length(n_list))

#################################
# single simulation check convergence
r_ite <- c()
prob_ite <- c()
cost_trt_ite <- c()
cost_c_ite <- c()
profit_ite <- c()
#################################

# specify start value
start_values <- 0.6

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
