# Set work directory, import data, and include relevant libraries
setwd("C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Codes")
report_path = "C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Report/images"
load("claims.RData")

library(stats4)
library(ggplot2)
library(dplyr)
library(actuar)


load("mle_summary_pareto.RData")

# Check initial conditions
mle_summary_pareto$dl_dalpha <- rep(17, nrow(mle_summary_pareto))
mle_summary_pareto$dl_dlambda <- rep(17, nrow(mle_summary_pareto))

for (row in 1:nrow(mle_summary_pareto)){
  alpha <- mle_summary_pareto$mle_alpha[row]
  lambda <- mle_summary_pareto$mle_lambda[row]
  excess <- claims$excess
  loss <- claims$loss
  mle_summary_pareto$dl_dalpha[row] <- n/alpha + 
    sum(log((lambda + excess)/(lambda + loss)))
  mle_summary_pareto$dl_dlambda[row] <- sum(alpha/(lambda + excess) - 
                                              (alpha + 1)/(lambda + loss))
}

# Function to calculate MLE alpha given a lambda
f_alpha <- function(lambda){
  return(-n/sum(log((lambda + claims$excess)/(lambda + claims$loss))))
}

## mean excess function
u_pareto <- function(deductible, alpha, lambda) {
  constant <- 1/(1-ppareto(deductible, alpha, lambda))
  integral <- deductible*(1-ppareto(deductible, alpha, lambda)) - 
    (lambda^alpha)/(1-alpha)*(deductible+lambda)^(-alpha + 1)
  print((lambda^alpha))
  return(constant * integral)
}

## Summarise data by excess level
claims_XS_mean <- claims %>%
  group_by(excess) %>%
  summarise(mean = mean(loss), 
            policies = n()) %>%
  mutate(proportion = policies/sum(policies)) %>%
  arrange(desc(policies))

## Overall theoretical average claim size function
u_pareto_aggregate <- function(alpha, lambda){
  overall_avg <- 0
  for (i in 1:nrow(claims_XS_mean)){
    excess <- claims_XS_mean$excess[i]
    weight <- claims_XS_mean$proportion[i]
    theo_average <- u_pareto(excess, alpha, lambda)
    overall_avg <- overall_avg + weight * theo_average
  }
  return(overall_avg)
}


claims_XS_mean <- claims_XS_mean %>%
  mutate(mle_alpha = 0,
         mle_lambda = 0)


# for(row in 1:3){
#   excess <- claims_XS_mean$excess[row]
#   mean <- claims_XS_mean$mean[row]
#   f <- function(lambda){
#     return(u_pareto(excess, f_alpha(lambda), lambda) - mean)
#   }
#   output <- uniroot(f, interval = c(10000, 20000))
#   lambda <- output$root
#   claims_XS_mean$mle_alpha[row] <- f_alpha(lambda)
#   claims_XS_mean$mle_lambda[row] <- lambda
# }

#claims_XS_mean$excess <- claims_XS_mean$excess %>% factor()
# Plot MME/MLE estimates
pareto_MLE_graph <- ggplot() + 
  geom_point(data = claims_XS_mean, aes(x = mle_alpha, y = mle_lambda, 
                                        col = excess, size = proportion)) + 
  labs(x = expression(widehat(alpha)), y = expression(widehat(lambda)),
       size = "Proportion of policies", col = "Deductible level ($)") +
  ggtitle("Pareto parameter estimates across different deductible levels")

# Find final estimates, weighted average
alpha_pareto <- sum(claims_XS_mean$mle_alpha * claims_XS_mean$proportion)
print(paste("The estimate for alpha is", alpha_pareto))
lambda_pareto <- sum(claims_XS_mean$mle_lambda * claims_XS_mean$proportion)
print(paste("The estimate for lambda is", lambda_pareto))

pareto_MLE_graph1 <- pareto_MLE_graph + geom_point(aes(x = alpha_pareto, y = lambda_pareto), 
                              shape = 4, size = 10, col = "blue", stroke = 2) +
  geom_label(aes(x = alpha_pareto, y = lambda_pareto - 100,
                 label ="Final estimate"))

# png(paste(report_path, "/Pareto_estimates.png", sep = ""))
# pareto_MLE_graph1
# dev.off()

# Find the (negative) log likelihood and AIC, BIC
## Create (negative) log likelihood function
f.nLL_pareto <- function(alpha, lambda){
  l <- sum(log(dpareto(claims$loss, alpha, lambda)) - 
             log(1-ppareto(claims$excess, alpha, lambda)))
  return(-l)
}
nLL_pareto <- f.nLL_pareto(alpha_pareto, lambda_pareto)
print(paste("The negative log likelihood is",nLL_pareto))

#AIC
AIC_pareto <- 2*nLL_pareto + 2 * 2
print(paste("The AIC is", AIC_pareto))
BIC_pareto <- 2*nLL_pareto + log(n)*2
print(paste("The BIC is", BIC_pareto))


# 
for (row in 1:nrow(mle_summary_pareto)){
  mle_summary_pareto$theo_mean[row] <- u_pareto(700, mle_summary_pareto$mle_alpha[row], mle_summary_pareto$mle_lambda[row])
}