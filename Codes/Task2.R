library(dplyr)
library(actuar)
library(ggplot2)
library(reshape2)

# Clear everything to prevent unintentional reuse of objects
rm(list=ls())

# Export parameters
report_path <-
  "C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Report/images/"
export <- TRUE

# Loss distribution parameters
alpha = 20
beta = 0.2
v_c0 = c(0, 10, 20, 50)
premium = 1.2 * alpha / beta

# Initialise dataframe to store approximations
all_psi2 <- data.frame(c0 = v_c0,
                       actual = rep(0, length(v_c0)),
                       h1m150 = rep(0, length(v_c0)),
                       h1m300 = rep(0, length(v_c0)),
                       h5m30 = rep(0, length(v_c0)),
                       h5m60 = rep(0, length(v_c0)))  %>%
  mutate(lower_d1 = 0,
         lower_d5 = 0,
         upper_d1 = 0,
         upper_d5 = 0)

row.names(all_psi2) <- v_c0
for (init_capital in v_c0) {
  c0 <- init_capital
  integrand <- function(y){
    (1 - pgamma(c0 + 2 * premium - y, shape = alpha, rate = beta)) *
      dgamma(y, shape = alpha, rate = beta)
  }
  integral_val <- integrate(integrand, lower = 0, upper = c0 + premium)$value
  psi_2 <- integral_val + 1 - pgamma(c0 + premium, shape = alpha, rate = beta)
  all_psi2[as.character(c0), "actual"] = psi_2
}

############################## e) ##################################
######## A: Method of rounding ####################


# Create a function for G*, the discretisation of G
G_discrete <- function(x, h, m, shape, rate){
  x_discrete <- floor(x/h)*h
  x_discrete <- min(h*m, x_discrete)
  if (x_discrete == h*m){
    G <- 1
  } else {
    G <- pgamma(x_discrete + h/2, shape = alpha, rate = beta)
  }
  return(G)
}

for (h in c(1, 5)){
  for (m in c(150, 300)/h){
    for (init_capital in v_c0){
      c0 <- init_capital
      
      g_discrete <- discretise(pgamma(x, shape = alpha, rate = beta),
                               from = 0, to = h*m,
                               step = h, method = "rounding")
      g_discrete[m + 1] <- 1- pgamma(h*m-0.5*h, shape = alpha, rate = beta)
      g_discrete[(m + 2):(m+30)] <- 0
      
      psi_1 <- function(c){
        return(1 - G_discrete(c + premium, h, m, shape = alpha, rate = beta))
      }
      
      sum_index <- 0:round(floor(c0 + premium)/h)
      psi_2 <- 0
      for (i in sum_index){
        psi_2 <- psi_2 + g_discrete[i + 1] * psi_1(c0 + premium - i*h)
      }
      psi_2 <- psi_2 + psi_1(c0)
      all_psi2[as.character(c0), paste("h", h, "m", m, sep = "")] <- psi_2
    }
    
  }
}

######## B: Method of bounds   ####################

# Use discretize function in actuar
upper_lim <- 1000
for (c0 in v_c0){
  for (method in c("lower", "upper")){
    for (d in c(1, 5)){
      discretised_gamma <- discretise(pgamma(x, shape = alpha, rate = beta),
                                      from = 0, to = upper_lim, 
                                      step = d,
                                      method = method)
      discrete_gamma_pmf <- function(x){
        if (((x/d)%%1) != 0){
          mass <- 0
        } else if ((x < 0)|(x/d + 1 > length(discretised_gamma))){
          mass <- 0
        } else {
          mass <- discretised_gamma[x/d + 1]
        }
        return(mass)
      }
      
      discrete_gamma_cdf <- function(x){
        index <- floor(x/d + 1)
        if (index > length(discretised_gamma)){
          cdf <- 1
        } else if (index < 1){
          cdf <- 0
        } else{
          cdf <- sum(discretised_gamma[1:index])
        }
        return(cdf)
      }
      
      psi_1 <- function(c){
        1 - discrete_gamma_cdf(c + premium)
      }
      
      sum_index <- 0:floor(c0 + premium)
      psi_2 <- 0
      for (i in sum_index){
        psi_2 <- psi_2 + discrete_gamma_pmf(i) * psi_1(c0 + premium - i)
      }
      psi_2 <- psi_2 + psi_1(c0)
      all_psi2[as.character(c0), paste(method, "_d", d, sep ="")] <- psi_2
    }
  }
}

# Plot differences
all_psi2_long <- melt(all_psi2, id = c("c0"), 
                      variable.name = "Method", 
                      value.name = "psi_2")

plot_approx_methods <- ggplot() + geom_line(data = all_psi2_long %>% 
                                              filter(Method != "actual"), 
                                            aes(x = c0, y = psi_2, col = Method)) +
  geom_line(data = all_psi2_long %>% filter(Method == "actual"),
            aes(x = c0, y = psi_2), linetype = "dotted",
            size = 1.1) + theme_bw() +
  labs(title = "Comparison of approximation methods",
       x = expression(c[0]),
       y = expression(psi[2](c[0])))

if (export == TRUE){
  ggsave(paste(report_path, "Approximations.pdf", sep = ""), plot_approx_methods)
}