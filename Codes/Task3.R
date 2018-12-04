library(actuar)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Clear everything to prevent unintentional reuse of objects
rm(list=ls())

# Toggle to export plots
report_path <- 
  "C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Report/images/"
export <- TRUE

# Points to evaluate CDF at
x <- 10*(1:5)

# Gamma distribution parameter
alpha <- 5
beta <- 0.8

# Binomial distribution parameters
n <- 5
p <- 0.997
q <- 1 - p

# Mass dispersal approximation parameters
h <- 1
m <- 70

################# Part A ################
round(pgamma(x, shape = n*alpha, rate = beta),7)

# Panjer values a and b
a <- -p/(1-p)
b <- (n+1)*p/(1-p)

# Discretisation of gamma
gamma_discrete <- discretise(pgamma(x, shape = alpha, rate = beta), 
                             from = 0, 
                             to = h*m, 
                             step = h,
                             method = "rounding")

# Function to output value of discretisation
f_gamma_discrete <- function(x){
    mass <- gamma_discrete[x + 1]
  return(mass)
}

# Initial values
g0 <- f_gamma_discrete(0)
f0 <- (q + p * g0)^n

f_y1 <- aggregateDist(method = "recursive", model.freq = "binomial", 
                      model.sev = gamma_discrete,
                      size = 5, prob = 0.99, convolve = 0, x.scale = 1)

# Run through panjer recursion
f_y <- f0
for (s in 1:max(x)){
  j <- 1:s
  f_y[s + 1] <- 1/(1 - a * g0) * sum((a + b * j / s) *
                                 f_gamma_discrete(j) * f_y[s - j + 1])
}

# Store values as dataframe
panjer_df <- data.frame(y = 0:(length(f_y)-1), f = f_y)
for (row in 1:nrow(panjer_df)){
  panjer_df$F_y[row] <- sum(panjer_df$f[1:row])
}
panjer_df <- panjer_df %>% 
  mutate(True_density = dgamma(y, shape = n*alpha, rate = beta),
         True_distribution = pgamma(y, shape = n*alpha, rate = beta))

panjer_plot_f <- ggplot(data = panjer_df, aes(y, f_y)) + 
  geom_bar(stat = 'identity') +
  labs(title = "Panjer's recursion mass function", y = expression(f[Y](y))) +
  theme_bw() + geom_line(aes(x = y, y = True_density), size = 0.7, color = 'red',
                         alpha= 0.7)
panjer_plot_f

interleaved_Fy <- rep(0, 2*nrow(panjer_df))
interleaved_Fy[seq(from = 1, by = 2, length.out = nrow(panjer_df))] <- panjer_df$F_y
interleaved_Fy[seq(from = 2, by = 2, length.out = nrow(panjer_df))] <- panjer_df$F_y

interleaved_y <- interleaved_Fy
interleaved_y[seq(from = 1, by = 2, length.out = nrow(panjer_df))] <- panjer_df$y
interleaved_y[seq(from = 2, by = 2, length.out = nrow(panjer_df))] <- 
  c(panjer_df$y[-1], tail(panjer_df$y, 1) + 0.1)

panjer_df_F <- data.frame(y = interleaved_y, F_y = interleaved_Fy)

panjer_plot_F <- ggplot(data = panjer_df_F, aes(y, F_y))  + geom_line() +
  labs(title = "Panjer's recursion distribution function", 
       y = expression(F[Y](y))) +
  theme_bw() +
  geom_line(data = panjer_df, aes(x = y, y = True_distribution), color = 'red', 
            size = 0.7, alpha = 0.8)
panjer_plot_F

panjer_combined <- grid.arrange(panjer_plot_f, panjer_plot_F, ncol = 2)

panjer_comparison <- panjer_df %>%
  filter(y %in% x) %>%
  select(-f, -True_density) %>%
  mutate(Difference = F_y - True_distribution) %>%
  mutate(perc_diff = round(100*Difference/True_distribution, 3)) %>% 
  round(6)

# Export graph
if (export == TRUE){
  ggsave(paste(report_path, "Panjer_combined.pdf", sep = ""), panjer_combined,
         width = 10, height = 5)
}