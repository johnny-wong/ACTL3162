library(stats4)
library(ggplot2)
library(dplyr)
library(actuar)
library(reshape2)
library(profvis)
library(microbenchmark)
profvis({
# Clear everything to prevent unintentional reuse of objects
rm(list=ls())

# Set work directory, import data, and include relevant libraries
setwd("C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Codes")
report_path <-
  "C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Report/images/"
load("claims.RData")
output <- FALSE


####### Summarise by excess level #######
claims_XS_mean <- claims %>%
  group_by(excess) %>%
  summarise(mean = mean(loss),
            policies = n()) %>%
  mutate(proportion = policies/sum(policies)) %>%
  arrange(desc(policies))

###### PARETO #############

# Function to calculate MLE alpha given a lambda
f_alpha <- function(lambda){
  return(-n/sum(log((lambda + claims$excess)/(lambda + claims$loss))))
}

## Create (negative) log likelihood function
f.nLL_pareto <- function(alpha, lambda){
  l <- sum(log(dpareto(claims$loss, alpha, lambda)) -
             log(1-ppareto(claims$excess, alpha, lambda)))
  return(-l)
}

init_lambda_pareto <- 200000

mle_pareto <- mle(f.nLL_pareto,
                  start = list(alpha = f_alpha(init_lambda_pareto),
                               lambda = init_lambda_pareto))

####### LOG NORMAL ########
# Function to calculuate the negative log likelihood
f.nLL_LN <- function(mu, sigma){
  return(-sum(log(dlnorm(claims$loss, mu, sigma)) -
                log(1-plnorm(claims$excess, mu, sigma))))
}

mle_LN <- mle(f.nLL_LN, start = list(mu = 9, sigma =0.2))
mle_LN

###### GAMMA #############
# function to calculate the negative log likelihood
f.nLL_gamma <- function(alpha, lambda){
  return(-sum(log(dgamma(claims$loss, shape = alpha, rate = lambda)) -
         log(1- pgamma(claims$excess, shape = alpha, rate = lambda))))
}

mean_loss <- mean(claims$loss)
alpha_MME <- n*mean_loss^2/sum((claims$loss - mean_loss)^2)
lambda_MME <- alpha_MME/mean_loss

mle_gamma <- mle(f.nLL_gamma, start = list(alpha = alpha_MME,
                                           lambda = lambda_MME),
      control=list(ndeps=c(1e-6, 1e-6)))


######################### Plotting #########################
############################################################

## Plotting PDF
f_pareto <- function(x, alpha, lambda){
  i <- 1:nrow(claims_XS_mean)
  prop <- claims_XS_mean$proportion
  deductible <- claims_XS_mean$excess
  result <- sum(prop[i]*dpareto(x, shape = alpha, scale = lambda)*
                  ifelse(x>=deductible, 1, 0)/
                  (1-ppareto(deductible[i], shape = alpha, scale = lambda)))
  return(result)
}

f_LN <- function(x, mu, sigma){
  i <- 1:nrow(claims_XS_mean)
  prop <- claims_XS_mean$proportion
  deductible <- claims_XS_mean$excess
  result <- sum(prop[i]*dlnorm(x, mu, sigma)*
                  ifelse(x>=deductible, 1, 0)/
                  (1-plnorm(deductible[i], mu, sigma)))
  return(result)
}

f_gamma <- function(x, alpha, lambda){
  i <- 1:nrow(claims_XS_mean)
  prop <- claims_XS_mean$proportion
  deductible <- claims_XS_mean$excess
  result <- sum(prop[i]*dgamma(x, shape = alpha, rate = lambda)*
        ifelse(x>=deductible, 1, 0)/
        (1-pgamma(deductible[i], shape = alpha, rate = lambda)))
  return(result)
}

x <- claims_XS_mean$excess[1]:ceiling(max(claims$loss))
y_pareto <- x
y_gamma <- x
y_LN <- x
for (ii in 1:length(x)){
  y_pareto[ii] <- f_pareto(x[ii],
                           coef(mle_pareto)[1],
                           coef(mle_pareto)[2])
  y_gamma[ii] <- f_gamma(x[ii],
                         coef(mle_gamma)[1],
                         coef(mle_gamma)[2])
  y_LN[ii] <- f_LN(x[ii],
                   coef(mle_LN)[1],
                   coef(mle_LN)[2])
}


fitted_pdf <- data.frame(loss = x,
                         pdf_pareto = 550000*y_pareto,
                         pdf_LN = 550000*y_LN,
                         pdf_gamma = 570000*y_gamma)
fitted_pdf_melted <- fitted_pdf %>%
  melt(id = c("loss"), variable.name = "Distribution")


plot_PDF_comp <- ggplot() + geom_histogram(data = claims, aes(loss)) +
  geom_line(data = fitted_pdf_melted, aes(x = loss, y = value,
                                          color=Distribution)) +
  labs(x = "Loss ($)", y = "Number of policies",
       title = "Fitted densities")+theme_bw() +
  geom_vline(xintercept = claims_XS_mean$excess, alpha = 0.5,
             linetype = "dotted") + theme_bw()


# Plotting CDF
F_pareto <- function(x, alpha, lambda){
  i <- 1:nrow(claims_XS_mean)
  prop <- claims_XS_mean$proportion
  deductible <- claims_XS_mean$excess
  result <- sum(prop * (ppareto(x, shape = alpha, scale = lambda) -
                  ppareto(deductible, shape = alpha, scale = lambda))*
                  ifelse(x>=deductible, 1, 0) /
                  (1-ppareto(deductible, shape = alpha, scale = lambda)))
  return(result)
}

F_LN <- function(x, mu, sigma){
  i <- 1:nrow(claims_XS_mean)
  prop <- claims_XS_mean$proportion
  deductible <- claims_XS_mean$excess
  result <- sum(prop*(plnorm(x, mu, sigma)-plnorm(deductible, mu, sigma))*
                  ifelse(x>=deductible, 1, 0)/
                  (1-plnorm(deductible, mu, sigma)))
  return(result)
}

F_gamma <- function(x, alpha, lambda){
  i <- 1:nrow(claims_XS_mean)
  prop <- claims_XS_mean$proportion
  deductible <- claims_XS_mean$excess
  result <- sum(prop*(pgamma(x, shape = alpha, rate = lambda) -
                        pgamma(deductible, shape = alpha, rate = lambda)) *
                  ifelse(x>=deductible, 1, 0)/
                  (1-pgamma(deductible, shape = alpha, rate = lambda)))
  return(result)
}

Y_pareto <- x
Y_gamma <- x
Y_LN <- x
for (ii in 1:length(x)){
  Y_pareto[ii] <- F_pareto(x[ii],
                           coef(mle_pareto)[1],
                           coef(mle_pareto)[2])
  Y_gamma[ii] <- F_gamma(x[ii],
                         coef(mle_gamma)[1],
                         coef(mle_gamma)[2])
  Y_LN[ii] <- F_LN(x[ii],
                   coef(mle_LN)[1],
                   coef(mle_LN)[2])

}


fitted_cdf <- data.frame(loss = x,
                         cdf_pareto = Y_pareto,
                         cdf_LN = Y_LN,
                         cdf_gamma = Y_gamma)
fitted_cdf_melted <- fitted_cdf %>%
  melt(id = c("loss"), variable.name = "Distribution")

names(fitted_cdf)[names(fitted_cdf)=="cdf_pareto"] <- "Pareto"
names(fitted_cdf)[names(fitted_cdf)=="cdf_LN"] <- "Log Normal"
names(fitted_cdf)[names(fitted_cdf)=="cdf_gamma"] <- "Gamma"

plot_CDF_comp <- ggplot() + stat_ecdf(data = claims, aes(loss)) +
  geom_line(data = fitted_cdf_melted, aes(x = loss, y = value,
                                          color=Distribution)) +
  labs(x = "Loss ($)", y = "Probability",
       title = "Fitted distributions against ECDF") +
  geom_vline(xintercept = claims_XS_mean$excess, alpha = 0.5,
             linetype = "dotted") + theme_bw()

if (output){
  ggsave(paste(report_path, "T1_PDF_comparison.pdf", sep = ""), plot_PDF_comp)
  ggsave(paste(report_path, "T1_CDF_comparison.pdf", sep = ""), plot_CDF_comp)
}

############# PP plot #################################
ordered_losses <- claims %>%
  select(loss) %>%
  arrange(loss) %>%
  mutate(theo_pareto = 0,
         theo_LN = 0,
         theo_gamma = 0)

ordered_losses$ecdf <- (1:n)/(n+1)
for (i in 1:n){
  ordered_losses$theo_pareto[i] <- F_pareto(ordered_losses$loss[i],
                                          coef(mle_pareto)[1],
                                          coef(mle_pareto)[2])
  ordered_losses$theo_LN[i] <- F_LN(ordered_losses$loss[i],
                                            coef(mle_LN)[1],
                                            coef(mle_LN)[2])
  ordered_losses$theo_gamma[i] <- F_gamma(ordered_losses$loss[i],
                                            coef(mle_gamma)[1],
                                            coef(mle_gamma)[2])
}

pp_pareto <- ggplot(data = ordered_losses, aes(x = theo_pareto, y = ecdf)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Theoretical CDF", y = "Empirical CDF", title = "P-P plot Pareto") +
  theme_bw()
pp_LN <- ggplot(data = ordered_losses, aes(x = theo_LN, y = ecdf)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Theoretical CDF", y = "Empirical CDF", title = "P-P plot Log Normal") +
  theme_bw()
pp_gamma <- ggplot(data = ordered_losses, aes(x = theo_gamma, y = ecdf)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Theoretical CDF", y = "Empirical CDF", title = "P-P plot Gamma") +
  theme_bw()
if (output){
  ggsave(paste(report_path, "PP_pareto.pdf", sep = ""), pp_pareto)
  ggsave(paste(report_path, "PP_LN.pdf", sep = ""), pp_LN)
  ggsave(paste(report_path, "PP_gamma.pdf", sep = ""), pp_gamma)
}

############# Q-Q plot ################################
merged_df <- ordered_losses %>% select(loss, ecdf) %>%
  merge(fitted_cdf %>% rename(quantile = loss), by = NULL)%>%
  mutate(diff_pareto = abs(ecdf - Pareto),
         diff_LN = abs(ecdf - `Log Normal`),
         diff_gamma = abs(ecdf - Gamma))

pareto_quantiles <- merged_df %>% select(loss, ecdf, diff_pareto) %>%
  group_by(loss, ecdf) %>%
  summarise(diff_pareto = min(diff_pareto))

pareto_quantiles_1 <- merge(pareto_quantiles,
                            merged_df %>% select(loss, quantile, diff_pareto),
                            by = c("diff_pareto", "loss"))

LN_quantiles <- merged_df %>% select(loss, ecdf, diff_LN) %>%
  group_by(loss, ecdf) %>%
  summarise(diff_LN = min(diff_LN))

LN_quantiles_1 <- merge(LN_quantiles,
                        merged_df %>% select(loss, quantile, diff_LN),
                        by = c("diff_LN", "loss"))

gamma_quantiles <- merged_df %>% select(loss, ecdf, diff_gamma) %>%
  group_by(loss, ecdf) %>%
  summarise(diff_gamma = min(diff_gamma))

gamma_quantiles_1 <- merge(gamma_quantiles,
                           merged_df %>% select(loss, quantile, diff_gamma),
                           by = c("diff_gamma", "loss"))


pareto_qq <- ggplot(data = pareto_quantiles_1, aes(x = quantile, y = loss)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, col = "red") +
  labs(title="Q-Q plot Pareto", x = "Theoretical quantile", y = "Observed quantile") +
  theme_bw()

LN_qq <- ggplot(data = LN_quantiles_1, aes(x = quantile, y = loss)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, col = "red") +
  labs(title="Q-Q plot Log Normal", x = "Theoretical quantile", y = "Observed quantile") +
  theme_bw()

gamma_qq <- ggplot(data = gamma_quantiles_1, aes(x = quantile, y = loss)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, col = "red") +
  labs(title="Q-Q plot Gamma", x = "Theoretical quantile", y = "Observed quantile") +
  theme_bw()

if (output){
  ggsave(paste(report_path, "QQ_pareto.pdf", sep = ""), pareto_qq)
  ggsave(paste(report_path, "QQ_LN.pdf", sep = ""), LN_qq)
  ggsave(paste(report_path, "QQ_gamma.pdf", sep = ""), gamma_qq)
}


########### Comparing statistics ######################
distribution_summary_stat <-
  data.frame(Measure = c("Parameter 1 estimate",
                        "Parameter 2 estimate",
                        "Negative log likelihood",
                        "AIC",
                        "BIC"),
            Pareto = c(coef(mle_pareto)[1],
                       coef(mle_pareto)[2],
                       f.nLL_pareto(coef(mle_pareto)[1],
                                    coef(mle_pareto)[2]),
                       AIC(mle_pareto),
                       AIC(mle_pareto, k = log(n))) %>%
              round(3),
            Log.Normal = c(coef(mle_LN)[1],
                             coef(mle_LN)[2],
                             f.nLL_LN(coef(mle_LN)[1],
                                          coef(mle_LN)[2]),
                             AIC(mle_LN),
                             AIC(mle_LN, k = log(n)))%>%
              round(3),
            Gamma = c(coef(mle_gamma)[1],
                      coef(mle_gamma)[2],
                      f.nLL_gamma(coef(mle_gamma)[1],
                                   coef(mle_gamma)[2]),
                      AIC(mle_gamma),
                      AIC(mle_gamma, k = log(n)))%>%
              round(3))

############################### Hypothesis tests ###########
##################### Chi-squared ##########################
# Choose bin points
bin_endpoints <- c(1000, seq(from = 2000, to = 10000, by = 500), c(20000))

# degrees of freedom = #bins - #parameters - 1
df <- (length(bin_endpoints) - 1) - 2 - 1
# Put losses into bins
grouped_losses <- claims %>%
  select(loss)
grouped_losses$binned <- cut(grouped_losses$loss, bin_endpoints)

# Count the number in each bin
grouped_losses_count <- grouped_losses %>%
  group_by(binned) %>%
  summarise(observations = n()) %>%
  mutate(pareto_theoretical = 0,
         ln_theoretical = 0,
         gamma_theoretical = 0)

# Compare with theoretical proportion in each bin
for (row in 1:nrow(grouped_losses_count)){
  # Pareto
  grouped_losses_count$pareto_theoretical[row] <- F_pareto(bin_endpoints[row + 1],
                                                           coef(mle_pareto)[1],
                                                           coef(mle_pareto)[2]) -
    F_pareto(bin_endpoints[row],
             coef(mle_pareto)[1],
             coef(mle_pareto)[2])

  # Log-Normal
  grouped_losses_count$ln_theoretical[row] <- F_LN(bin_endpoints[row + 1],
                                                           coef(mle_LN)[1],
                                                           coef(mle_LN)[2]) -
    F_LN(bin_endpoints[row],
         coef(mle_LN)[1],
         coef(mle_LN)[2])

  # Gamma
  grouped_losses_count$gamma_theoretical[row] <- F_gamma(bin_endpoints[row + 1],
                                                   coef(mle_gamma)[1],
                                                   coef(mle_gamma)[2]) -
    F_gamma(bin_endpoints[row],
            coef(mle_gamma)[1],
            coef(mle_gamma)[2])
}

# Find chi - squared test statistic for each distribution
X2_hypothesis_test <- grouped_losses_count %>%
  mutate(X2_pareto = (observations -
                        n * pareto_theoretical)^2/(n * pareto_theoretical),
         X2_ln = (observations - n * ln_theoretical)^2/(n * ln_theoretical),
         X2_gamma = (observations - n * gamma_theoretical)^2/
           (n * gamma_theoretical)) %>%
  select(X2_pareto, X2_ln, X2_gamma) %>%
  melt(variable.name = "Distribution") %>%
  group_by(Distribution) %>%
  summarise(Test_statistic = sum(value)) %>%
  mutate(p_value = 1 - pchisq(Test_statistic, df))


################## Kolmogorov Smirnoff #################
loss_ecdf <- ecdf(claims$loss)

### Pareto
KS_hypothesis_test <- data.frame(loss = claims$loss) %>%
  mutate(Pareto = 0,
         ln = 0,
         gamma = 0)

for (row in 1:nrow(KS_hypothesis_test)){
  KS_hypothesis_test$Pareto[row] <- abs(loss_ecdf(KS_hypothesis_test$loss[row]) -
                                          F_pareto(KS_hypothesis_test$loss[row],
                                                   coef(mle_pareto)[1],
                                                   coef(mle_pareto)[2]))

  KS_hypothesis_test$ln[row] <- abs(loss_ecdf(KS_hypothesis_test$loss[row]) -
                                          F_LN(KS_hypothesis_test$loss[row],
                                                   coef(mle_LN)[1],
                                                   coef(mle_LN)[2]))

  KS_hypothesis_test$gamma[row] <- abs(loss_ecdf(KS_hypothesis_test$loss[row]) -
                                          F_gamma(KS_hypothesis_test$loss[row],
                                                   coef(mle_gamma)[1],
                                                   coef(mle_gamma)[2]))
}

KS_statistics <- KS_hypothesis_test %>%
  select(-loss) %>%
  melt(variable.name = "Distribution") %>%
  group_by(Distribution) %>%
  summarise(KS_test_statistic = max(value))

# Critical values for alpha = 20%, 15%, 10%, 5%, 1%, 0.1%
KS_critical <- c(1.073, 1.138, 1.224, 1.358, 1.628, 1.94947)/sqrt(n)
KS_critical
})
