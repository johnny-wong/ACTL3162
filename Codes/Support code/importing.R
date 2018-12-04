# Imports data and saves as an R image
setwd("C:/Users/johnn/Dropbox/Uni/2017 S2/ACTL3162/Assignment/Codes")

claims <- read.csv(file = "data.csv.csv")

# General data manipulation
default_excess <- 700
x <- claims$paid
claims$excess <- default_excess + claims$elect_XS
claims$loss <- claims$paid + claims$excess
n <- nrow(claims)

save.image(file = "claims.RData")
