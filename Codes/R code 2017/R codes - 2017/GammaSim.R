# Copyright E.A. Valdez, Aug 2006
# Simulating from Cook-Johnson copula and Normal copula with Gamma(5,1) marginals.
# Same correlation assumption, but two different dependence structure

# set the number of pairs to generate
n.sim <- 5000
tau <- 0.5398931
x.cook <- rep(0,n.sim)
y.cook <- rep(0,n.sim)
x.norm <- rep(0,n.sim)
y.norm <- rep(0,n.sim)
# first, we simulate from a Cook-Johnson
theta <- 2*tau/(1-tau)
for (i in 1:n.sim){
# generate independent Exponential(1)
y.exp <- rexp(2)
# generate a Gamma(1/theta,1)
z.gam <- rgamma(1,shape=1/theta)
ucook.vector <- (1+y.exp/z.gam)^(-1/theta)
xcook.vector <- qgamma(ucook.vector,shape=5)
x.cook[i] <- xcook.vector[1]
y.cook[i] <- xcook.vector[2]
}
x95.cook <- quantile(x.cook,0.95)
y95.cook <- quantile(y.cook,0.95)

# next, we simulate from a Normal copula
rho <- sin(pi*tau/2)
print(rho)
corr.mat <- matrix(c(1,rho,rho,1),nrow=2,ncol=2)
print(corr.mat)
b.mat <- chol(corr.mat)
print(b.mat)
for (i in 1:n.sim){
# generate independent Normal
z.norm <- rnorm(2)
w.norm <- t(b.mat) %*% z.norm
unorm.vector <- pnorm(w.norm)
xnorm.vector <- qgamma(unorm.vector,shape=5)
x.norm[i] <- xnorm.vector[1]
y.norm[i] <- xnorm.vector[2]
}
x95.norm <- quantile(x.norm,0.95)
y95.norm <- quantile(y.norm,0.95)
 
par(mfrow=c(1,2))
plot(x.cook,y.cook,xlab="x~Gamma(5,1)",ylab="y~Gamma(5,1)",xlim=c(0,16),ylim=c(0,16),main="Cook-Johnson copula",pch=20)
abline(h=y95.cook,col="blue")
abline(v=x95.cook,col="blue")
plot(x.norm,y.norm,xlab="x~Gamma(5,1)",ylab="y~Gamma(5,1)",xlim=c(0,16),ylim=c(0,16),main="Normal copula",pch=20)
abline(h=y95.norm,col="blue")
abline(v=x95.norm,col="blue")
