# Fitting the Gamma distribution to insurance claims data.
#
# Parameters are called alpha (parm[1]) and beta (parm[2]).

# read the data file

InsClaim <- read.csv("Valdez-DataSetA.csv")
attach(InsClaim)

source("negll_gamma_DataA.R")


samp.mu <- mean(claimj)
samp.var <- sd(claimj)^2
beta.tilde <- samp.mu/samp.var
alpha.tilde <- samp.mu*beta.tilde
init.est <- c(alpha.tilde,beta.tilde)

fit.gamma <- constrOptim(init.est, negll_gamma_DataA, NULL, ui=rbind(c(1,0),c(0,1)), ci=c(0,0), x=claimj)
parm.hat <- fit.gamma$par
loglik.gamma <- -fit.gamma$value
alpha.hat <- parm.hat[1]
beta.hat <- parm.hat[2]

# next estimate the standard errors.
library(nlme)
negll.gamma.Hess <- fdHess(parm.hat, negll_gamma_DataA, x=claimj)
inv.gamma.Hess <- solve(negll.gamma.Hess$Hessian)
parm.se <- sqrt(diag(inv.gamma.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=4)
rownames(output)<- c("alpha","beta")
colnames(output)<- c("estimate","std error")
print(output)

# do some graphical tests
par(mfrow=c(2,2))

# the density plot
hist(claimj,breaks=45,prob=T)
lines(sort(claimj),dgamma(sort(claimj),shape=alpha.hat,rate=beta.hat))

# the empirical distribution function
samp.pct <- (1:length(claimj)-0.5)/length(claimj)
plot(sort(claimj),samp.pct,type="s",xlab="claimj",ylab="cdf",main="Empirical CDF")
abline(1,0,lty=3)
abline(0,0,lty=3)
lines(sort(claimj),pgamma(sort(claimj),shape=alpha.hat,rate=beta.hat))

# the quantile-quantile plot
plot(qgamma(samp.pct,shape=alpha.hat,rate=beta.hat),sort(claimj),xlim=c(0,6500),ylim=c(0,6500),
   xlab="theoretical quantiles",ylab="sample quantiles",main="Q-Q plot",cex=0.55)
abline(0,1)

# the probability-probability plot
plot(pgamma(sort(claimj),shape=alpha.hat,rate=beta.hat),samp.pct,
   xlab="theoretical probability",ylab="sample probability",main="P-P plot",cex=0.55)
abline(0,1)

# now compute test statistics

# the Kolmogorov-Smirnoff statistic
samp.pct <- (1:length(claimj))/length(claimj)
theor.pct <- pgamma(sort(claimj),shape=alpha.hat,rate=beta.hat)
D.stats <- abs(samp.pct-theor.pct)
KStest.stat <- max(D.stats)
names(KStest.stat) <- c("Kolmogorov Smirnoff test statistic")
print(KStest.stat)

# the Anderson-Darling test statistic
tmp1 <- (1-samp.pct[-length(samp.pct)])^2
tmp2 <- log(1-theor.pct[-length(theor.pct)])-log(1-theor.pct[-1])
tmp3 <- (samp.pct[-length(samp.pct)]^2)*log(theor.pct[-1]/theor.pct[-length(theor.pct)])
tmp4 <- theor.pct[-length(-1)]-theor.pct[-1]
ADtest.stat <- sum(tmp1*tmp2+tmp3-tmp4)
names(ADtest.stat) <- c("Anderson-Darling test statistic")
print(ADtest.stat)

# the chi-square test statistic
# subjective choice of intervals
c.interval <- seq(0,6000,by=300)
numb.par <- 2
deg.f <- length(c.interval)-1-numb.par-1
claimj.freq <- table(cut(sort(claimj),breaks=c.interval,dig.lab=5))
prob.claimj <- pgamma(c.interval,shape=alpha.hat,rate=beta.hat)
exp.claimj <- round(diff(prob.claimj)*length(claimj),4)
chisq.claimj <- (claimj.freq-exp.claimj)^2/exp.claimj
print(cbind(claimj.freq,exp.claimj,chisq.claimj))
chisqtest.stat <- sum(chisq.claimj)
chi.pvalue <- pchisq(chisqtest.stat, df=deg.f, lower.tail=FALSE)
names(chisqtest.stat) <- c("Chi-square test statistic")
print(chisqtest.stat)
names(chi.pvalue) <- c("Chi-square test p-value")
print(chi.pvalue)

# the SBC (Schwarz Bayesian Criterion) test
SBC.stat <- loglik.gamma - (numb.par/2)*log(length(claimj))
out <- rbind(loglik.gamma,SBC.stat)
colnames(out) <- c("value")
rownames(out) <- c("Neg Log Likelihood","SBC criterion")
print(out)

detach(InsClaim)
