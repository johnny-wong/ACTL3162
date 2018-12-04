# Fitting the log-normal distribution to insurance claims data (using Dataset B).
#
# Parameters are called alpha (parm[1]) and beta (parm[2]).

# read the data file
InsClaim <- read.csv("Klugman-DataSetB-modified.csv")
attach(InsClaim)

source("negll_gamma_DataB.R")
source("KMestimates.R")
source("interpol_fn.R")


# now find the parameter estimates using constrained optimization
# first set the initial parameter estimates

samp.mu <- mean(claimj)
samp.var <- sd(claimj)^2
beta.tilde <- samp.mu/samp.var
alpha.tilde <- samp.mu*beta.tilde
init.est <- c(alpha.tilde,beta.tilde)

fit.gamma <- constrOptim(init.est, negll_gamma_DataB, NULL, ui=c(0,1), ci=c(0,0), x=claimj, excess=excessj, censor=censorj)
parm.hat <- fit.gamma$par
loglik.gamma <- -fit.gamma$value
alpha.hat <- parm.hat[1]
beta.hat <- parm.hat[2]

# next estimate the standard errors.
library(nlme)
negll.gamma.Hess <- fdHess(parm.hat, negll_gamma_DataB, x=claimj, excess=excessj, censor=censorj)
inv.gamma.Hess <- solve(negll.gamma.Hess$Hessian)
parm.se <- sqrt(diag(inv.gamma.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=4)
rownames(output)<- c("alpha","beta")
colnames(output)<- c("estimate","std error")
print(output)

# do some graphical tests
par(mfrow=c(2,2))

# do the Kaplan-Meier estimates
KMestimates <- KMestimates(x=claimj, excess=excessj, censor=censorj)
attach(KMestimates)

# the density plot
hist(claimj,breaks=30,prob=T)
lines(sort(claimj),dgamma(sort(claimj),shape=alpha.hat,rate=beta.hat))

# the empirical distribution function (uses Kaplan-Meier estimates)
plot(KM.y,KM.F,type="s",xlab="claimj",ylab="cdf",main="Kaplan-Meier CDF")
abline(1,0,lty=3)
abline(0,0,lty=3)
lines(KM.y,pgamma(KM.y,shape=alpha.hat,rate=beta.hat))

# the quantile-quantile plot
plot(qgamma(KM.F,shape=alpha.hat,rate=beta.hat),sort(KM.y),xlim=c(0,6500),ylim=c(0,6500),
   xlab="theoretical quantiles",ylab="sample quantiles",main="Q-Q plot",cex=0.55)
abline(0,1)

# the probability-probability plot
plot(pgamma(sort(KM.y),shape=alpha.hat,rate=beta.hat),KM.F,
   xlab="theoretical probability",ylab="sample probability",main="P-P plot",cex=0.55)
abline(0,1)

# now compute test statistics

# the Kolmogorov-Smirnoff statistic
# there are 3 truncation points
th <- c(100,250,500)
count.th <- c(30,40,30)
KM.th <- interpol_fn(th,KM.y[which(KM.F>=0)],KM.F[which(KM.F>=0)])
nh <- cumsum(c(30*1,40*0.93333,30*0.85098))
# nh <- cumsum((1-KM.th)*count.th)
n <- max(nh)
nj <- rep(0,length(KM.t[-1]))
for(i in 1:length(KM.t[-1])) {
  nj[i] <- nh[1]*(KM.t[i+1]==th[1])+nh[2]*(KM.t[i+1]==th[2])+nh[3]*(KM.t[i+1]==th[3])
}
samp.pct <- KM.F[which(KM.F>=0)]
# theor.pct <- pgamma(KM.y[which(KM.F>=0)],shape=alpha.hat,rate=beta.hat)/(1-pgamma(KM.t[which(KM.F>=0)],shape=alpha.hat,rate=beta.hat))
theor.pct <- pgamma(KM.y[which(KM.F>=0)],shape=alpha.hat,rate=beta.hat)
D.stats <- sqrt(n)*abs(samp.pct-theor.pct) + 0.19/sqrt(n)
KStest.stat <- max(D.stats)
names(KStest.stat) <- c("Kolmogorov Smirnoff test statistic")
print(KStest.stat)

# the Anderson-Darling test statistic
theor.pct <- c(theor.pct,max(theor.pct))
tmp1 <- nj*((1-samp.pct)^2)*(log(1-theor.pct[-length(theor.pct)])-log(1-theor.pct[-1]))
tmp2 <- nj*(samp.pct^2)*log(theor.pct[-1]/theor.pct[-length(theor.pct)])
tmp3 <- nj*(theor.pct[-1]-theor.pct[-length(theor.pct)])
ADtest.stat <- sum(tmp1) + sum(tmp2) - sum(tmp3)
names(ADtest.stat) <- c("Anderson-Darling test statistic")
print(ADtest.stat)

# the chi-square test statistic
# subjective choice of intervals
c.interval <- c(0,100,184,424,561,707,904,1028,1138,1215,1379,1557,1742,1899,2675,3634,4409)
nj.interval <- c(nh[1],nh[2],rep(nh[3],14))
numb.par <- 2
deg.f <- length(c.interval)-1-numb.par-1
claimj.freq <- table(cut(claimj,breaks=c.interval,dig.lab=5))
probj <- pgamma(c.interval,shape=alpha.hat,rate=beta.hat)
pj <- diff(probj)
pj <- c(pj[-1],1-sum(pj[-1]))
pj.hat <- rep(1/length(c.interval[-1]),length(c.interval[-1]))
chisq.claimj <- nj.interval*(pj-pj.hat)^2/pj
print(cbind(claimj.freq,nj.interval,pj,pj.hat,chisq.claimj))
chisqtest.stat <- sum(chisq.claimj)
chi.pvalue <- pchisq(chisqtest.stat, df=deg.f, lower.tail=FALSE)
names(chisqtest.stat) <- c("Chi-square test statistic")
print(chisqtest.stat)
names(chi.pvalue) <- c("Chi-square test p-value")
print(chi.pvalue)

# the SBC (Schwarz Bayesian Criterion) test
numb.par <- 2
SBC.stat <- loglik.gamma - (numb.par/2)*log(length(claimj))
out <- rbind(loglik.gamma,SBC.stat)
colnames(out) <- c("value")
rownames(out) <- c("Neg Log Likelihood","SBC criterion")
print(out)

detach(InsClaim)