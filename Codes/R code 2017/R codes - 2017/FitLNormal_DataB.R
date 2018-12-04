# Fitting the log-normal distribution to insurance claims data (using Dataset B).
#
# Parameters are called mulog (parm[1]) and sigmalog (parm[2]).

# read the data file
setwd("C:/Users/lenovo/Dropbox/ACTL3003 stuff/Lecture slides/Module 3 (with 2010)/11_IRM_3_Rcode")

InsClaim <- read.csv("Klugman-DataSetB-modified.csv")
attach(InsClaim)

source("negll_lognorm_DataB.R")
source("KMestimates.R")
source("interpol_fn.R")

# now find the parameter estimates using constrained optimization
# first set the initial parameter estimates

samp.mulog <- mean(log(claimj))
samp.sdlog <- sd(log(claimj))
init.est <- c(samp.mulog,samp.sdlog)

fit.lognorm <- constrOptim(init.est, negll_lognorm_DataB, NULL, ui=c(0,1), ci=c(0,0), x=claimj, excess=excessj, censor=censorj)
parm.hat <- fit.lognorm$par
loglik.lognorm <- -fit.lognorm$value
mulog.hat <- parm.hat[1]
sdlog.hat <- parm.hat[2]

# next estimate the standard errors.
library(nlme)
negll.lognorm.Hess <- fdHess(parm.hat, negll_lognorm_DataB, x=claimj, excess=excessj, censor=censorj)
inv.lognorm.Hess <- solve(negll.lognorm.Hess$Hessian)
parm.se <- sqrt(diag(inv.lognorm.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=4)
rownames(output)<- c("mulog","sdlog")
colnames(output)<- c("estimate","std error")
print(output)

# do some graphical tests
par(mfrow=c(2,2))

# do the Kaplan-Meier estimates
KMestimates <- KMestimates(x=claimj, excess=excessj, censor=censorj)
attach(KMestimates)

# the density plot
hist(claimj,breaks=30,prob=T)
lines(sort(claimj),dlnorm(sort(claimj),meanlog=mulog.hat,sdlog=sdlog.hat))

# the empirical distribution function (uses Kaplan-Meier estimates)
plot(KM.y,KM.F,type="s",xlab="claimj",ylab="cdf",main="Kaplan-Meier CDF")
abline(1,0,lty=3)
abline(0,0,lty=3)
lines(KM.y,plnorm(KM.y,meanlog=mulog.hat,sdlog=sdlog.hat))

# the quantile-quantile plot
plot(qlnorm(KM.F,meanlog=mulog.hat,sdlog=sdlog.hat),sort(KM.y),xlim=c(0,6500),ylim=c(0,6500),
   xlab="theoretical quantiles",ylab="sample quantiles",main="Q-Q plot",cex=0.55)
abline(0,1)

# the probability-probability plot
plot(plnorm(sort(KM.y),meanlog=mulog.hat,sdlog=sdlog.hat),KM.F,
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
# theor.pct <- plnorm(KM.y[which(KM.F>=0)],meanlog=mulog.hat,sdlog=sdlog.hat)/(1-plnorm(KM.t[which(KM.F>=0)],meanlog=mulog.hat,sdlog=sdlog.hat))
theor.pct <- plnorm(KM.y[which(KM.F>=0)],meanlog=mulog.hat,sdlog=sdlog.hat)
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
probj <- plnorm(c.interval,meanlog=mulog.hat,sdlog=sdlog.hat)
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
SBC.stat <- loglik.lognorm - (numb.par/2)*log(length(claimj))
out <- rbind(loglik.lognorm,SBC.stat)
colnames(out) <- c("value")
rownames(out) <- c("Neg Log Likelihood","SBC criterion")
print(out)

detach(InsClaim)
detach(KMestimates)
