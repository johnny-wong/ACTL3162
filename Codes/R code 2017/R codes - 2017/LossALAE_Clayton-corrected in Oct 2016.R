# Data analysis of the Loss-ALAE data - fitting the Clayton copula.
#
# read the data file
Loss.ALAE <- read.csv("LossData-FV.csv")
attach(Loss.ALAE)

source("DataSummStats.R")
source("dPareto.R")
source("pPareto.R")


# understand and summarize the data.
DataSummStats(LOSS)
DataSummStats(ALAE)
DataSummStats(LIMIT)
DataSummStats(LOSS[which(CENSOR==0)])
DataSummStats(LOSS[which(CENSOR==1)])

# plot LOSS vs ALAE
plot(log(LOSS),log(ALAE),main="LOSS vs ALAE on a log scale",pch=19)

# to fit the Clayton copula, need C.Clayton, C1.Clayton, C2.Clayton and C12.Clayton functions
"C.Clayton" <- function(u1,u2,d.alpha)
{
	tmp1 <- u1^(-d.alpha)
	tmp2 <- u2^(-d.alpha)
	return((tmp1+tmp2-1)^(-1/d.alpha))
}
"C2.Clayton" <- function(u1,u2,d.alpha)
{
	C <- C.Clayton(u1,u2,d.alpha)
	return((C/u2)^(d.alpha+1))
}
"C12.Clayton" <- function(u1,u2,d.alpha)
{
	C <- C.Clayton(u1,u2,d.alpha)
	tmp1 <- (d.alpha+1)*(C^d.alpha)
	tmp2 <- (C/(u1*u2))^(d.alpha+1)
	return(tmp1*tmp2)
}

"neg.loglik" <- function(parm,x1,x2,censor)
{
	lambda1 <- parm[1]
	theta1 <- parm[2]
	lambda2 <- parm[3]
	theta2 <- parm[4]
	d.alpha <- parm[5]
	f1 <- dPareto(x1,alpha=theta1,x0=lambda1)
	f2 <- dPareto(x2,alpha=theta2,x0=lambda2)
	u1 <- pPareto(x1,alpha=theta1,x0=lambda1)
	u2 <- pPareto(x2,alpha=theta2,x0=lambda2)
	f12 <- C12.Clayton(u1,u2,d.alpha)
	n.censor <- f1*f2*f12
	c2 <- C2.Clayton(u1,u2,d.alpha)
	y.censor <- f2*(1-c2)
	log.ll <- (1-censor)*log(n.censor)+censor*log(y.censor)
	return(sum(-log.ll))
}
init.est <- c(14000,1.1,16000,2.3,1.5)

fit.Clayton <- optim(init.est, neg.loglik,NULL,x1=LOSS,x2=ALAE,censor=CENSOR)
parm.hat <- fit.Clayton$par
loglik.Clayton <- -fit.Clayton$value
AIC.Clayton <- (-2*loglik.Clayton+2*length(parm.hat))/length(LOSS)
print(loglik.Clayton)
print(fit.Clayton)
print(AIC.Clayton)

# next estimate the standard errors.
library(nlme)
negll.Clayton.Hess <- fdHess(parm.hat, neg.loglik, x1=LOSS,x2=ALAE,censor=CENSOR)
inv.Clayton.Hess <- solve(negll.Clayton.Hess$Hessian)
parm.se <- sqrt(diag(inv.Clayton.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=3)
rownames(output)<- c("lambda1","theta1","lambda2","theta2","d.alpha")
colnames(output)<- c("estimate","std error")
print(output)
