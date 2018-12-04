# Data analysis of the Loss-ALAE data - fitting the Independence copula.
#
# read the data file

Loss.ALAE <- read.csv("LossData-FV.csv")
attach(Loss.ALAE)

source("DataSummStats")
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

# to fit the Independent copula
"neg.loglik" <- function(parm,x1,x2,censor)
{
	lambda1 <- parm[1]
	theta1 <- parm[2]
	lambda2 <- parm[3]
	theta2 <- parm[4]
	f1 <- dPareto(x1,alpha=theta1,x0=lambda1)
	f2 <- dPareto(x2,alpha=theta2,x0=lambda2)
	u1 <- pPareto(x1,alpha=theta1,x0=lambda1)
	u2 <- pPareto(x2,alpha=theta2,x0=lambda2)
	n.censor <- f1*f2
	c2 <- u1
	y.censor <- f2*(1-c2)
	log.ll <- (1-censor)*log(n.censor)+censor*log(y.censor)
	return(sum(-log.ll))
}
init.est <- c(14000,1.1,16000,2.3)

fit.Indep <- optim(init.est, neg.loglik,NULL,x1=LOSS,x2=ALAE,censor=CENSOR)
parm.hat <- fit.Indep$par
loglik.Indep <- -fit.Indep$value
AIC.Indep <- (-2*loglik.Indep+2*length(parm.hat))/length(LOSS)
print(loglik.Indep)
print(fit.Indep)
print(AIC.Indep)

# next estimate the standard errors.
library(nlme)
negll.Indep.Hess <- fdHess(parm.hat, neg.loglik, x1=LOSS,x2=ALAE,censor=CENSOR)
inv.Indep.Hess <- solve(negll.Indep.Hess$Hessian)
parm.se <- sqrt(diag(inv.Indep.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=3)
rownames(output)<- c("lambda1","theta1","lambda2","theta2")
colnames(output)<- c("estimate","std error")
print(output)
