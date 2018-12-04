# Data analysis of the Loss-ALAE data - fitting the Frank copula.
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

# to fit the Frank copula, need C.Frank, C1.Frank, C2.Frank and C12.Frank functions
"C12.Frank" <- function(u1,u2,d.alpha)
{
	tmp1 <- exp(d.alpha)-1
	tmp2 <- exp(d.alpha*(u1+u2))
	tmp3 <- (exp(d.alpha*u1)-1)*(exp(d.alpha*u2)-1)
	num <- d.alpha*tmp1*tmp2
	den <- (tmp1+tmp3)^2
	return(num/den)
}
"C2.Frank" <- function(u1,u2,d.alpha)
{
	tmp1 <- exp(d.alpha)-1
	tmp2 <- (exp(d.alpha*u1)-1)*(exp(d.alpha*u2)-1)
	num <- exp(d.alpha*u2)*(exp(d.alpha*u1)-1)
	den <- tmp1+tmp2
	return(num/den)
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
	f12 <- C12.Frank(u1,u2,d.alpha)
	n.censor <- f1*f2*f12
	c2 <- C2.Frank(u1,u2,d.alpha)
	y.censor <- f2*(1-c2)
	log.ll <- (1-censor)*log(n.censor)+censor*log(y.censor)
	return(sum(-log.ll))
}
init.est <- c(14000,1.1,16000,2.3,-3.1)

fit.Frank <- optim(init.est, neg.loglik,NULL,x1=LOSS,x2=ALAE,censor=CENSOR)
parm.hat <- fit.Frank$par
loglik.Frank <- -fit.Frank$value
AIC.Frank <- (-2*loglik.Frank+2*length(parm.hat))/length(LOSS)
print(loglik.Frank)
print(fit.Frank)
print(AIC.Frank)

# next estimate the standard errors.
library(nlme)
negll.Frank.Hess <- fdHess(parm.hat, neg.loglik, x1=LOSS,x2=ALAE,censor=CENSOR)
inv.Frank.Hess <- solve(negll.Frank.Hess$Hessian)
parm.se <- sqrt(diag(inv.Frank.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=3)
rownames(output)<- c("lambda1","theta1","lambda2","theta2","d.alpha")
colnames(output)<- c("estimate","std error")
print(output)
