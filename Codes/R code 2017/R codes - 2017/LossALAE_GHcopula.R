# Data analysis of the Loss-ALAE data - fitting the Gumbel-Hougaard copula.
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

# to fit the Gumbel-Hougaard copula, need C.Gumbel, C1.Gumbel, C2.Gumbel and C12.Gumbel functions
"C.Gumbel" <- function(u1,u2,d.alpha)
{
	temp <- (((-log(u1))^d.alpha)+((-log(u2))^d.alpha))^(1/d.alpha)
	return(exp(-temp))
}
"C1.Gumbel" <- function(u1,u2,d.alpha)
{
	C <- C.Gumbel(u1,u2,d.alpha)
	temp <- (log(u1)/log(C))^(d.alpha-1)
	return(temp*(C/u1))
}
"C2.Gumbel" <- function(u1,u2,d.alpha)
{
	C <- C.Gumbel(u1,u2,d.alpha)
	temp <- (log(u2)/log(C))^(d.alpha-1)
	return(temp*(C/u2))
}
"C12.Gumbel" <- function(u1,u2,d.alpha)
{
	C <- C.Gumbel(u1,u2,d.alpha)
	C1 <- C1.Gumbel(u1,u2,d.alpha)
	C2 <- C2.Gumbel(u1,u2,d.alpha)
	temp <- 1+(d.alpha-1)/-log(C)
	return(temp*C1*C2/C)
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
	f12 <- C12.Gumbel(u1,u2,d.alpha)
	n.censor <- f1*f2*f12
	c2 <- C2.Gumbel(u1,u2,d.alpha)
	y.censor <- f2*(1-c2)
	log.ll <- (1-censor)*log(n.censor)+censor*log(y.censor)
	return(sum(-log.ll))
}
init.est <- c(14000,1.1,14000,2.1,1.4)

fit.Gumbel <- optim(init.est, neg.loglik,NULL,x1=LOSS,x2=ALAE,censor=CENSOR)
parm.hat <- fit.Gumbel$par
loglik.Gumbel <- -fit.Gumbel$value
AIC.Gumbel <- (-2*loglik.Gumbel+2*length(parm.hat))/length(LOSS)
print(loglik.Gumbel)
print(fit.Gumbel)
print(AIC.Gumbel)

# next estimate the standard errors.
library(nlme)
negll.Gumbel.Hess <- fdHess(parm.hat, neg.loglik, x1=LOSS,x2=ALAE,censor=CENSOR)
inv.Gumbel.Hess <- solve(negll.Gumbel.Hess$Hessian)
parm.se <- sqrt(diag(inv.Gumbel.Hess))
output <- cbind(parm.hat,parm.se)
output <- round(output,digits=3)
rownames(output)<- c("lambda1","theta1","lambda2","theta2","d.alpha")
colnames(output)<- c("estimate","std error")
print(output)
