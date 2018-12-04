# negative log-likelihood of a log-normal distribution for the Dataset B
"negll_lognorm_DataB" <- function(parm, x, excess, censor) {
	numcens0 <- log(dlnorm(x,meanlog=parm[1],sdlog=parm[2]))
	numcens1 <- log(1-plnorm(x,meanlog=parm[1],sdlog=parm[2]))
	den <- log(1-plnorm(excess,meanlog=parm[1],sdlog=parm[2]))
	tmp1 <- (1-censor)*(numcens0 - den)
	tmp2 <- censor*(numcens1 - den)
	result <- -sum(tmp1+tmp2)
	return(result)
}
