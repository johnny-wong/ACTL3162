# negative log-likelihood of a gamma distribution for the Dataset B
"negll_gamma_DataB" <- function(parm, x, excess, censor) {
	numcens0 <- log(dgamma(x,shape=parm[1],rate=parm[2]))
	numcens1 <- log(1-pgamma(x,shape=parm[1],rate=parm[2]))
	den <- log(1-pgamma(excess,shape=parm[1],rate=parm[2]))
	tmp1 <- (1-censor)*(numcens0 - den)
	tmp2 <- censor*(numcens1 - den)
	result <- -sum(tmp1+tmp2)
	return(result)
}
