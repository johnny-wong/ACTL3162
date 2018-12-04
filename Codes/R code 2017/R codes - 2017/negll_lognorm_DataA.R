# negative log-likelihood of a log-normal distribution
# for fitting Valdez-DataSetA
"negll_lognorm_DataA" <- function(parm, x) {
	tmp1 <- ((log(x)-parm[1])^2)/(2*parm[2]^2)
	tmp2 <- log(x*sqrt(2*pi)*parm[2])
	result <- sum(tmp1+tmp2)
	return(result)
}
