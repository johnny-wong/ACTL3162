# negative log-likelihood of a Gamma distribution
# for fitting Valdez-DataSetA
"negll_gamma_DataA" <- function(parm, x) {
	tmp1 <- parm[2]*x + log(gamma(parm[1]))
	tmp2 <- parm[1]*log(parm[2]) + (parm[1]-1)*log(x)
	result <- sum(tmp1-tmp2)
	return(result)
}
