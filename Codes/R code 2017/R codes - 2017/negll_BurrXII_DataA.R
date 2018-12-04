# negative log-likelihood of a Burr XII(alpha,gam,theta) distribution
# for fitting Valdez-DataSetA
"negll_BurrXII_DataA" <- function(parm, x) {
	tmp1 <- log(x)+(parm[1]+1)*log(1+(x/parm[3])^parm[2])
	tmp2 <- log(parm[1]*parm[2]) + parm[2]*log(x/parm[3])
	tmp3 <- tmp1-tmp2
	result <- sum(tmp3)
	return(result)
}
