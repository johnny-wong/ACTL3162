# distribution function of the Burr XII(alpha,gam,theta) random variable
"pBurrXII" <- function(x,alpha,gam,theta)
{
	tmp <- 1/(1+(x/theta)^gam)
	result <- 1-tmp^alpha
	return(result)
}
