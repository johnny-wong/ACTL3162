# quantile function of the Burr XII(alpha,gam,theta) random variable
"qBurrXII" <- function(q,alpha,gam,theta)
{
	tmp <- ((1-q)^(-1/alpha))-1
	result <- theta*tmp^(1/gam)
	return(result)
}
