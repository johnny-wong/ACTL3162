# distribution function of the Pareto(alpha,x0) random variable
"pPareto" <- function(x,alpha,x0)
{
	result <- 1-(x0/(x+x0))^alpha
	return(result)
}
