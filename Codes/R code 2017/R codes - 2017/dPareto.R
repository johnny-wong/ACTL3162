# density of the Pareto(alpha,x0) random variable
"dPareto" <- function(x,alpha,x0)
{
	num <- alpha*(x0^alpha)
	den <- (x+x0)^(alpha+1)
	result <- num/den
	return(result)
}
