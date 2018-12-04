# density of the Burr XII(alpha,gam,theta) random variable
"dBurrXII" <- function(x,alpha,gam,theta)
{
	num <- alpha*gam*((x/theta)^gam)
	den <- x*((1+(x/theta)^gam)^(alpha+1))
	result <- num/den
	return(result)
}
