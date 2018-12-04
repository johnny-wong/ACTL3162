# function for interpolating from a non-smooth function y of x.
"interpol_fn" <- function(c,x,y)
{
	sort.x <- sort(x)
	sort.y <- sort(y)
	c.length <- length(c)
	out <- rep(0,c.length)
	for(i in 1:c.length){
	if (c[i]==0) out[i]=0 else
	if (c[i]<min(x)) out[i]=min(y) else
	if (c[i]>max(x)) out[i]=max(y) else
	{	loc1 <- length(sort.x[which(sort.x<=c[i])])
	out[i] <- sort.y[loc1] + (c[i] - sort.x[loc1])*((sort.y[loc1+1]-sort.y[loc1])/(sort.x[loc1+1]-sort.x[loc1]))
	}
	}
	return(out)
}
