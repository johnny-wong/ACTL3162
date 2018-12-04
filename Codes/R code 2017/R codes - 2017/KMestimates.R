# Kaplan-Meier estimates for truncated/censored data
"KMestimates" <- function(x,excess,censor)
{
	y0 <- min(excess)
	y.uncensored <- sort(c(x[which(censorj==0)],max(x)))
	r <- rep(0,length(y.uncensored))
	s <- rep(1,length(y.uncensored))
	for(i in 1:length(y.uncensored)) {
	  r[i] <- sum(excess < y.uncensored[i]) - sum(x < y.uncensored[i])
	}
	KM.F <- c(0,1-cumprod((r-s)/r))
	KM.F[length(KM.F)] <- ifelse(censor[which(x==max(x))]==1,NA,1)
	KM.y <- c(y0,y.uncensored)
	KM.r <- c(NA,r)
	KM.t <- c(y0,excess[which(censorj==0)],excess[which(x==max(x))])
	output <- cbind(KM.y,KM.t,KM.r,KM.F)
	col.titles <- c("KM.y","KM.t","KM.r","KM.F")
	colnames(output) <- col.titles
	return(data.frame(output))
}

		

	  
	