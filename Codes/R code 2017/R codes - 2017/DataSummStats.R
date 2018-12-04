"DataSummStats" <- function(data)
{
	# gives the following summary statistics
	# count, mean, median,Q1/Q3, variance, sd, min, max, skewness, kurtosis
	row.titles <- c("Number", "Mean", "5th Q", "25th Q", "Median", "75th Q" , "95th Q",
			"Variance", "StdDev", "Minimum", "Maximum", "Skewness", "Kurtosis")
	output <- c(1:13)
	output[1] <- length(data)
	output[2] <- mean(data)
	output[3] <- quantile(data,0.05)
	output[4] <- quantile(data,0.25)
	output[5] <- quantile(data,0.5)
	output[6] <- quantile(data,0.75)
	output[7] <- quantile(data,0.95)
	output[8] <- var(data)
	output[9] <- sd(data)
	output[10] <- min(data)
	output[11] <- max(data)
	output[12] <- (sum((data-mean(data))^3)/var(data)^1.5) /output[1]
	output[13] <- (sum((data-mean(data))^4)/var(data)^2) /output[1] - 3
	output <- round(output,2)
	output<-cbind(output)
	rownames(output)<-c(row.titles)
	colnames(output)<-c("Value")
	print(output,print.gap="3")
}
