#aicc.R part of marginalmodelplots	
#by Andrew Redd
#Lisenced under GPL version 2 or newer.
aicc<-function (x, ..., pen=2) {
	m <- match.call()
	if (is.numeric(x)){ 
		m[[1]] <- as.name("locfit.raw")
	}else {
		m[[1]] <- as.name("locfit")
		names(m)[2] <- "formula"
	}
	m$pen <- NULL
	fit <- eval(m, sys.frame(sys.parent()))
	dp <- fit$dp
	z <- dp[c("lk", "df1", "df2")]
	n <- fit$mi["n"]
	z <- c(z, -2 * z[1] + pen*n*(z[2]+1)/(n-(z[2]+2)))
	names(z) <- c("lik", "infl", "vari", "aicc")
	z
}