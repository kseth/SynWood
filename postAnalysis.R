library(MASS)
source("extrapol_field.R")

nameSimul <- "name_LLJ100"
daisyChainSeeds <- 1000:1048

outfiles <- paste0("completethetasamples_all", daisyChainSeeds, ".txt")
allRuns <- read.table(file = outfiles[1], header = TRUE)
allLengths <- dim(allRuns)[1]

for(nums in 2:length(outfiles)){
	allRuns <- rbind(allRuns, read.table(file = outfiles[nums], header = TRUE))
	allLengths <- c(allLengths, dim(allRuns)[1])
}

names(allRuns) <- { if(dim(allRuns)[2] == 4) c("LL", "LP", "rateMove", "weightJumpInMove") else c("LL", "LP", "rateMove", "weightJumpInMove", "detectRate") }
traces(allRuns)

realMeans <- c(0.04, 0.10, 1)
names(realMeans) <- c("rateMove", "weightJumpInMove", "detectRate")
lims <- data.frame("rateMove" = c(0, .08), "weightJumpInMove" = c(0, 1), "detectRate" = c(0, 1))

dev.new()
par(mfrow=c(dim(allRuns)[2],1))

for(param in names(realMeans)){
	if(param %in% names(allRuns))
		get.estimate(allRuns[[param]],true.val=realMeans[param], name=param, xlim = lims[, param])
}

totalWeight <- realMeans["weightJumpInMove"] + 1
true.val2 <- c(1/totalWeight, realMeans["weightJumpInMove"]/totalWeight)
names(true.val2) <- c("rateHop", "rateJump")
xlim <- c(0, 1)
totalWeight <- allRuns$weightJumpInMove + 1
rateHop <- 1/totalWeight
rateJump <- allRuns$weightJumpInMove/totalWeight

get.estimate(rateHop, true.val = true.val2["rateHop"], name = "rateHop", xlim = xlim)
get.estimate(rateJump, true.val = true.val2["rateJump"], name = "rateJump", xlim = xlim)

dev.copy2pdf(file = paste0(nameSimul, "_comb_post.pdf"))
graphics.off()


## determine mean(sd(eachRun))/sd(allRuns)
## determine mean(zscore(eachRun))/zscore(allRuns)
sd_all <- apply(allRuns[, -(1:2)], 2, sd)
mean_all <- apply(allRuns[, -(1:2)], 2, mean)
zscore_all <- (realMeans[names(realMeans) %in% names(allRuns)] - mean_all)/sd_all

#determine the standard deviation, mean, median, zscores and quantiles (0.025, 0.975) for each of the runs

sd_each <- apply(allRuns[1:allLengths[1], -(1:2)], 2, sd)
sd_each <- data.frame(t(sd_each))

mean_each <- apply(allRuns[1:allLengths[1], -(1:2)], 2, mean)
mean_each <- data.frame(t(mean_each))

median_each <- apply(allRuns[1:allLengths[1], -(1:2)], 2, median)
median_each <- data.frame(t(median_each))

zscore_each <- (realMeans[names(realMeans) %in% names(allRuns)] - mean_each)/sd_each
zscore_each <- data.frame((zscore_each))

quantile_temp <- apply(allRuns[1:allLengths[1], -(1:2)], 2, quantile, probs = c(0.025, 0.975))

quantile_each <- t(as.vector(quantile_temp))
quantile_each <- data.frame(quantile_each)

#using the counts, determine how many times we fall outside of the confidence interval
counts <- rep(0, length(realMeans))
names(counts) <- names(realMeans)

for(param in names(realMeans)[names(realMeans)%in%names(allRuns)]){
	if(realMeans[param] >= quantile_temp[1, param] && realMeans[param] <= quantile_temp[2, param])
		counts[param] <- counts[param]+1
}


#do Cook's p
prior_q_distribution <- list()
length(prior_q_distribution) <- length(realMeans)
names(prior_q_distribution) <- names(realMeans)  

for(param in names(realMeans)[names(realMeans)%in%names(allRuns)]){
	qdist <- ecdf(allRuns[1:allLengths[1], param]) # find the emperical cdf
	quant <- qdist(realMeans[param]) # find the quantile of the real value
	if(quant == 0)
		quant <- 1/(allLengths[1]+1)
	if(quant == 1)
		quant <- 1 - 1/(allLengths[1]+1)
 
	prior_q_distribution[[param]] <- quant # put into the q_distribution
}



for(run in 2:length(allLengths)){

	sd_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, sd)
	mean_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, mean)
	median_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, median)
	zscore_temp <- (realMeans[names(realMeans) %in% names(allRuns)] - mean_temp)/sd_temp
	
	quantile_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, quantile, probs = c(0.025, 0.975))
	
	for(param in names(realMeans)[names(realMeans)%in%names(allRuns)]){
		if(realMeans[param] >= quantile_temp[1, param] && realMeans[param] <= quantile_temp[2, param])
			counts[param] <- counts[param]+1
	}

	#do Cook's p
	for(param in names(realMeans)[names(realMeans)%in%names(allRuns)]){
		qdist <- ecdf(allRuns[(allLengths[run-1]+1):allLengths[run], param]) # find the emperical cdf
		quant <- qdist(realMeans[param]) # find the quantile of the real value

 		if(quant == 0)
			quant <- 1/(allLengths[run]-allLengths[run-1]+1)
		if(quant == 1)
			quant <- 1 - 1/(allLengths[run]-allLengths[run-1]+1)
	
		prior_q_distribution[[param]] <- c(prior_q_distribution[[param]], quant) # put into the q_distribution
	}

	sd_each <- rbind(sd_each, t(sd_temp))
	mean_each <- rbind(mean_each, t(mean_temp))
	median_each <- rbind(median_each, t(median_temp))
	zscore_each <- rbind(zscore_each, zscore_temp)
	quantile_temp <- as.numeric(t(as.vector(quantile_temp)))
	names(quantile_temp) <- names(quantile_each)
	quantile_each <- rbind(quantile_each, quantile_temp)
}

# find the mean of the sds
mean_sd_each <- apply(sd_each, 2, mean)

# find the mean of the zscores
mean_zscore_each <- apply(zscore_each, 2, mean)

cat("mean(sd by run)/sd: \n", mean_sd_each/sd_all, "\n mean(zscore by run): \n", mean_zscore_each, "\n zscore overall: \n", zscore_all, "\n")


# determine the percentbad - percent that the parameters are outside confidence interval
percentgood <- counts/length(allLengths)
percentbad <- 1-percentgood
print(percentbad[1:2])

# determine the log percent that the median is off from the realmeans
percentoff_all <- abs((log(median_each) - log(realMeans[1:2]))/log(realMeans[1:2]))
percentoff <- apply(percentoff_all, 2, mean)

## caterpillar plot
dev.new()
par(mfrow=c(2, 1))
plot(1:length(allLengths), log(median_each[, 1]), xlab = "MCMC Chain Index", ylab = "log parameter interval", main = paste0("log rate of movement (%off (log) ", signif(percentoff[1], 3)," ) (%out ", signif(percentbad[1],2), " )"), pch = 16, ylim = c(-4.25, -2.25))
abline(h = log(realMeans["rateMove"]), col = "green")
arrows(1:length(allLengths), log(quantile_each[, 1]), 1:length(allLengths), log(quantile_each[, 2]), code = 3, angle=90, length=0.01)

plot(1:length(allLengths), log(median_each[, 2]), xlab = "MCMC Chain Index", ylab = "log parameter interval", main = paste0("log weight of jump (%off (log) ", signif(percentoff[2], 3)," ) (%out ", signif(percentbad[2],2), " )"), pch = 16, ylim = log(c(1e-3, 1)))
abline(h = log(realMeans["weightJumpInMove"]), col = "green")
arrows(1:length(allLengths), log(quantile_each[, 3]), 1:length(allLengths), log(quantile_each[, 4]), code = 3, angle=90, length=0.01)

dev.copy2pdf(file = paste0(nameSimul, "_cp4.pdf"))

#convert q_dist to Cook's p
#transform cook's q to a chisq then to a p value:
cookp <- rep(0, length(realMeans)) 
names(cookp) <- names(realMeans)

for(param in names(realMeans)[names(realMeans)%in%names(allRuns)]){

	cookchisq <- sum(qnorm(prior_q_distribution[[param]])^2)
	cookp[param] <- pchisq(cookchisq, df = length(daisyChainSeeds))
}

cookp <- 1-cookp

#plot Cook's q, p
dev.new()
par(mfrow = c(3, 2))
hist(prior_q_distribution[["rateMove"]], xlab = "quantiles", main = paste0("rate move (cook's p ", signif(cookp["rateMove"], 3), ")"))
hist(prior_q_distribution[["weightJumpInMove"]], xlab = "quantiles", main = paste0("weight jump (cook's p ", signif(cookp["weightJumpInMove"], 3), ")"))

get.estimate(prior_q_distribution[["rateMove"]], name = "rate move", xlim = c(0, 1))
abline(h = 1, col = "green")
get.estimate(prior_q_distribution[["weightJumpInMove"]], name = "weight jump", xlim = c(0, 1))
abline(h = 1, col = "green")

rmDens <- density(prior_q_distribution[["rateMove"]], from = 0.01, to = 0.99, kernel = "gaussian", adj = 0.5)
plot(rmDens, xlim = c(0.001, 0.999), main = "rateMove", xlab = "quantiles", ylab = "density")
abline(h = 1, col = "green")
wjDens <- density(prior_q_distribution[["weightJumpInMove"]], from = 0.01, to = 0.99, kernel = "gaussian", adj = 0.5)
plot(wjDens, xlim = c(0.001, 0.999), main = "weight jump", xlab = "quantile", ylab = "density")
abline(h = 1, col = "green")

dev.copy2pdf(file = paste0(nameSimul, "_cookstest.pdf"))



## auxiliary functions
dtrunc <- function(x, spec, a = -Inf, b = Inf, ...)
{
	tt <- rep(0, length(x))
	g <- get(paste("d", spec, sep = ""), mode = "function")
	G <- get(paste("p", spec, sep = ""), mode = "function")
	tt[x>=a & x<=b] <- g(x[x>=a&x<=b], ...)/(G(b, ...) - G(a, ...))
	return(tt)
}

