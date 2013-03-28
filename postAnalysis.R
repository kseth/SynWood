source("extrapol_field.R")

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

## determine mean(sd(eachRun))/sd(allRuns)
## determine mean(zscore(eachRun))/zscore(allRuns)
sd_all <- apply(allRuns[, -(1:2)], 2, sd)
mean_all <- apply(allRuns[, -(1:2)], 2, mean)
zscore_all <- (realMeans[names(realMeans) %in% names(allRuns)] - mean_all)/sd_all

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

counts <- c(0, 0, 0)
names(counts) <- names(realMeans)

for(param in names(realMeans)[names(realMeans)%in%names(allRuns)])
if(realMeans[param] >= quantile_temp[1, param] && realMeans[param] <= quantile_temp[2, param])
counts[param] <- counts[param]+1

for(run in 2:length(allLengths)){

sd_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, sd)
mean_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, mean)
median_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, median)
zscore_temp <- (realMeans[names(realMeans) %in% names(allRuns)] - mean_temp)/sd_temp

quantile_temp <- apply(allRuns[(allLengths[run-1]+1):allLengths[run], -(1:2)], 2, quantile, probs = c(0.025, 0.975))

for(param in names(realMeans)[names(realMeans)%in%names(allRuns)])
if(realMeans[param] >= quantile_temp[1, param] && realMeans[param] <= quantile_temp[2, param])
counts[param] <- counts[param]+1


sd_each <- rbind(sd_each, t(sd_temp))
mean_each <- rbind(mean_each, t(mean_temp))
median_each <- rbind(median_each, t(median_temp))
zscore_each <- rbind(zscore_each, zscore_temp)
quantile_temp <- as.numeric(t(as.vector(quantile_temp)))
names(quantile_temp) <- names(quantile_each)
quantile_each <- rbind(quantile_each, quantile_temp)
}

mean_sd_each <- apply(sd_each, 2, mean)
mean_zscore_each <- apply(zscore_each, 2, mean)
cat("mean(sd by run)/sd: \n", mean_sd_each/sd_all, "\n mean(zscore by run): \n", mean_zscore_each, "\n zscore overall: \n", zscore_all, "\n")

percentgood <- counts/length(allLengths)
percentbad <- 1-percentgood
print(percentbad[1:2])

percentoff_all <- abs((median_each - realMeans[1:2])/realMeans[1:2])
percentoff <- apply(percentoff_all, 2, mean)

graphics.off()

## caterpillar plot
dev.new()
par(mfrow=c(2, 1))
plot(1:length(allLengths), median_each[, 1], xlab = "MCMC Chain Index", ylab = "parameter interval", main = paste0("rate of movement (%off ", signif(percentoff[1], 3)," ) (%out ", signif(percentbad[1],2), " )"), pch = 16, ylim = c(0, 0.1))
abline(h = realMeans["rateMove"], col = "green")
arrows(1:length(allLengths), quantile_each[, 1], 1:length(allLengths), quantile_each[, 2],code = 3, angle=90, length=0.01)

plot(1:length(allLengths), log(median_each[, 2]), xlab = "MCMC Chain Index", ylab = "log parameter interval", main = paste0("log weight of jump (%off ", signif(percentoff[2], 3)," ) (%out ", signif(percentbad[2],2), " )"), pch = 16, ylim = log(c(1e-3, 1)))
abline(h = log(realMeans["weightJumpInMove"]), col = "green")
arrows(1:length(allLengths), log(quantile_each[, 3]), 1:length(allLengths), log(quantile_each[, 4]), code = 3, angle=90, length=0.01)

dev.copy2pdf(file = "caterpillar3.pdf")

