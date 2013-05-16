library(vioplot)
source("spatcontrol/spatcontrol.R", chdir = TRUE)

Monitored<-read.table("thetasamples_all314159.txt",header=TRUE)
Monitored<-Monitored[-(1:300), ]
params <- c("rateMove", "weightSkipInMove", "weightJumpInMove")
traces(Monitored)

dev.new()
par(mfrow=c(2,2))
plot(Monitored[, 2], xlab = "MCMC iteration", ylab = "Likelihood", pch = ".")
for(param in params){
plot(Monitored[[param]], xlab = "MCMC iteration", ylab = param, pch = ".")
}

lims <- data.frame("rateMove" = c(0, .08), "weightJumpInMove" = c(0, 1), "weightSkipInMove" = c(0, 1), "detectRate" = c(0, 1))

dev.new()
par(mfrow=c(3, 1))
for(param in params){
	get.estimate(Monitored[[param]], name=param, xlim = lims[, param])
}

weightLocal <- 1-Monitored$weightJumpInMove
weightSkip <- weightLocal * Monitored$weightSkipInMove
weightHop <- weightLocal-weightSkip

dev.new()
par(mfrow=c(3, 1))
get.estimate(weightHop, name="percent Hops", xlim = c(0, 1))
get.estimate(weightSkip, name="percent Skips", xlim = c(0, 1))
get.estimate(Monitored$weightJumpInMove, name="percent Jumps", xlim = c(0, 1))

dev.new()
vioplot(weightHop, weightSkip, Monitored$weightJumpInMove, names = c("Hops", "Skips", "Jumps"), ylim = c(0, 1), col = "violet")
title(xlab = "Move Type", ylab = "Percent of Moves")
