source("extrapol_field.R")

Monitored<-read.table("completethetasamples_all12000.txt",header=TRUE)

# obtain values from param.r
# delta = rateSkipInMove/rateHopInMove
# rateSkipInMove + rateHopInMove = 1 - rateJumpInMove
# delta * rateHopInMove = rateSkipInMove
# delta * rateHopInMove + rateHopInMove = 1 - rateJumpInMove
# rateHopInMove = (1 - rateJumpInMove) / (1 + delta)
# rateSkipInMove = 1 - rateJumpInMove - (1 - rateJumpInMove)/(1+delta)
# rateSkipInMove = (1 - rateJumpInMove)(1 - 1/(1+delta))
# delta = rateSkipInMove(1-rateJumpInMove-rateSkipInMove)
rateMove <- 0.04
weightHopInMove <- 1
weightJumpInMove <- 0.10
detectRate <- 0.7

true.val<-c(rateMove, weightJumpInMove, detectRate)
names(true.val) <- c("rateMove", "weightJumpInMove", "detectRate")
names(Monitored) <- c("LL", "LP", names(true.val)[1:(dim(Monitored)[2]-2)])
keep <- which(names(true.val) %in% names(Monitored))
true.val <- true.val[keep]
traces(Monitored)

dev.new()
par(mfrow=c(4,1))
plot(Monitored[, 2], xlab = "MCMC iteration", ylab = "Likelihood", pch = ".")
for(param in names(true.val)){
plot(Monitored[[param]], xlab = "MCMC iteration", ylab = param, pch = ".")
}

lims <- data.frame("rateMove" = c(0, .08), "weightJumpInMove" = c(0, 1), "detectRate" = c(0, 1))

dev.new()
par(mfrow=c(4,1))
for(param in names(true.val)){
if(param != "weightJumpInMove")
get.estimate(Monitored[[param]],true.val=true.val[param],name=param, xlim = lims[, param])
}

 
totalWeight <- weightJumpInMove + 1
true.val2<-c(1/totalWeight, weightJumpInMove/totalWeight)
names(true.val2) <- c("rateHop", "rateJump")
xlim = c(0, 1)
totalWeight <- Monitored$weightJumpInMove + 1
rateHop <- 1/totalWeight
rateJump <- Monitored$weightJumpInMove/totalWeight

get.estimate(rateHop, true.val = true.val2["rateHop"], name = "rateHop", xlim = xlim)
get.estimate(rateJump, true.val = true.val2["rateJump"], name = "rateJump", xlim = xlim)

