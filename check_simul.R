source("spatcontrol/spatcontrol.R", chdir = TRUE)

Monitored<-read.table("completethetasamples_all201000.txt",header=TRUE)
names(Monitored)<-c("LL", "LP", "rateMove", "weightJumpInMove") 
params <- c("rateMove", "weightJumpInMove")
traces(Monitored)

dev.new()
par(mfrow=c(2, 2))
plot(Monitored[, 2], xlab = "MCMC iteration", ylab = "Likelihood", pch = ".")
for(param in params){
plot(Monitored[[param]], xlab = "MCMC iteration", ylab = param, pch = ".")
}

lims <- data.frame("rateMove"=c(0, .08), "weightJumpInMove"=c(0, 1), "weightSkipInMove"=c(0, 1), "detectRate"=c(0, 1))

dev.new()
par(mfrow=c(3, 1))
for(param in params){
	get.estimate(Monitored[[param]], name=param, xlim = lims[, param])
}

weightHop <- 1-Monitored$weightJumpInMove
get.estimate(weightHop, name="weightHopInMove", xlim = c(0, 1))

print(mean(exp(Monitored$LP)))
print(quantile(Monitored$rateMove, c(0.025, 0.5, 0.975)))
print(quantile(1/Monitored$rateMove, c(0.025, 0.5, 0.975))/52)
print(quantile(Monitored$weightJumpInMove, c(0.025, 0.5, 0.975)))
print(quantile((1-Monitored$weightJumpInMove)*Monitored$weightSkipInMove, c(0.025, 0.5, 0.975)))
print(quantile((1-Monitored$weightJumpInMove)*(1-Monitored$weightSkipInMove), c(0.025, 0.5, 0.975)))


