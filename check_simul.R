source("spatcontrol/spatcontrol.R", chdir = TRUE)

Monitored<-read.table("thetasamples_all109000.txt",header=TRUE)
names(Monitored) <- c("LL", "LP","rateMove", "rateJump", "rateIntro")
traces(Monitored)

dev.new()
par(mfrow=c(5,1))
plot(Monitored[, 2], xlab = "MCMC iteration", ylab = "Likelihood", pch = ".")
for(param in names(Monitored)){
plot(Monitored[[param]], xlab = "MCMC iteration", ylab = param, pch = ".")
}

lims <- data.frame("rateMove" = c(0, .04), "rateJump" = c(0, 1), "rateIntro" = c(0, 0.20))

dev.new()
par(mfrow=c(3,1))
for(param in names(Monitored)[3:5]){
	get.estimate(Monitored[[param]], name=param, xlim = lims[, param])
}
