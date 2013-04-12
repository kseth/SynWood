## set the seed for the simulation
set.seed(seedSimul)

#===================
# Run simulation 
# Generate data 
#===================

## plot initially infested houses
dev.new()
par(mfrow = c(2, 2))
infested <- rep(0, length(maps$X))
infested[startInfestH] <- 1
plot_reel(maps$X, maps$Y, infested, base = 0, top = 1)

## run 1 gillespie simulation to give second timepoint data 
start <- Sys.time()
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)
print(Sys.time() - start)

## plot results of gillespie
binomEndInfested <- secondTimePointSimul$infestedDens
cat("starting # infested:", length(startInfestH), " ending # infested:", length(which(binomEndInfested!=0)), "\n")
plot_reel(maps$X, maps$Y, binomEndInfested, base = 0, top = 1)

## remove data
binomEndInfestedR <- simulObserved(binomEndInfested, detectRate, 1)
endInfestHR <- which(binomEndInfestedR == 1)

#plot the results
plot_reel(maps$X, maps$Y, binomEndInfestedR, base = 0, top = 1)

#calculate statistics
secondTimePointStatsR <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestHR, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)

## obtain stats from the second gillespie simulation now messed up via observation error
if(!is.vector(secondTimePointStatsR$statsTable)){
	statsData <- secondTimePointStatsR$statsTable[, 1]
}else{
	statsData <- secondTimePointStatsR$statsTable
}

## close the device so it prints
dev.off()

#==================
# Priors (also the place to change the parameters)
#==================
priorMeans<-c(0.04, 0.05, 0.80)
priorSd <- c(1, 0.5, 0.20)
priorType <- c("lnorm", "noPrior", "boundednorm")
priorIntervals <- list(c(0, 1), c(0, 10), c(0, 1)) # only considered if bounded function
realMeans<-c(rateMove, weightJumpInMove, detectRate)
sampling <- c("lnorm", "boundednorm", "boundednorm")
sdProposal <- c(0.4, 0.2, 0.2)
names(priorMeans)<-c("rateMove" , "weightJumpInMove", "detectRate") 
names(sampling)<-names(priorMeans)
names(realMeans)<-names(priorMeans)
names(priorSd)<-names(priorMeans)
names(priorType)<-names(priorMeans)
names(priorIntervals) <- names(priorMeans) 

initValues<-priorMeans #the values with which to initialize the sampler
initValues["rateMove"]<-0.05
initValues["weightJumpInMove"]<-0.15
initValues["detectRate"]<-0.70

if(!sampleDR){ #if don't want to sample over detectRate
	priorMeans<-priorMeans[-3]
	priorSd<-priorSd[-3]
	priorType<-priorType[-3]
	priorIntervals <- priorIntervals[-3]
	realMeans<-realMeans[-3]
	sampling<-sampling[-3]
	sdProposal<-sdProposal[-3]
	initValues<-initValues[-3]
}


#=================
# List of data to pass to model + sampler
#=================

MyDataFullSample <- list(y={if(useBinLik) binomEndInfestedR else statsData},
	     trans=NULL,
	     stratHopSkipJump = stratHopSkipJump,
	     blockIndex=blockIndex,
	     dist_out = {if(!useBinLik && ("semivariance" %in% useStats)) bin_dist_out else NULL},
	     map.partitions = {if(!useBinLik && ("grid" %in% useStats)) map.partitions else NULL}, 
	     conc.circs = {if(!useBinLik && ("circles" %in% useStats)) circles else NULL}, 
	     useStats = useStats,
	     infestH=startInfestH,
	     timeH=timeH,
	     endTime=nbit,
	     maps=maps,
	     nbit=nbit,
	     Nrep=Nrep,
	     priorMeans=priorMeans,
	     priorSd=priorSd,
	     priorType=priorType,
	     priorIntervals=priorIntervals,
	     initValues=initValues,
	     defaultDR=defaultDR,
	     genIntervals=genIntervals,
	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
	     sampling=sampling # method of sampling parameters
		)

#=================
## Test modelToUse to make sure something meaningful comes out
#=================
start<-Sys.time()
modelToUse<-{if(useBinLik) binomNoKernelModel else noKernelModel}
ModelOutGood<-modelToUse(priorMeans,MyDataFullSample)
cat(Sys.time()-start, "\n")
start<-Sys.time()
ModelOutBest<-modelToUse(realMeans,MyDataFullSample)
cat(Sys.time()-start, "\n")

# weibull order plotting to check if circle stats are normal
# for(i in 1:dim(secondTimePointSimul$circle.statsTable)[1])
# {
# 
# order <- order(secondTimePointSimul$circle.statsTable[i, ])
# plot((1:dim(secondTimePointSimul$circle.statsTable)[2])/(dim(secondTimePointSimul$circle.statsTable)[2]+1), secondTimePointSimul$circle.statsTable[i, order], main = i, type = "l")
# readline()
#  
# }
# 
# weibull order plotting to check if grid stats are normal
# for(i in 1:dim(secondTimePointSimul$grid.statsTable)[1])
# {
# 
# order <- order(secondTimePointSimul$grid.statsTable[i, ])
# plot((1:dim(secondTimePointSimul$grid.statsTable)[2])/(dim(secondTimePointSimul$grid.statsTable)[2]+1), secondTimePointSimul$grid.statsTable[i, order], main = i, type = "l")
# readline()
#  
# }

# good should be worse than best (ideally, need -4 because simulations may not be ideal)
# expect_true(ModelOutGood$Dev>ModelOutBest$Dev-4)
#=================
## Make call to MCMC
#=================
MCMC(MyDataFullSample, Model=modelToUse, sdprop=sdProposal, monitor.file=monitor.file)

