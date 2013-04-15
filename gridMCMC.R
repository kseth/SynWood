#===================
# Run simulation 
# Generate data 
#===================
# set the seed for the simulation
set.seed(seedSimul)

# create a dummy timeH
timeH <- rep(-2, length(startInfestH))

# plot initially infested houses
dev.new()
par(mfrow = c(2, 2))
infested <- rep(0, length(maps$X))
infested[startInfestH] <- 1
plot_reel(maps$X, maps$Y, infested, base = 0, top = 1)

# use randominitdays to generate starting point
if(randominitdays == 0){ # no random init, seeding points are starting points
	startInfestH2 <- startInfestH
}else{ # run gillespie to generate a random starting set from seeding points
	randominitout <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = randominitdays, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = FALSE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)
	infested2 <- randominitout$infestedDens
	startInfestH2 <- which(infested2 == 1)
}

plot_reel(maps$X, maps$Y, infested2, base = 0, top = 1)

# run 1 gillespie simulation to give second timepoint data 
start <- Sys.time()
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH2, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)
print(Sys.time() - start)

# plot results of gillespie
binomEndInfested <- secondTimePointSimul$infestedDens
cat("starting # infested:", length(startInfestH2), " ending # infested:", length(which(binomEndInfested!=0)), "\n")
plot_reel(maps$X, maps$Y, binomEndInfested, base = 0, top = 1)

# remove data
binomEndInfested2 <- simulObserved(binomEndInfested, detectRate, 1)
endInfestH2 <- which(binomEndInfested2 == 1)

# plot the results
plot_reel(maps$X, maps$Y, binomEndInfested2, base = 0, top = 1)

# calculate statistics
secondTimePointStats2 <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestH2, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)

# obtain stats from the second gillespie simulation now messed up via observation error
if(!is.vector(secondTimePointStats2$statsTable)){
	statsData <- secondTimePointStats2$statsTable[, 1]
}else{
	statsData <- secondTimePointStats2$statsTable
}

# close the device so it prints
dev.off()

#=================
# Model and
# List of data to pass to model + sampler
#=================
modelToUse<-{if(useBinLik) binomNoKernelModel else noKernelModel}

MyDataFullSample <- list(y={if(useBinLik) binomEndInfested2 else statsData},
	     trans=NULL,
	     stratHopSkipJump = stratHopSkipJump,
	     blockIndex=blockIndex,
	     dist_out = {if(!useBinLik && ("semivariance" %in% useStats)) bin_dist_out else NULL},
	     map.partitions = {if(!useBinLik && ("grid" %in% useStats)) map.partitions else NULL}, 
	     conc.circs = {if(!useBinLik && ("circles" %in% useStats)) circles else NULL}, 
	     useStats = useStats,
	     infestH=startInfestH2,
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
	     sampling=sampling
		)

#=================
# Test modelToUse to make sure something meaningful comes out
#=================
start<-Sys.time()
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
# Make call to MCMC
#=================
MCMC(MyDataFullSample, Model=modelToUse, sdprop=sdProposal, monitor.file=monitor.file)

