#===================
# Run simulation 
# Generate data for a second time point from startingInfested
# possibly with observation error
#===================
# set the seed for the simulation
set.seed(seedSimul)

# create a dummy timeH
timeH <- rep(-2, length(startingInfested))

# plot initially infested houses
dev.new()
par(mfrow = c(2, 2))
infested <- rep(0, length(maps$X))
infested[startingInfested] <- 1
plot_reel(maps$X, maps$Y, infested, base = 0, top = 1)

# run 1 gillespie simulation to give second timepoint data 
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startingInfested, timeH=timeH, endTime = nbit, rateMove = rateMove, rateHopInMove = weightHop, rateSkipInMove = weightSkip, rateJumpInMove = weightJump, Nrep = 1, coords = maps[, c("X", "Y")], simul=TRUE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro,iPartLMoments=iPartLMoments)


# plot results of gillespie
binomEndInfested_prenoise <- secondTimePointSimul$infestedDens
cat("starting # infested:", length(startingInfested), " ending # infested:", length(which(binomEndInfested_prenoise!=0)), "\n")
plot_reel(maps$X, maps$Y, binomEndInfested_prenoise, base = 0, top = 1)

# remove data
binomEndInfested_noisy <- simulObserved(binomEndInfested_prenoise, detectRate, 1)
endingInfested_noisy <- which(binomEndInfested_noisy == 1)

# plot the results
plot_reel(maps$X, maps$Y, binomEndInfested_noisy, base = 0, top = 1)

# binomial (1 or 0) for all units infested or not at the end
binomialEndInfested <- binomEndInfested_noisy

#final endInfestedHouses
endInfestedHouses <- endingInfested_noisy

# calculate statistics on possibly noisy presence/absence map
secondTimePointStats2 <- noKernelMultiGilStat(stratHopSkipJump=stratHopSkipJump, blockIndex=blockIndex, infestH =endInfestedHouses, timeH=timeH, endTime=nbit,  coords = maps[, c("X", "Y")], simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro, whichPairwise=whichPairwise,iPartLMoments=iPartLMoments)
	
# obtain stats from the second gillespie simulation now messed up via observation error
if(!is.vector(secondTimePointStats2$statsTable)){
	statsData <- secondTimePointStats2$statsTable[, 1]
}else{
	statsData <- secondTimePointStats2$statsTable
}
