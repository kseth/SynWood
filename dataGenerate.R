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
	randominitout <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = randominitdays, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = FALSE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro)
	infested2 <- randominitout$infestedDens
	startInfestH2 <- which(infested2 == 1)
	circles <- conc.circles(maps$X, maps$Y, circleRadii, startInfestH2) #new circles centered at new infested
        timeH <- rep(-2, length(startInfestH2))	
 
}

plot_reel(maps$X, maps$Y, infested2, base = 0, top = 1)

# run 1 gillespie simulation to give second timepoint data 
start <- Sys.time()
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH2, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro)
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
secondTimePointStats2 <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestH2, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro)

# obtain stats from the second gillespie simulation now messed up via observation error
if(!is.vector(secondTimePointStats2$statsTable)){
	statsData <- secondTimePointStats2$statsTable[, 1]
}else{
	statsData <- secondTimePointStats2$statsTable
}

# binomial (1 or 0) for all units infested or not
binomialEndInfested <- binomEndInfested2

# close the device so it prints
dev.off()
