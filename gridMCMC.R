#library(testthat)
#library(verification)
#source("RanalysisFunctions.R")
#source("logLik.r") # override some sl functions and add synLik
#source("functions_sampling.R")
#source("convergence_test.r")
#source("extrapol_field.R")
#source("functions_migration.R")
#source("models.R")
#source("MCMC.R")

#==================
# set parameters for simulation
#==================
## name the simulation!
## nameSimul <- "GRID_36x36_Hop_Jump_BinLik_9914_Messy_FitNoise"

## set the seed for the simulation
## seedSimul <- 9914 
set.seed(seedSimul)

## name the MCMC file output
## monitor.file <- "thetasamples_all.txt"

## set spam memory options
## spam.options(nearestdistnnz=c(13764100,400))

## how many gillespie repetitions per iteration
## Nrep <- 800 

## Make simulation messy or not messy
## detectRate <- 0.7 # true detection rate 
## sampleDR <- TRUE # if true, MCMC will sample over error rates
## defaultDR <- 1 # DR assumed by multiGilStat (should be 1 if detectRate==1, can set to 0.7, only used if not sampling over DR)
 
## size of grid
## num.rows <- 36
## num.cols <- 36
## row.dist <- 10
 
## parameters for uniform hop/skip/jump model
## limitHopSkip <- 40
## limitJump <- 200
## lowerLimitJump <- 100 
## rateMove <- 0.04

## the noKernelMultiGilStat normalizes these weights
## weightHopInMove <- 1	# must always be 1!
## weightSkipInMove <- 0.0
## weightJumpInMove <- 0.1 

## which likelihood to use? 
## useBinLik <- FALSE

## which statistics to use?
## useStats <- c("grid", "circles") # disregarded if useBinLik == TRUE

## make a map with just x, y
## maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)

## distance classes for the general variogram
## genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30))
## genIntervals <- seq(10, 40, 15)  if want to combine with grid

## bin the map into different distances classes
## bin_dist_out <- makeDistClasses(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals)

## partition the map
## map.partitions <- list()
## length(map.partitions) <- 6 #6 different grid partitions will be used

## map.partitions[[1]] <- partitionMap(maps$X, maps$Y, 12) #into 12 by 12 (each cell 3 by 3)
## map.partitions[[2]] <- partitionMap(maps$X, maps$Y, 9)  #into 9 by 9 (each cell 4 by 4)
## map.partitions[[3]] <- partitionMap(maps$X, maps$Y, 6)  #into 6 by 6 (each cell 6 by 6)
## map.partitions[[4]] <- partitionMap(maps$X, maps$Y, 4)  #into 4 by 4 (each cell 9 by 9)
## map.partitions[[5]] <- partitionMap(maps$X, maps$Y, 3)  #into 3 by 3 (each cell 12 by 12)
## map.partitions[[6]] <- partitionMap(maps$X, maps$Y, 2)  #into 2 by 2 (each cell 18 by 18)

## set blockIndex to NULL
## no blocks!
## blockIndex = NULL

## make stratified matrix (no skips, set blockIndex to NULL)
## stratHopSkipJump <- generate_stratified_mat(coords=maps[, c("X", "Y")], limitHopSkip, limitJump, lowerLimitJump=lowerLimitJump, blockIndex=blockIndex)

## pick starting point for simulations
## startInfestH <- ceiling(num.rows*(num.rows/2) + num.rows/2)
## startInfestH <- c(startInfestH, startInfestH + 1, startInfestH - 1) 

## make the concentric circles
## circleRadii <- c(0, 20, 35, 50, 80, 110, 155, 200)
## circles <- conc.circles(maps$X, maps$Y, circleRadii, startInfestH) 

## create a dummy timeH
## timeH <- rep(-2, length(startInfestH))

## let the infestation spread for two years, nbit <- 52 * 2
## nbit <- 104

#===================
# Run simulation and get end point statistics 
#===================

## plot initially infested houses
dev.new()
par(mfrow = c(2, 2))
infested <- rep(0, length(maps$X))
infested[startInfestH] <- 1
plot_reel(maps$X, maps$Y, infested, base = 0, top = 1)

## run 1 gillespie simulation to give second timepoint data 
start <- Sys.time()
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = rateMove, weightHopInMove = weightHopInMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)
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
secondTimePointStatsR <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestHR, timeH=timeH, endTime = nbit, rateMove = rateMove, weightHopInMove = weightHopInMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)

## do the same thing in C	
## endInfestHC <- which(binomEndInfested == 1)

## calculate statistics
## secondTimePointStatsC <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestHC, timeH=timeH, endTime = nbit, rateMove = rateMove, weightHopInMove = weightHopInMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, detectRate = detectRate)

#plot the results
## plot_reel(maps$X, maps$Y, secondTimePointStatsC$infestedDens, base = 0, top = 1)

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
priorMeans<-c(0.045, 0.05, 0.80)
priorSd <- c(1, 1, 0.20)
priorType <- c("boundedlnorm", "boundedlnorm", "boundednorm")
priorIntervals <- list(c(0, 1), c(0, 10), c(0, 1))
realMeans<-c(rateMove, weightJumpInMove, detectRate)
sampling <- c("lnorm", "lnorm", "boundednorm")
sdProposal <- c(0.4, 0.4, 0.2)
names(priorMeans)<-c("rateMove" , "weightJumpInMove", "detectRate") 
names(sampling)<-names(priorMeans)
names(realMeans)<-names(priorMeans)
names(priorSd)<-names(priorMeans)
names(priorType)<-names(priorMeans)
names(priorIntervals) <- names(priorMeans) 

if(!sampleDR){ #if don't want to sample over detectRate
	priorMeans<-priorMeans[-3]
	priorSd<-priorSd[-3]
	priorType<-priorType[-3]
	priorIntervals <- priorIntervals[-3]
	realMeans<-realMeans[-3]
	sampling<-sampling[-3]
	sdProposal<-sdProposal[-3]
}

#=================
# List of data to pass to model + sampler
#=================

MyDataFullSample <- list(y={if(useBinLik) binomEndInfestedR else statsData},
	     trans=NULL,
	     stratHopSkipJump = stratHopSkipJump,
	     blockIndex=blockIndex,
	     dist_out = {if("semivariance" %in% useStats) bin_dist_out else NULL},
	     map.partitions = {if("grid" %in% useStats) map.partitions else NULL}, 
	     conc.circs = {if("circles" %in% useStats) circles else NULL}, 
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

