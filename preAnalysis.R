library(testthat)
library(verification)
library(msm)
source("../spatcontrol/spatcontrol.R", chdir=TRUE) #has all the spatial analysis + convergence functions
source("logLik.R") # override some sl functions and add synLik
source("functions_migration.R")

# set spam memory options
spam.options(nearestdistnnz=c(13764100,400))

# set seed
seedSimul <- 100000

#======================
# Set details of grid (simulation environment)
#====================== 
# size of grid
num.rows <- 48 
num.cols <- 48
row.dist <- 10
 
## make a grid map 
maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)

# parameters for uniform hop/skip/jump model
limitHopSkip <- 30
limitJump <- 100
lowerLimitJump <- 30 

# set blockIndex to NULL
# no blocks!
blockIndex = NULL

# make stratified matrix (no skips, set blockIndex to NULL)
# matrix contains where each house can hop to, where can jump to
stratHopSkipJump <- generate_stratified_mat(coords=maps[, c("X", "Y")], limitHopSkip, limitJump, lowerLimitJump=lowerLimitJump, blockIndex=blockIndex)

#=======================
# Pick spatial locations to seed the starting infestation
#=======================
# pick starting point for simulations
startInfestH <- ceiling(num.rows*(num.rows/2) + num.rows/2)
startInfestH <- c(startInfestH + 1, startInfestH - 1) 

#=======================
# Pick length of simulation
#=======================
# let the infestation spread for two years, nbit <- 52 * 2
nbit <- 104

# dummy thing to keep track of times of infestation
timeH <- rep(-2, length(startInfestH))

#=======================
# Priors
#=======================
detectRate <- 1 # true detection rate (1 == don't remove data) 
sampleDR <- FALSE # if true, MCMC will sample over error rates
defaultDR <- 1 # DR assumed by multiGilStat (1 if detectRate==1, can set to 0.7, only used if !sampleDR)

realMeans <- c(0.04, 0.10, 0.80)
priorMeans<-c(0.04, 0.10, 0.80) 
priorSd <- c(1, 0.4, 0.20)
priorType <- c("lnorm", "boundednorm", "boundednorm")
priorIntervals <- list(NULL, c(0, 1), c(0, 1)) # only considered if priorType is bounded

names(priorMeans)<-c("rateMove" , "weightJumpInMove", "detectRate") 
names(priorSd)<-names(priorMeans)
names(priorType)<-names(priorMeans)
names(priorIntervals)<-names(priorMeans)
names(realMeans) <- names(priorMeans)

if(!sampleDR){ #if don't want to sample over detectRate
	priorMeans<-priorMeans[-3]
	realMeans<-realMeans[-3]
	priorSd<-priorSd[-3]
	priorType<-priorType[-3]
	priorIntervals <- priorIntervals[-3]
}

#=========================
# Set fitting method
#=========================
# which statistics to use? 
# choices: "grid", "circles", "semivariance"
useStats <- c("grid", "circles", "semivariance") # disregarded if useBinLik

# distance classes for the general variogram
genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30))
# genIntervals <- seq(10, 40, 15)  # if want to combine with grid

# partition sizes for the map
# i.e. 12x12 cells (each cell 3x3), 8x8 cells (each cell 6x6), ..., 3x3 cells (each cell 16x16)
partitionSizes <- c(16, 12, 8, 6, 4, 3)

# radii for the concentric circles
# first interval 0-20, 20-35, etc.
circleRadii <- c(0, 20, 35, 50, 80, 110, 155, 200)

#=========================
# Construct the things needed for semivariance, grid, circles
#=========================
# bin the map into different distances classes
bin_dist_out <- makeDistClasses(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals)

# partition the map
map.partitions <- list()
length(map.partitions) <- length(partitionSizes) #number different grid partitions will be used

for(part in 1:length(partitionSizes))
	map.partitions[[part]] <- partitionMap(maps$X, maps$Y, partitionSizes[part]) 

# make the concentric circles
circles <- conc.circles(maps$X, maps$Y, circleRadii, startInfestH) 

#========================
# 10,000 repetitions from real values
#========================
real_reps <- 10000

ts <- Sys.time()
real_sims <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = realMeans["rateMove"], weightSkipInMove = 0, weightJumpInMove = realMeans["weightJumpInMove"], Nrep = real_reps, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats=TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)
real_sims_stats <- real_sims$statsTable
te <- Sys.time()

print(te-ts)
stop()

#========================
# Num draws from prior, num simulations
#========================
num_draws <- 10000
sim_results <- data.frame()

for(draw in 1:num_draws){
	print(draw)
	sim_vals <- rep(0, length(priorMeans))
	names(sim_vals) <- names(priorMeans)
	for(name in names(sim_vals)){
		if(priorType[name] == "lnorm")
			sim_vals[name] <- rlnorm(1, meanlog = log(priorMeans[name]), sdlog = priorSd[name])
		else if(priorType[name] == "norm")
			sim_vals[name] <- rnorm(1, mean = priorMeans[name], sd = priorSd[name])
		else if(priorType[name] == "boundednorm")
			sim_vals[name] <- rtnorm(1, mean = priorMeans[name], sd = priorSd[name], lower = priorIntervals[[name]][1], upper = priorIntervals[[name]][2])
		else sim_vals[name] <- runif(1, priorIntervals[[name]][1], priorIntervals[[name]][2])
	}

	simulation <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump,
					    blockIndex = blockIndex, 
			    infestH = startInfestH, timeH = timeH, endTime = nbit, 
			    rateMove = sim_vals["rateMove"],
			    weightSkipInMove = 0,
			    weightJumpInMove = sim_vals["weightJumpInMove"], 
			    Nrep = 1, 
			    coords = maps[, c("X", "Y")], 
			    breaksGenVar = genIntervals,
			    seed=seedSimul,
			    simul=TRUE,
			    getStats=TRUE,
			    dist_out=bin_dist_out, map.partitions=map.partitions, conc.circs=circles,
			    typeStat=useStats,
			    detectRate=ifelse("detectRate" %in% names(sim_vals), sim_vals["detectRate"], defaultDR))

	stats <- simulation$statsTable
	ll <- synLik(real_sims_stats, stats, trans = NULL)
	sim_vals <- c(sim_vals, ll)
	sim_results <- rbind(sim_results, sim_vals)
}

library(rgl)
plot3d(sim_results[, 1], sim_results[, 2], sim_results[, 3], xlim = c(0, 0.1), ylim = c(0, 0.5))

goodLL <- which(sim_results[, 3] > 0)
plot3d(sim_results[goodLL, 1], sim_results[goodLL, 2], sim_results[goodLL, 3])
