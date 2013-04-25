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

randominitdays <- 52

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

realMeans <- c(0.04, 0.1, 0.80)
priorMeans<-c(0.04, 0.1, 0.80) #should be realMeans (unless want some sort of skewed prior) 
priorSd <- c(0.5, 0.3, 0.20)
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
useStats <- c("grid", "circles", "semivariance")

# distance classes for the general variogram
# genIntervals <- c(seq(10, 190, 30), seq(240, 340, 50))
genIntervals <- seq(0,72,9)
genIntervals <- c(0, cumsum(genIntervals) + 15)

# partition sizes for the map
# i.e. 12x12 cells (each cell 3x3), 8x8 cells (each cell 6x6), ..., 3x3 cells (each cell 16x16)
partitionSizes <- c(16, 12, 8, 6, 4, 3)

# radii for the concentric circles
# first interval 0-20, 20-35, etc.
# circleRadii <- c(0, 20, 35, 50, 80, 110, 155, 200)
circleRadii <- genIntervals

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
# get startInfestH
#========================

# use randominitdays to generate starting point
if(randominitdays == 0){ # no random init, seeding points are starting points
	startInfestH <- startInfestH
}else{ # run gillespie to generate a random starting set from seeding points
	randominitout <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = randominitdays, rateMove = realMeans["rateMove"], weightSkipInMove = 0, weightJumpInMove = realMeans["weightJumpInMove"], Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = FALSE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)
	startInfestH <- which(randominitout$infestedDens == 1) 
	circles <- conc.circles(maps$X, maps$Y, circleRadii, startInfestH) 
}

#========================
# 10,000 repetitions from real values
#========================
real_reps <- 10000

ts <- Sys.time()

real_sims <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = realMeans["rateMove"], weightSkipInMove = 0, weightJumpInMove = realMeans["weightJumpInMove"], Nrep = real_reps, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats=TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)

te <- Sys.time()

print(te-ts)

real_sims_stats <- real_sims$statsTable
rownames(real_sims_stats) <- paste0("stats",1:dim(real_sims_stats)[1])

#========================
# Simulations and their statistics
#========================
num_draws <- 10000 #number iterations (draws from prior)
sim_params <- data.frame() #data frame to hold simulated thetas
sim_stats <- data.frame() #data frame to hold stats from simulated thetas

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

	sim_params <- rbind(sim_params, sim_vals)
	sim_stats <- rbind(sim_stats, stats)
}

names(sim_params) <- c("rateMove", "rateJump") #assign names
sim_stats <- as.matrix(sim_stats) #convert dataframe to matrix
colnames(sim_stats) <- paste0("stats", 1:dim(sim_stats)[2]) #rename columns

# assign which statistics are which (may have to be done manually if change statistics)
grid_stats <- paste0("stats", 2:(length(partitionSizes)*2)-1)
grid_var_stats <- paste0("stats", c((1:5) * 2, 11))
grid_count_stats <- paste0("stats", (1:5) * 2 - 1)
circ_stats <- paste0("stats", length(partitionSizes)*2-1+1:(length(circleRadii)*2-2))
semivar_stats <- paste0("stats", length(partitionSizes)*2-1+length(circleRadii)*2-2+1:(length(genIntervals)*2-2))
num_inf_stats <- paste0("stats", dim(sim_stats)[2])

#========================
# Calculating likelihoods
#========================
sim_ll <- data.frame()
for(draw in 1:num_draws){

	print(draw)

	ll <- rep(0, 8)
	ll[1] <- synLik(real_sims_stats, sim_stats[draw, c(grid_stats, circ_stats, semivar_stats, num_inf_stats)], trans = NULL)
	ll[2] <- synLik(real_sims_stats[grid_stats, ], sim_stats[draw, grid_stats], trans = NULL)
	ll[3] <- synLik(real_sims_stats[circ_stats, ], sim_stats[draw, circ_stats], trans = NULL)
	ll[4] <- synLik(real_sims_stats[semivar_stats, ], sim_stats[draw, semivar_stats], trans = NULL)
	ll[5] <- log(density(real_sims_stats[num_inf_stats, ], from=sim_stats[draw, num_inf_stats], to=sim_stats[draw, num_inf_stats], n=1)$y)
	ll[6] <- dnorm(sim_stats[draw, num_inf_stats], mean=mean(real_sims_stats[num_inf_stats, ]), sd=sd(real_sims_stats[num_inf_stats, ]), log=TRUE)
	ll[7] <- synLik(real_sims_stats[grid_var_stats, ], sim_stats[draw, grid_var_stats], trans = NULL)
	ll[8] <- synLik(real_sims_stats[grid_count_stats, ], sim_stats[draw, grid_count_stats], trans = NULL)
	sim_ll <- rbind(sim_ll, ll)
}

names(sim_ll) <- c("all_stats", "grid_stats", "circ_stats", "semivar_stats", "num_inf_stats", "normal_num_inf_stats", "grid_var_stats", "grid_count_stats")

stop()

#==================
# Plotting
#==================

## determine which quantile to cut off the likelihoods
cutoff <- c(1, 0.9, 0.8, 0.7, 0.60, 0.50)
q_ll <- apply(sim_ll, 2, quantile, cutoff)
print(q_ll)

## 60% cut off means keep the top 2/5ths of results
cutoff_ll <- q_ll["60%", ]

## 3D plotting
library(rgl)

goodLL <- which(sim_ll[, "all_stats"] > cutoff_ll["all_stats"])
plot3d(sim_params[goodLL, 1], sim_params[goodLL, 2], sim_ll[goodLL, "all_stats"], col = xtocolors(sim_ll[goodLL, "all_stats"]))

## 2D plotting
#uniform xlims, ylims
xlim = c(0.01, 0.06)
ylim = c(0, 1)

dev.new()
par(mfrow = c(2, 4))

for(stat in names(sim_ll)){
	goodLL <- which(sim_ll[, stat] > cutoff_ll[stat])

	plot(sim_params[goodLL, 1], sim_params[goodLL, 2], col = xtocolors(sim_ll[goodLL, stat], crp=colorRampPalette(c("black", "red", "orange", "yellow"))), xlab = "rateMove", ylab = "rateJump", main = stat, xlim = xlim, ylim = ylim)
	abline(h = realMeans["weightJumpInMove"])
	abline(v = realMeans["rateMove"])
}

## 2D contour plotting
dev.new()
par(mfrow = c(2, 4))

for(stat in names(sim_ll)){
	goodLL <- which(sim_ll[, stat] > cutoff_ll[stat])

	plot(sim_params[goodLL, 1], sim_params[goodLL, 2], col = xtocolors(sim_ll[goodLL, stat], crp=colorRampPalette(c("black", "red", "orange", "yellow"))), xlab = "rateMove", ylab = "rateJump", main = stat, xlim = xlim, ylim = ylim, pch = 1)
	abline(h = realMeans["weightJumpInMove"])
	abline(v = realMeans["rateMove"])

	gridfromsample <- grid.from.sample(sim_params[goodLL, 1], sim_params[goodLL, 2], sim_ll[goodLL, stat], steps = 20, tr = 0.01)
	contour(gridfromsample$xs, gridfromsample$ys, gridfromsample$zs, add = TRUE)


}

#======================
# Finding credible intervals
#======================
library(pracma)

##ratemove crI
order_rm <- order(sim_params$rateMove)
rateMoves <- sim_params$rateMove[order_rm]
density_rm <- exp(sim_ll$all_stats)
density_rm <- density_rm[order_rm]
area <- trapz(rateMoves, density_rm)
density_rm <- density_rm/area

##rateJump crI
order_rj <- order(sim_params$rateJump)
rateJumps <- sim_params$rateJump[order_rj]
density_rj <- exp(sim_ll$all_stats)
density_rj <- density_rj[order_rj]
area <- trapz(rateJumps, density_rj)
density_rj <- density_rj/area

dev.new()
par(mfrow = c(2, 1))
plot(rateMoves, density_rm, type = "l")
abline(v = realMeans["rateMove"], col = "green")
plot(rateJumps, density_rj, type = "l")
abline(v = realMeans["weightJumpInMove"], col = "green")

rm_post_stats <- data.frame()
rj_post_stats <- data.frame()

for(stat in names(sim_ll)){

	##ratemove crI
	order_rm <- order(sim_params$rateMove)
	rateMoves <- sim_params$rateMove[order_rm]
	density_rm <- exp(sim_ll[, stat])
	density_rm <- density_rm[order_rm]
	area <- trapz(rateMoves, density_rm)
	density_rm <- density_rm/area
	
	##rateJump crI
	order_rj <- order(sim_params$rateJump)
	rateJumps <- sim_params$rateJump[order_rj]
	density_rj <- exp(sim_ll[, stat])
	density_rj <- density_rj[order_rj]
	area <- trapz(rateJumps, density_rj)
	density_rj <- density_rj/area
	
	rmvals <- rep(0, 3)
	rjvals <- rep(0, 3)
	rmvals[1] <- weighted.mean(rateMoves, density_rm)
	rjvals[1] <- weighted.mean(rateJumps, density_rj)

	rm_post_stats <- rbind(rm_post_stats, rmvals)	
	rj_post_stats <- rbind(rj_post_stats, rjvals)	
}


