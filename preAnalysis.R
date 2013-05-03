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

realMeans <- c(0.035, 0.500, 0.80)
priorMeans<-c(0.035, 0.5, 0.80) #should be realMeans (unless want some sort of skewed prior) 
priorSd <- c(0.75, 0.25, 0.20)
priorType <- c("lnorm", "unif", "boundednorm")
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

#========================
# get startInfestH
#========================
# use randominitdays to generate starting point
if(randominitdays == 0){ # no random init, seeding points are starting points
	startInfestH <- startInfestH
}else{ # run gillespie to generate a random starting set from seeding points
	randominitout <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = randominitdays, rateMove = realMeans["rateMove"], weightSkipInMove = 0, weightJumpInMove = realMeans["weightJumpInMove"], Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = NULL, simul=TRUE, getStats = FALSE, seed = seedSimul, dist_out = NULL, map.partitions = NULL, conc.circs = NULL)
	startInfestH <- which(randominitout$infestedDens == 1) 
}


#=========================
# Set fitting method
#=========================
# which statistics to use? 
# choices: "grid", "circles", "semivariance"
useStats <- c("grid", "circles", "semivariance")

# distance classes for the general variogram
# genIntervals <- c(seq(10, 190, 30), seq(240, 340, 50))
genIntervals <- seq(0, 70, 9)
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
# 5000 repetitions from real values
#========================
real_reps <- 5000 

ts <- Sys.time()

real_sims <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = realMeans["rateMove"], weightSkipInMove = 0, weightJumpInMove = realMeans["weightJumpInMove"], Nrep = real_reps, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats=TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles)

te <- Sys.time()

print(te-ts)

real_sims_stats <- real_sims$statsTable

#========================
# Simulations and their statistics
#========================
num_draws <- 2000 #number iterations (draws from prior)
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

# assign which statistics are which (may have to be done manually if change statistics)
grid_stats <- 1:(length(partitionSizes)*6)
grid_var_stats <- grid_stats[which(grid_stats %% 6 == 1)]
grid_count_stats <- grid_stats[which(grid_stats %% 6 == 2)]
grid_regression_stats <- grid_stats[which(grid_stats %% 6 %in% c(3:5, 0))]
circ_stats <- grid_stats[length(grid_stats)] + 1:((length(circleRadii)-1)*2)
semivar_stats <- circ_stats[length(circ_stats)] + 1:((length(genIntervals)-1)*4)
semivar_newnew_stats <- semivar_stats[which(semivar_stats %% 4 %in% 1:2)] 
semivar_oldnew_stats <- semivar_stats[which(semivar_stats %% 4 %in% c(3, 0))] 
num_inf_stats <- semivar_stats[length(semivar_stats)]+1

stop()
#========================
# Calculating likelihoods
#========================
sim_ll <- data.frame()
sim_ll2 <- data.frame()
tll1 <- 0
tll2 <- 0
for(draw in 1:num_draws){

	print(draw)
	
	ll <- rep(0, 10)
	#ll[1] <- synLik(real_sims_stats, sim_stats[draw, ], trans = NULL)
	#ll[2] <- synLik(real_sims_stats[grid_stats, ], sim_stats[draw, grid_stats], trans = NULL)
	ll[3] <- synLik(real_sims_stats[grid_var_stats, ], sim_stats[draw, grid_var_stats], trans = NULL)
	ll[4] <- synLik(real_sims_stats[grid_count_stats, ], sim_stats[draw, grid_count_stats], trans = NULL)
	#ll[5] <- synLik(real_sims_stats[grid_regression_stats, ], sim_stats[draw, grid_regression_stats], trans=NULL)

	start1 <- Sys.time()
	ll[6] <- synLik(real_sims_stats[circ_stats, ], sim_stats[draw, circ_stats], trans = NULL)
	ll[7] <- synLik(real_sims_stats[semivar_stats, ], sim_stats[draw, semivar_stats], trans = NULL)
	ll[8] <- synLik(real_sims_stats[semivar_newnew_stats, ], sim_stats[draw, semivar_newnew_stats], trans = NULL)	
	ll[9] <- synLik(real_sims_stats[semivar_oldnew_stats, ], sim_stats[draw, semivar_oldnew_stats], trans = NULL)	
	end1 <- Sys.time()
	ll[10] <- log(density(real_sims_stats[num_inf_stats, ], from=sim_stats[draw, num_inf_stats], to=sim_stats[draw, num_inf_stats], n=1)$y)
	sim_ll <- rbind(sim_ll, ll)


	ll2 <- rep(0, 10)
	ll2[1] <- synLik.modified(real_sims_stats, sim_stats[draw, ], trans = NULL)
	ll2[2] <- synLik.modified(real_sims_stats[grid_stats, ], sim_stats[draw, grid_stats], trans = NULL)
	ll2[3] <- synLik.modified(real_sims_stats[grid_var_stats, ], sim_stats[draw, grid_var_stats], trans = NULL)
	ll2[4] <- synLik.modified(real_sims_stats[grid_count_stats, ], sim_stats[draw, grid_count_stats], trans = NULL)
	ll2[5] <- synLik.modified(real_sims_stats[grid_regression_stats, ], sim_stats[draw, grid_regression_stats], trans=NULL)
	start2 <- Sys.time()
	ll2[6] <- synLik.modified(real_sims_stats[circ_stats, ], sim_stats[draw, circ_stats], trans = NULL)
	ll2[7] <- synLik.modified(real_sims_stats[semivar_stats, ], sim_stats[draw, semivar_stats], trans = NULL)
	ll2[8] <- synLik.modified(real_sims_stats[semivar_newnew_stats, ], sim_stats[draw, semivar_newnew_stats], trans = NULL)	
	ll2[9] <- synLik.modified(real_sims_stats[semivar_oldnew_stats, ], sim_stats[draw, semivar_oldnew_stats], trans = NULL)
	end2 <- Sys.time()
	ll2[10] <- log(density(real_sims_stats[num_inf_stats, ], from=sim_stats[draw, num_inf_stats], to=sim_stats[draw, num_inf_stats], n=1)$y)

	tll1 <- tll1 + difftime(end1, start1, units = "secs")
	tll2 <- tll2 + difftime(end2, start2, units = "secs")

	sim_ll2 <- rbind(sim_ll2, ll2)

}

names(sim_ll) <- c("all_stats", "grid_stats", "grid_var_stats", "grid_count_stats", "grid_regression_stats", "circ_stats", "semivar_stats", "semivar_new-new_stats", "semivar_old-new_stats", "num_inf_stats") 
names(sim_ll2) <- c("all_stats", "grid_stats", "grid_var_stats", "grid_count_stats", "grid_regression_stats", "circ_stats", "semivar_stats", "semivar_new-new_stats", "semivar_old-new_stats", "num_inf_stats") 

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

## 2D plotting
#uniform xlims, ylims
xlim = c(0.01, 0.06)
ylim = c(0, 1)

dev.new()
par(mfrow = c(2, 5))

for(stat in names(sim_ll)){
	goodLL <- which(sim_ll[, stat] > cutoff_ll[stat])
	orderLL <- order(sim_ll[goodLL, stat])

	plot(sim_params[goodLL, 1][orderLL], sim_params[goodLL, 2][orderLL], col = xtocolors(sim_ll[goodLL, stat][orderLL], crp=colorRampPalette(c("black", "red", "orange", "yellow"))), xlab = "rateMove", ylab = "rateJump", main = stat, xlim = xlim, ylim = ylim, pch = 1)
	abline(h = realMeans["weightJumpInMove"])
	abline(v = realMeans["rateMove"])
}

## 2D contour plotting
dev.new()
par(mfrow = c(2, 5))

for(stat in names(sim_ll)){
	goodLL <- which(sim_ll[, stat] > cutoff_ll[stat])
	orderLL <- order(sim_ll[goodLL, stat])

	plot(sim_params[goodLL, 1][orderLL], sim_params[goodLL, 2][orderLL], col = xtocolors(sim_ll[goodLL, stat][orderLL], crp=colorRampPalette(c("black", "red", "orange", "yellow"))), xlab = "rateMove", ylab = "rateJump", main = stat, xlim = xlim, ylim = ylim, pch = 1)
	abline(h = realMeans["weightJumpInMove"])
	abline(v = realMeans["rateMove"])

	gridfromsample <- grid.from.sample(sim_params[goodLL, 1], sim_params[goodLL, 2], sim_ll[goodLL, stat], steps = 20, tr = 0.01)
	contour(gridfromsample$xs, gridfromsample$ys, gridfromsample$zs, add = TRUE)
}

#======================
# Finding credible intervals
#======================
library(pracma)

## x is the x coordinates
## y is the density at the x coordinates
## all y >= 0
density_crI <- function(x, y, alpha=0.05, sigfig=5){

	#sort (x, y) such that <x> is in order
	orderx <- order(x)
	x <- x[orderx]
	y <- y[orderx]

	alpha <- signif(alpha, sigfig)

	hi <- max(y)
	lo <- min(y)
	mid <- (hi+lo)/2

	ycheck <- y
	ycheck[y < mid] <- 0

	area <- signif(trapz(x, ycheck), sigfig)
	oldarea <- area
	newarea <- 1-alpha
	print(area)

	if(1-area == alpha)
		return(list(crI=quantile(x[y > mid], 0:1), ll=mid))
	if(1-area < alpha){
		mid <- c(mid, (mid+hi)/2)
		lo <- mid[length(mid)-1]
	}else{
		mid <- c(mid, (mid+lo)/2)
		hi <- mid[length(mid)-1]
	}
	while(lo < hi &&  oldarea!=newarea){
	
		oldarea <- area
		ycheck <- y
		ycheck[y < mid[length(mid)]] <- 0
		area <- signif(trapz(x, ycheck), sigfig)
		newarea <- area
		print(area)

		if(1-area == alpha)
			return(list(crI=quantile(x[y > mid[length(mid)]], 0:1), ll=mid[length(mid)]))
		if(1-area < alpha){
			mid <- c(mid, (mid[length(mid)]+hi)/2)
			lo <- mid[length(mid)-1]
		}else{
			mid <- c(mid, (mid[length(mid)]+lo)/2)
			hi <- mid[length(mid)-1]
		}
	}	

	return(list(crI=quantile(x[y > mid[length(mid)-1]], 0:1), ll=mid[length(mid)-1]))
}

rm_post_stats <- data.frame()
rj_post_stats <- data.frame()
smoothed_rm <- list()
smoothed_rj <- list()

for(stat in names(sim_ll)){

	print(stat)
	##ratemove crI
	order_rm <- order(sim_params$rateMove)
	rateMoves <- sim_params$rateMove[order_rm]
	density_rm <- exp(sim_ll[, stat])
	if(max(density_rm) < 1){ #if all sim_lls are too small, the log will need to be transformed
		density_rm <- exp(sim_ll[, stat] - max(sim_ll[, stat]))
	}
	density_rm <- density_rm[order_rm]
	smooth_rm <- smooth.spline(rateMoves, density_rm, spar=0.2)
	smooth_rm$y[smooth_rm$y < 0] <- 0
	area <- trapz(smooth_rm$x, smooth_rm$y)
	smooth_rm$y <- smooth_rm$y/area
	mean_rm <- trapz(smooth_rm$x, smooth_rm$y*smooth_rm$x)
	var_rm <- trapz(smooth_rm$x, smooth_rm$y*smooth_rm$x*smooth_rm$x) - mean_rm^2
	cr_rm <- density_crI(smooth_rm$x, smooth_rm$y)

	##rateJump crI
	order_rj <- order(sim_params$rateJump)
	rateJumps <- sim_params$rateJump[order_rj]
	density_rj <- exp(sim_ll[, stat])
	if(max(density_rj) < 1e-1){ #if all sim_lls are too small, the log will need to be transformed
		density_rj <- exp(sim_ll[, stat] - max(sim_ll[, stat]))
	}
	density_rj <- density_rj[order_rj]
	smooth_rj <- smooth.spline(rateJumps, density_rj, spar=0.2)
	smooth_rj$y[smooth_rj$y < 0] <- 0
	area <- trapz(smooth_rj$x, smooth_rj$y)
	smooth_rj$y <- smooth_rj$y/area
	mean_rj <- trapz(smooth_rj$x, smooth_rj$y*smooth_rj$x)
	var_rj <- trapz(smooth_rj$x, smooth_rj$y*smooth_rj$x*smooth_rj$x) - mean_rm^2
	cr_rj <- density_crI(smooth_rj$x, smooth_rj$y)

	rmvals <- rep(0, 4)
	rjvals <- rep(0, 4)
	rmvals[1] <- mean_rm
	rmvals[2] <- sqrt(var_rm)
	rmvals[3:4] <- cr_rm$crI
	rmvals[5] <- cr_rm$ll
	rjvals[1] <- mean_rj
	rjvals[2] <- sqrt(var_rj)
	rjvals[3:4] <- cr_rj$crI
	rjvals[5] <- cr_rj$ll

	rm_post_stats <- rbind(rm_post_stats, rmvals)	
	rj_post_stats <- rbind(rj_post_stats, rjvals)
	smoothed_rm[[stat]] <- smooth_rm
	smoothed_rj[[stat]] <- smooth_rj	
}

colnames(rm_post_stats) <- c("mean", "sd", "2.5%", "97.5%")
rownames(rm_post_stats) <- names(sim_ll)
colnames(rj_post_stats) <- c("mean", "sd", "2.5%", "97.5%")
rownames(rj_post_stats) <- names(sim_ll)

##plot rate moves
dev.new()
par(mfrow = c(5, 2))
for(stat in names(sim_ll)){
	plot(smoothed_rm[[stat]], main = paste0("ratemove ", stat), type = "l", xlim = c(0, 0.2))
	abline(v=rm_post_stats[stat, 1], col = "blue")
	abline(v=rm_post_stats[stat, 1]+rm_post_stats[stat, 2], col = "blue")
	abline(v=rm_post_stats[stat, 1]-rm_post_stats[stat, 2], col = "blue")
	abline(h=rm_post_stats[stat, 5], col = "red")
	abline(v=realMeans["rateMove"], col = "green")
}

##plot rate jumps
dev.new()
par(mfrow = c(5, 2))
for(stat in names(sim_ll)){
	plot(smoothed_rj[[stat]], main = paste0("ratejump ", stat), type = "l")
	abline(v=rj_post_stats[stat, 1], col = "blue")
	abline(v=rj_post_stats[stat, 1]+rm_post_stats[stat, 2], col = "blue")
	abline(v=rj_post_stats[stat, 1]-rm_post_stats[stat, 2], col = "blue")
	abline(h=rj_post_stats[stat, 5], col = "red")
	abline(v=realMeans["weightJumpInMove"], col = "green")
}

