library(testthat)
library(verification)
source("RanalysisFunctions.R")
source("logLik.r") # override some sl functions and add synLik
source("functions_sampling.R")
source("convergence_test.r")
source("extrapol_field.R")
source("functions_migration.R")
source("models.R")
source("MCMC.R")

## name the simulation
nameSimul <- "GRID_50X50_HopJump_BinLik_10:109*1000_J30_100"

## the file to store the log of the simulation (i.e. which seed currently on, time of simulation, etc.)
log.file <- "daisyChainLogFile.txt"

## pick the seeds for the simulation
daisyChainSeeds <- 10:109*1000
 
## set spam memory options
spam.options(nearestdistnnz=c(13764100,400))

## how many gillespie repetitions per iteration
Nrep <- 400 

## Make simulation messy or not messy
detectRate <- 1 # true detection rate 
sampleDR <- FALSE # if true, MCMC will sample over error rates
defaultDR <- 1 # DR assumed by multiGilStat (should be 1 if detectRate==1, can set to 0.7, only used if sampleDR is FALSE)
 
## size of grid
num.rows <- 50 
num.cols <- 50
row.dist <- 10
 
## parameters for uniform hop/skip/jump model
limitHopSkip <- 30
limitJump <- 100
lowerLimitJump <- 30 
rateMove <- 0.04

## the noKernelMultiGilStat normalizes these weights
weightSkipInMove <- 0.0
weightJumpInMove <- 0.10 

## which likelihood to use? 
useBinLik <- TRUE

## which statistics to use? 
## choices: "grid", "circles", "semivariance"
useStats <- c("circles") # disregarded if useBinLik == TRUE

## make a map with just x, y
maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)

## distance classes for the general variogram
genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30))
## genIntervals <- seq(10, 40, 15)  # if want to combine with grid

## bin the map into different distances classes
bin_dist_out <- makeDistClasses(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals)

## partition the map
map.partitions <- list()
length(map.partitions) <- 6 #6 different grid partitions will be used

map.partitions[[1]] <- partitionMap(maps$X, maps$Y, 12) #into 12 by 12 (each cell 3 by 3)
map.partitions[[2]] <- partitionMap(maps$X, maps$Y, 9)  #into 9 by 9 (each cell 4 by 4)
map.partitions[[3]] <- partitionMap(maps$X, maps$Y, 6)  #into 6 by 6 (each cell 6 by 6)
map.partitions[[4]] <- partitionMap(maps$X, maps$Y, 4)  #into 4 by 4 (each cell 9 by 9)
map.partitions[[5]] <- partitionMap(maps$X, maps$Y, 3)  #into 3 by 3 (each cell 12 by 12)
map.partitions[[6]] <- partitionMap(maps$X, maps$Y, 2)  #into 2 by 2 (each cell 18 by 18)

## set blockIndex to NULL
## no blocks!
blockIndex = NULL

## make stratified matrix (no skips, set blockIndex to NULL)
stratHopSkipJump <- generate_stratified_mat(coords=maps[, c("X", "Y")], limitHopSkip, limitJump, lowerLimitJump=lowerLimitJump, blockIndex=blockIndex)

### pick starting point for simulations
startInfestH <- ceiling(num.rows*(num.rows/2) + num.rows/2)
startInfestH <- c(startInfestH, startInfestH + 1, startInfestH - 1) 

## let the infestation spread for two years, nbit <- 52 * 2
nbit <- 104

## make the concentric circles
circleRadii <- c(0, 20, 35, 50, 80, 110, 155, 200)
circles <- conc.circles(maps$X, maps$Y, circleRadii, startInfestH) 

## create a dummy timeH
timeH <- rep(-2, length(startInfestH))
for(daisyChainNumber in daisyChainSeeds){
	seedSimul <- daisyChainNumber
	monitor.file <- paste0("thetasamples_all", seedSimul, ".txt")
	cat("current seed: ", daisyChainNumber, " monitor.file: ", monitor.file, "\n", file = log.file, append = TRUE)
	ts <- Sys.time()
	source("gridMCMC.R")
	te <- Sys.time()
	cat("time: ", as.numeric(difftime(te, ts, unit = "mins")), " mins\n\n", file = log.file, append = TRUE)

	## reset all global variables
	## reset minLLever
	source("models.R")
}
