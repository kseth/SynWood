library(testthat)
library(verification)
source("spatcontrol/spatcontrol.R", chdir=TRUE) #has all the spatial analysis + convergence functions
source("SynWood/logLik.R", chdir = TRUE) # override some sl functions and add synLik
source("SynWood/functions_sampling.R", chdir = TRUE)
source("SynWood/functions_migration.R", chdir = TRUE)
source("SynWood/models.R", chdir = TRUE)
source("SynWood/MCMC.R", chdir = TRUE)

#=====================
# Seeding, Simulation management, and memory/saving options
#=====================
# name the simulation
nameSimul <- "MarianaMelgar_Blocks_SynLik_Semivariance_seed109000_LLJ600"

# the file to store the log of the simulation (i.e. which seed currently on, time of simulation, etc.)
log.file <- "daisyChainLogFile.txt"

# pick the seeds for the simulation
daisyChainSeeds <- 109000 
 
# set spam memory options
spam.options(nearestdistnnz=c(13764100,400))

# is this going to be random data generation
generateData <- FALSE 

#======================
# Set details of grid (simulation environment)
#====================== 
# size of simulated grid
num.rows <- 48 
num.cols <- 48
row.dist <- 10
 
## make a grid map 
if(generateData){
	maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)
}else{ #using mariana melgar
	
	#random stuff for now
	load("byManzVig.img")

	# focus on Mariano Melgar/Paucarpata
	dat<-byManz[byManz$D %in% c(10,13), ]
	maps <- dat[, c("X", "Y")]
	# distances: should be hops up to 600, jumps up to 15000 (all map)
	
}

# parameters for uniform hop/skip/jump model
limitHopSkip <- 600
limitJump <- 15000
lowerLimitJump <- 600 

# set blockIndex to NULL
# no blocks!
blockIndex <- NULL

start <- Sys.time()
# make stratified matrix (no skips, set blockIndex to NULL)
# matrix contains where each house can hop to, where can jump to
stratHopSkipJump <- generate_stratified_mat(coords=maps[, c("X", "Y")], limitHopSkip, limitJump, lowerLimitJump=lowerLimitJump, blockIndex=blockIndex)
print(Sys.time() - start)

#=======================
# Pick spatial locations to seed the starting infestation
#=======================
# pick starting point for simulations
startInfestH <- ceiling(num.rows*(num.rows/2) + num.rows/2)
startInfestH <- c(startInfestH + 1, startInfestH - 1) 

#=======================
# Pick length of simulation
#=======================
# random initial generation, if set to 0, startInfestH will be used as starting points
# otherwise, will generate a random set of starting initial points by gillespie until randominitdays
randominitdays <- 52 

# let the infestation spread for two years, nbit <- 52 * 2
nbit <- 104

#=======================
# Choice of simulation parameters
#=======================
rateMove <- 0.04
weightSkipInMove <- 0.0 #(% of hops we skip)
weightJumpInMove <- 0.10 #(% of times we jump) 

detectRate <- 1 # true detection rate (1 == don't remove data) 
sampleDR <- FALSE # if true, MCMC will sample over error rates
defaultDR <- 1 # DR assumed by multiGilStat (1 if detectRate==1, can set to 0.7, only used if !sampleDR)

rateIntro <- 0.005 # rate of new introductions into the map (1/rateIntro = time to random introduction)
sampleRI <- TRUE
defaultRI <- 0 #RI assumed by multiGilStat (0 if rateIntro==0, can set to any value, only used if !sampleRI)

#=======================
# Priors && Sampling Methodologies
#=======================
priorMeans<-c(0.02, 0.10, 0.80, 0.005) 
priorSd <- c(1, 0.5, 0.20, 1)
priorType <- c("lnorm", "noPrior", "boundednorm", "lnorm")
priorIntervals <- list(NULL, NULL, c(0, 1), NULL) # only considered if priorType is bounded
realMeans<-c(rateMove, weightJumpInMove, detectRate, rateIntro)
sampling <- c("lnorm", "boundednorm", "boundednorm", "lnorm")
sdProposal <- c(0.4, 0.2, 0.2, 0.3)

names(priorMeans)<-c("rateMove" , "weightJumpInMove", "detectRate", "rateIntro") 
names(sampling)<-names(priorMeans)
names(realMeans)<-names(priorMeans)
names(priorSd)<-names(priorMeans)
names(priorType)<-names(priorMeans)
names(priorIntervals) <- names(priorMeans) 

# values with which to initialize the sampler
initValues<-priorMeans
initValues["rateMove"]<-0.02
initValues["weightJumpInMove"]<-0.10
initValues["detectRate"]<-0.80
initValues["rateIntro"]<-0.005

if(!sampleDR){ #if don't want to sample over detectRate
	rem <- which(names(priorMeans) == "detectRate")
	priorMeans<-priorMeans[-rem]
	priorSd<-priorSd[-rem]
	priorType<-priorType[-rem]
	priorIntervals <- priorIntervals[-rem]
	realMeans<-realMeans[-rem]
	sampling<-sampling[-rem]
	sdProposal<-sdProposal[-rem]
	initValues<-initValues[-rem]
}

if(!sampleRI){ #if don't want to sample over rateIntro
	rem <- which(names(priorMeans) == "rateIntro")
	priorMeans<-priorMeans[-rem]
	priorSd<-priorSd[-rem]
	priorType<-priorType[-rem]
	priorIntervals <- priorIntervals[-rem]
	realMeans<-realMeans[-rem]
	sampling<-sampling[-rem]
	sdProposal<-sdProposal[-rem]
	initValues<-initValues[-rem]

}


#=========================
# Set fitting method
# How many repetitions to use
# What kind of grid, circles, distance classes?
#=========================
# how many gillespie repetitions per iteration
Nrep <- 100 

# which likelihood to use? 
useBinLik <- FALSE 
modelToUse<-{if(useBinLik) binomNoKernelModel else noKernelModel}

# which statistics to use? 
# choices: "grid", "circles", "semivariance"
useStats <- c("semivariance") # disregarded if useBinLik

# distance classes for the general variogram
# genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30))
genIntervals <- c(0, cumsum(1:9*75))
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

#=======================
# Running the MCMC for each of the seeds in the chain
#=======================

for(daisyChainNumber in daisyChainSeeds){
	seedSimul <- daisyChainNumber
	monitor.file <- paste0("thetasamples_all", seedSimul, ".txt")
	cat("current seed: ", daisyChainNumber, " monitor.file: ", monitor.file, "\n", file = log.file, append = TRUE)
	ts <- Sys.time()

	if(generateData)
		source("SynWood/dataGenerate.R", chdir = TRUE)
	else{
		startingInfested <- NULL	
		binomialEndInfested <- dat$nbPos>0	
		endInfestH <- which(binomialEndInfested)
		timeH <- -2
		nbit <- 10*52 
		secondTimePointStats <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestH, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro)

		# obtain stats from the second gillespie simulation now messed up via observation error
		if(!is.vector(secondTimePointStats$statsTable)){
			statsData <- secondTimePointStats$statsTable[, 1]
		}else{
			statsData <- secondTimePointStats$statsTable
		}


	}

	#All the data to pass to the model + sampler
	MyDataFullSample <- list(y={if(useBinLik) binomialEndInfested else statsData},
	     trans=NULL,
	     stratHopSkipJump = stratHopSkipJump,
	     blockIndex=blockIndex,
	     dist_out = {if(!useBinLik && ("semivariance" %in% useStats)) bin_dist_out else NULL},
	     map.partitions = {if(!useBinLik && ("grid" %in% useStats)) map.partitions else NULL}, 
	     conc.circs = {if(!useBinLik && ("circles" %in% useStats)) circles else NULL}, 
	     useStats = useStats,
	     infestH=startingInfested,
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
	     defaultRI=defaultRI,
	     genIntervals=genIntervals,
	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
	     sampling=sampling
		)
	
	#Test modelToUse to make sure something meaningful comes out
	start<-Sys.time()
	ModelOutGood<-modelToUse(priorMeans,MyDataFullSample)
	cat(Sys.time()-start, "\n")
	start<-Sys.time()
	ModelOutBest<-modelToUse(realMeans,MyDataFullSample)
	cat(Sys.time()-start, "\n")
	#good should be worse than best (ideally, need -4 because simulations may not be ideal)
	#expect_true(ModelOutGood$Dev>ModelOutBest$Dev-4)

	#run the MCMC function
	MCMC(MyDataFullSample, Model=modelToUse, sdprop=sdProposal, monitor.file=monitor.file)

	te <- Sys.time()
	cat("time: ", as.numeric(difftime(te, ts, unit = "mins")), " mins\n\n", file = log.file, append = TRUE)

	# reset all global variables
	# reset minLLever
	source("SynWood/models.R", chdir = TRUE)
}

