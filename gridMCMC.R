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

#==================
# set parameters for simulation
#==================
## name the simulation!
nameSimul <- "GRID_33x33_Hop_Jump_BinLik_genSemivar_seed2"
set.seed(2)

## set spam memory options
spam.options(nearestdistnnz=c(13764100,400))

## how many gillespie repetitions per iteration
Nrep <- 600

## size of grid
num.rows <- 33
num.cols <- 33
row.dist <- 10

## distance classes for the general variogram
genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30)) 

## parameters for uniform hop/skip/jump model
limitHopSkip <- 40
limitJump <- 160
rateMove <- 0.04

## the noKernelMultiGilStat normalizes these weights
weightHopInMove <- 1
weightSkipInMove <- 0.0
weightJumpInMove <- 0.1 

# make a map with just x, y
maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)
# set blockIndex to NULL
# no blocks!
blockIndex = NULL

# make stratified matrix (no skips, set blockIndex to NULL)
stratHopSkipJump <- generate_stratified_mat(coords=maps[, c("X", "Y")], limitHopSkip, limitJump, blockIndex=blockIndex)

#===================
# Prep geospatial/coordinate/household data for simulations
#===================
### starting point for simulations
startInfestH <- ceiling(33*(33/2))
startInfestH <- c(startInfestH, startInfestH + 1, startInfestH - 1) 

## plot initially infested houses
dev.new()
par(mfrow = c(3, 1))
infested <- rep(0, length(maps$X))
infested[startInfestH] <- 1
plot_reel(maps$X, maps$Y, infested, base = 0, top = 1)

## create a dummy timeH
timeH <- rep(-2, length(startInfestH))

## let the infestation spread for two years, nbit <- 52 * 2
nbit <- 104

## run 1 gillespie simulation to give second timepoint data 
start <- Sys.time()
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = rateMove, weightHopInMove = weightHopInMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = FALSE, seed = 2)
print(Sys.time() - start)

## obtain stats from the second gillespie simulation
## if(!is.vector(secondTimePointSimul$statsTable)){
## 	statsData <- secondTimePointSimul$statsTable[, 1]
## }else{
## 	statsData <- secondTimePointSimul$statsTable
## }

## plot results of gillespie
binomEndInfested <- secondTimePointSimul$infestedDens
cat("starting # infested:", length(startInfestH), " ending # infested:", length(which(binomEndInfested!=0)), "\n")
plot_reel(maps$X, maps$Y, binomEndInfested, base = 0, top = 1)

## close the device so it prints
dev.off()

# weibull order plotting to check if stats are normal
# for(i in 1:dim(secondTimePointSimul$statsTable)[1])
# {
# 
# order <- order(secondTimePointSimul$statsTable[i, ])
# plot((1:dim(secondTimePointSimul$statsTable)[2])/(dim(secondTimePointSimul$statsTable)[2]+1), secondTimePointSimul$statsTable[i, order], main = i)
# readline()
# 
# }

#==================
# Priors (also the place to change the parameters)
#==================
priorMeans<-c(0.045, 0.05)
priorSdlog <- c(1, 1)
realMeans<-c(rateMove, weightJumpInMove)
sampling <- c("lnorm", "lnorm")
names(priorMeans)<-c("rateMove" , "weightJumpInMove") 
names(sampling)<-c("rateMove", "weightJumpInMove")
names(realMeans)<-names(priorMeans)
names(priorSdlog)<-names(priorMeans)

### Intervals of definition for the parameters
### No longer used, leave in for backward consistency
paramInf<-c(0.002,0.0001)
paramSup<-c(0.30,1000)
names(paramInf)<-names(sampling)
names(paramSup)<-names(sampling)

### LD formalism for Data (no setting should be made past this declaration)
PGF<-function(Data){ # parameters generating functions (for init etc...)
	
	priorMeans<-Data$priorMeans
	priorSdlog<-Data$priorSdlog
	
	values<-rlnorm(length(priorMeans),mean=log(priorMeans), sd=priorSdlog)
	return(values)
}

#=================
# List of data to pass to model + sampler
#=================

MyDataFullSample <- list(y=as.integer(binomEndInfested),
	     trans=NULL,
	     stratHopSkipJump = stratHopSkipJump,
	     blockIndex=blockIndex,
	     dist_out = NULL, # dist_out = makeDistClasses(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals), 
	     infestH=startInfestH,
	     timeH=timeH,
	     endTime=nbit,
	     maps=maps,
	     nbit=nbit,
	     Nrep=Nrep,
	     priorMeans=priorMeans,
	     priorSdlog=priorSdlog,
	     genIntervals=genIntervals,
	     paramInf=paramInf,
	     paramSup=paramSup,
	     PGF=PGF,
	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
	     sampling=sampling # method of sampling parameters
		)

#=================
## Test binomNoKernelModel to make sure something meaningful comes out
#=================
start<-Sys.time()
ModelOutGood<-binomNoKernelModel(priorMeans,MyDataFullSample)
cat(Sys.time()-start, "\n")
start<-Sys.time()
ModelOutBest<-binomNoKernelModel(realMeans,MyDataFullSample)
cat(Sys.time()-start, "\n")

# good should be worse than best (by a fudge factor of 4)
expect_true(ModelOutGood$Dev>ModelOutBest$Dev-4)

#=================
## Make call to MCMC
#=================
MCMC(MyDataFullSample, binomNoKernelModel)

