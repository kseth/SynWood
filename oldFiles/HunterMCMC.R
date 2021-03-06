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
## set spam memory options
spam.options(nearestdistnnz=c(13764100,400))

## distance classes for the general variogram
genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30)) 

nameSimul<-"HunterBinomMessyData_notopen_notdetected" # used by sec_launch.sh to give a name to the output folder
Nrep=500
set.seed(1)
monitor.file<- "thetasamples_all.txt"

## parameters for uniform hop/skip/jump model
limitHopSkip <- 60
limitJump <- 1000
rateMove <- 0.04

## the noKernelMultiGilStat normalizes these weights
weightHopInMove <- 1
weightSkipInMove <- 0.35
weightJumpInMove <- 0.10 

## csv of hunter map
maps.tot<-read.csv("maps_hunter_blocks.csv")

# focusing on hunter for now
maps<-maps.tot[which(maps.tot$D==7 & maps.tot$X>226000 & maps.tot$Y>8179000 & maps.tot$Y<8180000),]
blockIndex <- maps$block_num

# if any NAs in block index (corrects by making NA = new, different blocks)
if(any(is.na(blockIndex))){
	bad <- which(is.na(blockIndex))
	blockIndex[bad] <- max(blockIndex[-bad])+(1:length(bad))
}

stratHopSkipJump <- generate_stratified_mat(coords = maps[, c("X", "Y")], limitHopSkip, limitJump, blockIndex)

#===================
# Prep geospatial/coordinate/household data for simulations
#===================
### starting point for simulations
startInfestH <- c(3200, 1, 10, 3210, 8, 15, 14, 3199, 3220, 16, 20, 25)

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
secondTimePointSimul <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = startInfestH, timeH=timeH, endTime = nbit, rateMove = rateMove, weightHopInMove = weightHopInMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE, getStats = FALSE)
print(Sys.time() - start)

## plot results of gillespie
binomEndInfested <- as.integer(secondTimePointSimul$infestedDens)
cat("starting # infested:", length(startInfestH), " ending # infested:", length(which(binomEndInfested!=0)), "\n")
plot_reel(maps$X, maps$Y, binomEndInfested, base = 0, top = 1)

binomEndInfested <- simulObserved(binomEndInfested, 0.7, 0.9)
print(length(which(binomEndInfested == 1)))
plot_reel(maps$X, maps$Y, binomEndInfested, base = 0, top = 1)

#==================
# Priors (also the place to change the parameters)
#==================
priorMeans<-c(0.045, 0.05, 0.25)
priorSdlog <- c(1, 1, 1)
realMeans<-c(rateMove, weightJumpInMove, weightSkipInMove)
sampling <- c("lnorm", "lnorm", "lnorm")
names(priorMeans)<-c("rateMove" , "weightJumpInMove", "weightSkipInMove") 
names(sampling)<-names(priorMeans)
names(realMeans)<-names(priorMeans)
names(priorSdlog)<-names(priorMeans)

### Intervals of definition for the parameters
### No longer used, leave in for backward consistency
paramInf<-c(0.002,0.0001,0.0001)
paramSup<-c(0.30,1000, 1000)
names(paramInf)<-names(priorMeans)
names(paramSup)<-names(priorMeans)

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

MyDataFullSample <- list(y=binomEndInfested,
	     trans=NULL,
	     stratHopSkipJump = stratHopSkipJump,
	     blockIndex=blockIndex,
	     dist_out = NULL, 
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
# ModelOutBest<-Model(realMeans,MyDataFullSample)
# expect_true(ModelOutGood$Dev>ModelOutBest$Dev-4)

# theta<-realMeans
# theta["rateMove"]<-log(0.5)
# ModelOutBest<-Model(theta,MyDataFullSample)

# weibull order plotting to check if stats are normal
# for(i in 1:dim(outBase$statsTable)[1])
# {
#  
# order <- order(outBase$statsTable[i, ])
# plot((1:dim(outBase$statsTable)[2])/(dim(outBase$statsTable)[2]+1), outBase$statsTable[i, order], main = i)
# Sys.sleep(0.5)
# 
# }

#=================
## Make call to MCMC
#=================
MCMC(MyDataFullSample, binomNoKernelModel)

