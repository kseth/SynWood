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
# import params of two maps from Jerusalen
# set parameters for simulation
#==================
## set spam memory options
spam.options(nearestdistnnz=c(13764100,400))

nameSimul<-"JerusalenBinomNoKernelWithWeights" # used by sec_launch.sh to give a name to the output folder
Nrep=100
set.seed(1)
monitor.file<- "thetasamples_all.txt"

## parameters for uniform hop/skip/jump model
limitHopSkip <- 60
limitJump <- 1000
rateMove <- 0.04

## the noKernelMultiGilStat normalizes these weights
weightHopInMove <- 1
weightSkipInMove <- 0.25
weightJumpInMove <- 0.10

## csv of Jerusalen
## "Encuesta_melgar.csv" contains data from the Encuesta for all of Mariana Melgar ~ 2008
## "jerusalen_encuesta_2011_f2_collapsed_KHL.csv" contains the 2011 data for just jerusalen
encuesta_melgar <- read.csv(file = "Encuesta_Melgar.csv")
jerusalen_encuesta <- read.csv(file = "jerusalen_encuesta_2011_full_f2_collapsed_KHL.csv")

keep <- which(encuesta_melgar$unicode %in% jerusalen_encuesta$Unicode)
match <- match(jerusalen_encuesta$Unicode, encuesta_melgar[keep, "unicode"])
maps <- cbind(jerusalen_encuesta[, -which(names(jerusalen_encuesta) %in% c("EASTING", "NORTHING"))], OLDSTATUS = encuesta_melgar[keep, "status"][match])
maps <- cbind(maps, X = jerusalen_encuesta[, "EASTING"], Y = jerusalen_encuesta[, "NORTHING"])
## fields now contained in jerusalen_encuesta:
## "Unicode","POINT_X","POINT_Y", "STATUS","D.x","L.x","V.x","BLOCK_NUM","TOTAL_C","TOTAL_P", "OLDSTATUS", "X", "Y"

startInfestH <- which(maps$OLDSTATUS == 1)
endInfestH <- which(maps$STATUS == 1)
binomEndInfested <- rep(0, length(maps$STATUS))
binomEndInfested[endInfestH] <- 1
blockIndex <- maps$BLOCK_NUM

stratHopSkipJump <- generate_stratified_mat(coords = maps[, c("X", "Y")], limitHopSkip, limitJump, blockIndex)

#===================
# Prep geospatial/coordinate/household data for simulations
#===================

## create a dummy timeH
timeH <- rep(-2, length(startInfestH))

## time between start and end is 3 years, 52*3
nbit <- 156

#==================
# Priors (also the place to change the parameters)
#==================
priorMeans<-c(rateMove, weightJumpInMove, weightSkipInMove)
priorSdlog <- c(1, 1, 1)
realMeans<-c(NA, NA, NA)
sampling <- c("lnorm", "lnorm", "lnorm")
names(priorMeans)<-c("rateMove" , "weightJumpInMove", "weightSkipInMove") 
names(sampling)<-names(priorMeans)
names(realMeans)<-names(priorMeans)
names(priorSdlog)<-names(priorMeans)

### Intervals of definition for the parameters
### No longer used, leave in for backward consistency
paramInf<-c(0.002,0.0001,0.0001)
paramSup<-c(0.30, 1000, 1000)
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

