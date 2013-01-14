# functions for the migration/spatial stats
# data for the tests
library(geoR)
library(testthat)
source("RanalysisFunctions.R")
source("maps_basic_regression.R") # gives maps
genIntervals <- seq(0, 250, 15) # distance classes for the general variogram

#======================
# Generate hop/skip/jump matrixes
#======================

##
# returns matrixes hopMat, skipMat, jumpMat
# with 1 if can hop/skip/jump there
# only need to pass the @colindices and @rowpointers to the c methods
##

##
# tests, add to test-functions_migration when have time
# length(hopMat@entries) + length(skipMat@entries) == length(dist_mat_hop_skip@entries)
# diag(hopMat) == diag(skipMat) == diag(jumpMat) == 0
##
generate_stratified_mat <- function(coords, limitHopSkip, limitJump, blockIndex)
{

	# make same block matrix of households in the same block
	SB <- nearest.dist(x=cbind(blockIndex,rep(0,length(blockIndex))), method="euclidian", upper=NULL, delta=0.1)
	SB@entries <- rep(1,length(SB@entries))
	SB<-as.spam(SB)

	# make distance matrixes for two thresholds: limitHopSkip, limitJump
	dist_mat_hop_skip <- nearest.dist(coords, y=NULL, method = "euclidian", delta = limitHopSkip, upper = NULL)
	dist_mat_jump <- nearest.dist(coords, y=NULL, method = "euclidian", delta = limitJump, upper = NULL)
	
	# remove diagonal values by cleaning up distances of zero
	dist_mat_hop_skip <- cleanup(dist_mat_hop_skip)
	dist_mat_jump <- cleanup(dist_mat_jump)
	
	#define hopMat, skipMat
	#hopMat if in same block
	#skipMat if across streets
	hopMat <- dist_mat_hop_skip
	skipMat <- dist_mat_hop_skip
	
	hopMat@entries <- rep(1, length(hopMat@entries))
	skipMat@entries <- rep(1, length(skipMat@entries))
	
	hopMat <- hopMat * SB
	hopMat <- as.spam(hopMat)

	skipMat <- skipMat - hopMat
	skipMat <- as.spam(skipMat)

	#define jumpMat
	jumpMat <- dist_mat_jump
	jumpMat@entries <- rep(1, length(jumpMat@entries))
	jumpMat <- as.spam(jumpMat)

	return(list(hopMat = hopMat, skipMat = skipMat, jumpMat = jumpMat))

}

#======================
# Generating probability matrix
#======================

generate_prob_mat <- function(halfDistJ, 
			      halfDistH, 
			      # rateMove, 
			      useDelta, 
			      delta, 
			      rateHopInMove, 
			      rateSkipInMove, 
			      rateJumpInMove, 
			      threshold, 
			      sp, 
			      dist_mat,
			      SB,
			      cumul=FALSE # out cumulative proba for C
){

	### Generating hop, skip, jump matrices

	# decreasing spatial link hop/skip
	weightMat<-exp(-dist_mat/halfDistH)
	diag(weightMat)<- rep(0,dim(weightMat)[1]) # emigrants;  can't go from house to same house

	# hop
	hopMat<-SB*weightMat 

	# skip
	DB <- abs(SB-1)
	skipMat<-DB*weightMat

	# jump 
	jumpMat<-exp(-dist_mat/halfDistJ)
	diag(jumpMat)<- rep(0,dim(jumpMat)[1]) # emigrants; can't go from house to same house

	# normalize jumps
	rsumQ<- jumpMat %*% rep(1,num_obs)
	probMatJ <- as.matrix(jumpMat * as.vector(1/rsumQ))

	if(useDelta){
		hopMat<-delta*skipMat+hopMat
	}else{
		# normalize skips
		rsumQ <- skipMat %*% rep(1,num_obs)
		probMatS <- as.matrix(skipMat * as.vector(1/rsumQ))
	}

	# normalize hops
	rsumQ<- hopMat %*% rep(1,num_obs)
	probMatH <- as.matrix(hopMat * as.vector(1/rsumQ))

	#======================
	# put together
	#======================
	if(useDelta){
		probMat <- (1-rateJumpInMove)*probMatH + rateJumpInMove*probMatJ
	}else{
		probMat <- rateHopInMove*probMatH + rateSkipInMove*probMatS + rateJumpInMove*probMatJ
	}

	diag(probMat) <- rep(0, dim(probMat)[1]) # emigrants; can't go from house to same house

	# #multiply the probMat by the probability that the infestations spreads
	# #rateMove from param.r
	# probMat <- probMat*rateMove

	if(cumul){
		#redefine probMat as the cumulative sum on each line
		for(row in 1:L)
		{
			probMat[row, ]<-cumsum(probMat[row, ])
		}

		# # transpose so that redeable by lines in C
		probMat<-t(probMat)
	}

	return(probMat)
}

# fast function keeping varibles for skip/hop/jumps but not as fix rates
fast_prob_mat <- function(halfDistJ, 
			      halfDistH, 
			      # rateMove, 
			      useDelta, 
			      delta, 
			      rateHopInMove, 
			      rateSkipInMove, 
			      rateJumpInMove, 
			      threshold, 
			      sp, 
			      dist_mat,
			      SB,
			      cumul=FALSE # out cumulative proba for C
){


	### Generating hop, skip, jump matrices

	# decreasing spatial link hop/skip
	weightMat<-exp(-dist_mat/halfDistH)
	diag(weightMat)<- rep(0,dim(weightMat)[1]) # emigrants;  can't go from house to same house

	# hop
	hopMat<-SB*weightMat 

	# skip
	DB <- abs(SB-1)
	skipMat<-DB*weightMat

	hopMat <- hopMat + skipMat*delta

	# jump 
	jumpMat<-exp(-dist_mat/halfDistJ)
	diag(jumpMat)<- rep(0,dim(jumpMat)[1]) # emigrants; can't go from house to same house

	probMat<-jumpMat+hopMat

	# # normalize hops
	# rsumQ<- hopMat %*% rep(1,num_obs)
	# probMatH <- as.matrix(hopMat * as.vector(1/rsumQ))

	# # normalize skips
	# rsumQ <- skipMat %*% rep(1,num_obs)
	# probMatS <- as.matrix(skipMat * as.vector(1/rsumQ))

	# # normalize jumps
	# rsumQ<- jumpMat %*% rep(1,num_obs)
	# probMatJ <- as.matrix(jumpMat * as.vector(1/rsumQ))

	# #======================
	# # put together
	# #======================
	# probMat <- rateHopInMove*probMatH + rateSkipInMove*probMatS + rateJumpInMove*probMatJ
	# diag(probMat) <- rep(0, dim(probMat)[1]) # emigrants; can't go from house to same house

	# #multiply the probMat by the probability that the infestations spreads
	# #rateMove from param.r
	# probMat <- probMat*rateMove

	if(cumul){
		#redefine probMat as the cumulative sum on each line
		for(col in 1:L){
			probMat[,col]<-cumsum(probMat[,col])
		}
	}

	# # transpose so that redeable by lines in C
	# probMat<-t(cumulProbMat)
	# }

	return(probMat)
}

#=============================
# Migration functions
#=============================
##  define the migration functions

basicMigration<-function(infestH, probMat){
  newInfestHouse <- infestH

  for(houseInit in infestH)
  {

    prob <- probMat[houseInit, ]
    rand <- runif(dim(probMat)[1])
    newInfestHouse <- union(newInfestHouse, which(prob>rand))
  }

  return(newInfestHouse)

}

# semi-gillespie: draw the next time of event for each house
#                 not nearly as fast as the gillespie but may 
#                 be easier to use with differing probabilities of events
#                 for different nodes
loopMethod <- function(events_time, events_loc, probMat, infestH, timeH, currentTime, endTime)
{
	if(currentTime >= endTime)
		return(list(infestH, timeH))
	
	#remove the house/time responsible for the current event from the queue
	events_time <- events_time[-1]
	current_loc <- events_loc[1]
	events_loc <- events_loc[-1]

	#need to pick the newly infested house from current_loc
	newInfested <- sample(x = dim(probMat)[1], size = 1, prob = probMat[current_loc, ])
	
		
	#put the new infested house and the current infested house back into the pile
	timeToNext <- currentTime + rexp(2, rate = 1/scale)
	#timeToNext <- currentTime + rweibull(2, shape = 1, scale = scale)

	if(!(newInfested %in% infestH))
	{
		infestH <- c(infestH, newInfested)
		timeH <- c(timeH, currentTime)

		events_time <- c(events_time, timeToNext)
		events_loc <- c(events_loc, current_loc, newInfested)
	}else
	{
		events_time <- c(events_time, timeToNext[1])
		events_loc <- c(events_loc, current_loc)
	}

	sortTime <- order(events_time)
	events_time <- events_time[sortTime]
	events_loc <- events_loc[sortTime]
	
	#cat("current time", currentTime, "infested", length(infestH), "\n")
	#cat("time: ", events_time, "\n")
	#cat("locations: ", events_loc, "\n\n")

	return(loopMethod(events_time, events_loc, probMat, infestH, timeH, events_time[1], endTime))
}

# # gillespie with multiple snapshots possible
# useless, just use ages from normal gillespie as demonstrated in basic_regression.R and infestSerieToMaps() in this file
# gillespie.mo<-function(probMat, # matrix with probability to end up in given house given "departure" from an initial house
# 		      cumulProbMat, # same but cumulative probability for one house of departure
# 		      infestH,  # vector of infestation
# 		      timeH,     # vector of time of infestation
# 		      snapTimes,# time of snapshots, the last one being the end time of simulation
# 		      scale,     # 1/(rate of departure per house)
# 		      toNextEvent=rweibull(1, shape = 1, scale = scale/length(infestH)), # time to the next event
# 		      currentTime=0 # initial time
# 		      ){
# 	snapshots<-list()
# 	for(itime in 1:length(snapTimes)){
# 		snapshots[[itime]]<-gillespie(probMat,cumulProbMat,infestH,timeH,snapTimes[itime],scale,toNextEvent,currentTime=currentTime)
# 		snapshots[[itime]]$time<-snapTimes[itime]
# 		currentTime <- snapTimes[itime]
# 		infestH<-snapshots[[itime]]$infestOrder
# 		timeH<-snapshots[[itime]]$infestTime
# 		toNextEvent<-snapshots[[itime]]$toNextEvent
# 		currentTime<-snapshots[[itime]]$time
# 	}
# 	return(snapshots)
# }
# test
# set.seed(777)
# out <- gillespie.mo(probMat, cumulProbMat, infestH, timeH, snapTimes=snapTimes, scale)
# 
# par(mfcol=c(2,length(snapTimes)))
# for(itime in 1:length(snapTimes)){
# 	outT<-out[[itime]]
# 	maps<-infestSerieToMaps(outT,sp)
# 	plot_reel(maps$X,maps$Y,maps$infest,base=0,main=paste(outT$time,"time unit (weeks)"))
# 	plot_reel(maps$X,maps$Y,log(maps$ages+1),base=0,top=max(log(maps$ages+1)))
# }



# with one start, one stop
gillespie <- function(probMat, # matrix with probability to end up in given house given "departure" from an initial house
		      cumulProbMat, # same but cumulative probability for one house of departure
		      infestH,  # vector of infestation
		      timeH,     # vector of time of infestation
		      endTime,  # end time of simulation
		      scale,     # 1/(rate of departure per house)
		      toNextEvent=rweibull(1, shape = 1, scale = scale/length(infestH)), # time to the next event
		      currentTime=0 # initial time
		      ){
	#cumulProbMat is not used in this function (is used in .C function)
	#pass as parameter for homogeneity	

	#toNextEvent - the time to the toNextEvent
	# cat("begin:", currentTime,"first event:",currentTime+toNextEvent,"end:",endTime,"\n")
	
	while(currentTime + toNextEvent <= endTime)
	{
		currentTime <- currentTime + toNextEvent

		#pick a location to be the infesting house
		current_loc <- as.integer(runif(1, min = 1, max = length(infestH)+1))
		current_loc <- infestH[current_loc]
		
		#pick a new house to become infested from the infesting house
		newInfested <- sample(x = dim(probMat)[1], size = 1, prob = probMat[current_loc, ])
		
		if(!(newInfested %in% infestH))
		{
			infestH <- c(infestH, newInfested)
			timeH <- c(timeH, currentTime)
		}
	
		#calculate the time to the next event again
		toNextEvent <- rweibull(1, shape = 1, scale = scale/length(infestH))
		
	}
	return(list(infestOrder=infestH,infestTime=timeH,toNextEvent=toNextEvent,time=endTime))
}

# import C functions if possible
importOk<-try(dyn.load("functions_migration.so"), silent=TRUE)

# define migration functions and 
# change all entries in A bigger than maxtobesetnull to 0
if(class(importOk)!="try-error"){

	gillespie<- function(probMat, cumulProbMat, infestH, timeH, endTime, scale,seed=runif(1, 1, 2^31-1)){
		
		#probMat is not used here, is used in R function
		#pass for homogeneity

		#for proper seeding of stochastic simulation	

		L<-dim(cumulProbMat)[1]
		indexInfest <- rep(-1, L)
		timeI <- rep(-1, L)
		infested <- rep(0, L)
		indexInfest[1:length(infestH)] <- infestH - 1
		timeI[1:length(timeH)] <- timeH
		infested[infestH] <- 1
		endIndex<-length(infestH)-1
		out<- .C("gillespie",
			 infested = as.integer(infested),
			 endIndex = as.integer(endIndex),
			 L = as.integer(L),
			 probMat = as.numeric(cumulProbMat),
			 endTime = as.numeric(endTime),
			 indexInfest = as.integer(indexInfest),
			 timeI = as.numeric(timeI),
			 scale = as.numeric(1/scale),
			 seed = as.integer(seed))

		infestH<-out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1 
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]

		# return(list(infestH, timeH))
		return(list(infestOrder=infestH,infestTime=timeH,time=out$endTime))
	}

	#### statistics functions
	####    with their tests
	# get the indexes of distance class for each pair

	makeDistClasses<-function(X,Y,breaks){
		cbin<-rep(0,length(breaks)-1)
		CClassIndex<-dists<-rep(0,length(X)^2)
		out<-.C("makeDistClasses",
			xc=as.numeric(X),
			L=as.integer(length(X)),
			yc=as.numeric(Y),
			cbin=as.integer(cbin),
			CClassIndex=as.integer(CClassIndex),
			dists=as.numeric(dists),
			nbbreaks=as.integer(length(breaks)),
			breaks=as.numeric(breaks),
			maxdist = as.numeric(max(breaks))
			)
		return(list(dists=out$dists,
			    CClassIndex=out$CClassIndex,
			    classSize=out$cbin))
	}

	# make distance classes taking into account streets	
	makeDistClassesWithStreets<-function(X,Y, breaks, blockIndex){
	
	  	# checks to avoid segfault
	  if(length(X)!=length(Y) || length(X) != length(blockIndex))
	    stop(paste("makeDistClassesWithStreets Abort\n
		       length(X)=",length(X),
		       "length(Y):",length(Y),
		       "length(blockIndex):",length(blockIndex)))

	  	# declare outputs
		cbin<-rep(0,length(breaks)-1)
		cbinas<-rep(0, length(breaks)-1)
		cbinsb<-rep(0, length(breaks)-1)
		CClassIndex<-dists<-rep(0,length(X)^2)
		
		out<-.C("makeDistClassesWithStreets",
			xc=as.numeric(X),
			L=as.integer(length(X)),
			yc=as.numeric(Y),
			cbin=as.integer(cbin),
			cbinas=as.integer(cbinas),
			cbinsb=as.integer(cbinsb),
			CClassIndex=as.integer(CClassIndex),
			dists=as.numeric(dists),
			nbbreaks=as.integer(length(breaks)),
			breaks=as.numeric(breaks),
			maxdist = as.numeric(max(breaks)),
			blockIndex = as.integer(blockIndex)
			)

		return(list(dists=out$dists,
			    CClassIndex=out$CClassIndex,
			    classSize=out$cbin,
			    classSizeSB=out$cbinsb,
			    classSizeAS=out$cbinas))
	}

	
	
	# NB: dist_indices are in C convention beginning at 0
	# see makeDistClasses()
	variogFromIndices<-function(CClassIndex,vectData,classSize){
		stats<-rep(0,2*length(classSize))
		out<-.C("modBinIt"
			,n=as.integer(length(vectData))
			,dist_index=as.integer(CClassIndex)
			,inf_data=as.numeric(vectData)
			,cbin=as.integer(classSize)
			,stats=as.double(stats)
			,nbins=as.integer(length(classSize)+1) # nb of breaks
			)
		# cat("stats",out$stats,"\n")
		variog<-out$stats[1:length(classSize)]
		sdvariog<-out$stats[(length(classSize)+1):length(stats)]
		return(list(variog=variog,sdvariog=sdvariog))
	}
		
	
	# pass cumulProbMat and dist_out to make this method faster	
	multiGilStat<-function(cumulProbMat, blockIndex, infestH, timeH, endTime, rateMove, Nrep, coords, breaksGenVar, seed=1, simul=TRUE, getStats=TRUE, halfDistJ = -1, halfDistH = -1, useDelta = -1, delta = -1, rateHopInMove = -1, rateSkipInMove = -1, rateJumpInMove = -1, breaksStreetVar = breaksGenVar, dist_out = NULL){
		
		# seed <- runif(1, 1, 2^31-1)
		#for random seeding of stochastic simulation	

	  	L<-dim(coords)[1]
		indexInfest <- rep(-1, L)
		timeI <- rep(-1, L)
		infested <- rep(0, L)
		indexInfest[1:length(infestH)] <- infestH - 1
		timeI[1:length(timeH)] <- timeH
		infested[infestH] <- 1
		infestedDens<-rep(0,length(infested))
		if(getStats){	
			# stats selection
			# need to implement system where we can add and remove stats
			###===================================
			## CURRENT STATS:
			## General Semivariance
			## General Semivariance Std. Dev.
			## Interblock Semivariance
			## Interblock Semivariance Std. Dev.
			## Intrablock Semivariance
			## Intrablock Semivariance Std. Dev.
			## By block Semivariance
			## By block Semivariance Std. Dev.
			##	= 8 * length(cbin)
			## Number Infested Houses
			## Number Infested Blocks
			## (Infested Houses)/(Infested Blocks)
			##	= 8 * length(cbin) + 3
			###===================================
 	      		sizeVvar<-8*length(cbin)
			nbStats<- sizeVvar + 3
			statsTable<-mat.or.vec(nbStats,Nrep)
		
			if(is.null(dist_out))	
				dist_out <- makeDistClassesWithStreets(as.vector(coords[, 1]), as.vector(coords[, 2]), breaksGenVar, blockIndex)
			
			dist_mat <- dist_out$dists
			dist_indices <- dist_out$CClassIndex
			cbin <- dist_out$classSize
			cbinas <- dist_out$classSizeAS
			cbinsb <- dist_out$classSizeSB
	
		}else{
 	      		sizeVvar<-0
			nbStats<-0
			statsTable<-0
			dist_mat <- 0
			dist_indices <- 0
			cbin <- 0
			cbinas <- 0
			cbinsb <- 0
		}

		# if cumulProbMat is not passed, create blank cumulProbMat for C computation
		# set useProbMat to FALSE	
		if(is.null(cumulProbMat)){
			cumulProbMat <- mat.or.vec(L, L)
			useProbMat <- FALSE
		}else{ # else pass dummy dist_mat so as not to take up memory space	
			dist_mat <- 0
            		useProbMat <- TRUE
		}

		out<- .C("multiGilStat",
			 # simulation parameters
			 probMat = as.numeric(cumulProbMat),
			 useProbMat = as.integer(useProbMat),
			 distMat = as.numeric(dist_mat),
		         halfDistJ = as.numeric(halfDistJ),
			 halfDistH = as.numeric(halfDistH), 
			 useDelta = as.integer(useDelta), 
			 delta = as.numeric(delta), 
			 rateHopInMove = as.numeric(rateHopInMove), 
			 rateSkipInMove = as.numeric(rateSkipInMove), 
			 rateJumpInMove = as.numeric(rateJumpInMove), 
			 blockIndex = as.integer(blockIndex),
			 simul = as.integer(simul),
			 infested = as.integer(infested),
			 infestedDens = as.numeric(infestedDens),
			 endIndex = as.integer(length(infestH) - 1),
			 L = as.integer(L),
			 endTime = as.numeric(endTime),
			 indexInfest = as.integer(indexInfest),
			 timeI = as.numeric(timeI),
			 rateMove = as.numeric(rateMove),
			 seed = as.integer(seed),
			 Nrep = as.integer(Nrep),
			 # stats
			 getStats=as.integer(getStats),
			 nbins = as.integer(length(breaksGenVar)),
			 cbin = as.integer(cbin),
			 cbinas = as.integer(cbinas),
			 cbinsb = as.integer(cbinsb),
			 indices = as.integer(dist_indices),
			 statsTable = as.numeric(statsTable),
			 nbStats = as.integer(nbStats),
             		 sizeVvar = as.integer(sizeVvar)
			 )

		out$infestedDens<-out$infestedDens/Nrep;
	
		# make matrix out of statsTable
		out$statsTable<-matrix(out$statsTable,byrow=FALSE,ncol=Nrep)
		
		# remove interblock and intrablock stats
		# keep:
		# general semivar + stdev
		# sameblock - acrossstreets semivar + stdev
		# inf.house, inf.block, and inf.house/inf.block count
		keepable<-c(1:(2*length(cbin)), 6*length(cbin)+1:(2*length(cbin)), sizeVvar+(1:3))
		# clean away NANs introduced in C
		notNAN <- which(!is.nan(out$statsTable[, 1]))
		
		keep<-intersect(notNAN,keepable)
		out$statsTable<-out$statsTable[keep, ]

		infestH <- out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1
		out$indexInfest <- infestH
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]
		out$timeI <- timeH

		return(out)
	}

	## stratHopSkipJump is the result of a call to generate_stratified_mat
	## list with spam matrices of hoppable, skippable, jumpable locations for each house
	## this method does not use the kernels or delta
	## instead need to pass weightHopInMove (set to 1), weightSkipInMove, weightJumpInMove
	##	these parameters are dealt with by the gillespie
	noKernelMultiGilStat <- function(stratHopSkipJump, blockIndex, infestH, timeH, endTime, rateMove, weightHopInMove, weightSkipInMove, weightJumpInMove, Nrep, coords, breaksGenVar, seed=1, simul=TRUE, getStats=TRUE, breaksStreetVar = breaksGenVar, dist_out = NULL){

		# convert weightHop, weightSkip, weightJump to rates by normalizing (rates needed by c code)
		rateHopInMove <- weightHopInMove/(weightHopInMove+weightSkipInMove+weightJumpInMove)	
		rateSkipInMove <- weightSkipInMove/(weightHopInMove+weightSkipInMove+weightJumpInMove)	
		rateJumpInMove <- weightJumpInMove/(weightHopInMove+weightSkipInMove+weightJumpInMove)	
	
		# seed <- runif(1, 1, 2^31-1)
		# for random seeding of stochastic simulation	

	  	L<-dim(coords)[1]
		indexInfest <- rep(-1, L)
		timeI <- rep(-1, L)
		infested <- rep(0, L)
		indexInfest[1:length(infestH)] <- infestH - 1
		timeI[1:length(timeH)] <- timeH
		infested[infestH] <- 1
		infestedDens<-rep(0,length(infested))


		# if stratHopSkipJump, throw error 	
		if(is.null(stratHopSkipJump)){
			stop("need to pass a stratified hop/skip/jump matrix; see generate_stratifed_mat")
		}

		if(getStats){

			if(is.null(dist_out))	
				dist_out <- makeDistClassesWithStreets(as.vector(coords[, 1]), as.vector(coords[, 2]), breaksGenVar, blockIndex)
			
			dist_indices <- dist_out$CClassIndex
			cbin <- dist_out$classSize
			cbinas <- dist_out$classSizeAS
			cbinsb <- dist_out$classSizeSB

			# stats selection
			# need to implement system where we can add and remove stats
			###===================================
			## CURRENT STATS:
			## General Semivariance
			## General Semivariance Std. Dev.
			## Interblock Semivariance
			## Interblock Semivariance Std. Dev.
			## Intrablock Semivariance
			## Intrablock Semivariance Std. Dev.
			## By block Semivariance
			## By block Semivariance Std. Dev.
			##	= 8 * length(cbin)
			## Number Infested Houses
			## Number Infested Blocks
			## (Infested Houses)/(Infested Blocks)
			##	= 8 * length(cbin) + 3
			###===================================
       			sizeVvar<-8*length(cbin)
			nbStats<- sizeVvar + 3
			statsTable<-mat.or.vec(nbStats,Nrep)
		}else{
			dist_indices <- 0
			cbin <- 0
			cbinas <- 0
			cbinsb <- 0
			sizeVvar <- 0
			nbStats <- 0
			statsTable <- 0
		}

		out<- .C("noKernelMultiGilStat",
			 # simulation parameters
			 hopColIndex = as.integer(stratHopSkipJump$hopMat@colindices-1),
			 hopRowPointer = as.integer(stratHopSkipJump$hopMat@rowpointers-1), 
			 skipColIndex = as.integer(stratHopSkipJump$skipMat@colindices-1),
			 skipRowPointer = as.integer(stratHopSkipJump$skipMat@rowpointers-1),
			 jumpColIndex = as.integer(stratHopSkipJump$jumpMat@colindices-1),
			 jumpRowPointer = as.integer(stratHopSkipJump$jumpMat@rowpointers-1),
			 rateHopInMove = as.numeric(rateHopInMove), 
			 rateSkipInMove = as.numeric(rateSkipInMove), 
			 rateJumpInMove = as.numeric(rateJumpInMove), 
			 blockIndex = as.integer(blockIndex),
			 simul = as.integer(simul),
			 infested = as.integer(infested),
			 infestedDens = as.numeric(infestedDens),
			 endIndex = as.integer(length(infestH) - 1),
			 L = as.integer(L),
			 endTime = as.numeric(endTime),
			 indexInfest = as.integer(indexInfest),
			 timeI = as.numeric(timeI),
			 rateMove = as.numeric(rateMove),
			 seed = as.integer(seed),
			 Nrep = as.integer(Nrep),
			 # stats
			 getStats=as.integer(getStats),
			 nbins = as.integer(length(breaksGenVar)),
			 cbin = as.integer(cbin),
			 cbinas = as.integer(cbinas),
			 cbinsb = as.integer(cbinsb),
			 indices = as.integer(dist_indices),
			 statsTable = as.numeric(statsTable),
			 nbStats = as.integer(nbStats),
             		 sizeVvar = as.integer(sizeVvar)
			 )

		out$infestedDens<-out$infestedDens/Nrep;
	
		# make matrix out of statsTable
		out$statsTable<-matrix(out$statsTable,byrow=FALSE,ncol=Nrep)
		
		# remove interblock and intrablock stats
		# keep:
		# general semivar + stdev
		# sameblock - acrossstreets semivar + stdev
		# inf.house, inf.block, and inf.house/inf.block count
		keepable<-c(1:(2*length(cbin)), 6*length(cbin)+1:(2*length(cbin)), sizeVvar+(1:3))
		# clean away NANs introduced in C
		notNAN <- which(!is.nan(out$statsTable[, 1]))
		
		keep<-intersect(notNAN,keepable)
		out$statsTable<-out$statsTable[keep, ]

		infestH <- out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1
		out$indexInfest <- infestH
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]
		out$timeI <- timeH

		return(out)
	}

    multiThetaMultiGilStat<-function(cumulProbMat, blockIndex, infestH, timeH, endTime, rateMove, Nrep, coords, breaksGenVar, seed=1, simul=TRUE, getStats=TRUE, halfDistJ = -1, halfDistH = -1, useDelta = -1, delta = -1, rateHopInMove = -1, rateSkipInMove = -1, rateJumpInMove = -1, breaksStreetVar = breaksGenVar, dist_out = NULL){
		
		# seed <- runif(1, 1, 2^31-1)
		#for random seeding of stochastic simulation	

	  	L<-dim(coords)[1]
		indexInfest <- rep(-1, L)
		timeI <- rep(-1, L)
		infested <- rep(0, L)
		indexInfest[1:length(infestH)] <- infestH - 1
		timeI[1:length(timeH)] <- timeH
		infested[infestH] <- 1
		infestedDens<-rep(0,length(infested))

		if(is.null(dist_out))	
		dist_out <- makeDistClassesWithStreets(as.vector(coords[, 1]), as.vector(coords[, 2]), breaksGenVar, blockIndex)

		dist_mat <- dist_out$dists
		dist_indices <- dist_out$CClassIndex
		cbin <- dist_out$classSize
		cbinas <- dist_out$classSizeAS
		cbinsb <- dist_out$classSizeSB

		# if cumulProbMat is not passed, create blank cumulProbMat for C computation
		# set useProbMat to FALSE	
		if(is.null(cumulProbMat)){
			cumulProbMat <- mat.or.vec(L, L)
			useProbMat <- FALSE
		}else{ # else pass dummy dist_mat so as not to take up memory space	
			dist_mat <- -1
            		useProbMat <- TRUE
		}

		# stats selection
		# need to implement system where we can add and remove stats
		###===================================
		## CURRENT STATS:
		## General Semivariance
		## General Semivariance Std. Dev.
		## Interblock Semivariance
		## Interblock Semivariance Std. Dev.
		## Intrablock Semivariance
		## Intrablock Semivariance Std. Dev.
		## By block Semivariance
		## By block Semivariance Std. Dev.
		##	= 8 * length(cbin)
		## Number Infested Houses
		## Number Infested Blocks
		## (Infested Houses)/(Infested Blocks)
		##	= 8 * length(cbin) + 3
		###===================================
       		sizeVvar<-8*length(cbin)
		nbStats<- sizeVvar + 3
		statsTable<-mat.or.vec(nbStats,Nrep)

		out<- .C("multiThetaMultiGilStat",
			 # simulation parameters
			 probMat = as.numeric(cumulProbMat),
			 distMat = as.numeric(dist_mat),
			 blockIndex = as.integer(blockIndex),
			 simul = as.integer(simul),
			 infested = as.integer(infested),
			 infestedDens = as.numeric(infestedDens),
			 endIndex = as.integer(length(infestH) - 1),
			 L = as.integer(L),
			 endTime = as.numeric(endTime),
			 indexInfest = as.integer(indexInfest),
			 timeI = as.numeric(timeI),
			 rateMove = as.numeric(rateMove),
			 sizeRateMove = as.integer(length(rateMove)),
			 seed = as.integer(seed),
			 Nrep = as.integer(Nrep),
			 # stats
			 getStats=as.integer(getStats),
			 nbins = as.integer(length(breaksGenVar)),
			 cbin = as.integer(cbin),
			 cbinas = as.integer(cbinas),
			 cbinsb = as.integer(cbinsb),
			 indices = as.integer(dist_indices),
			 sizeVvar = as.integer(sizeVvar),
			 statsTable = as.numeric(statsTable),
			 nbStats = as.integer(nbStats)
			 )

		out$infestedDens<-out$infestedDens/(Nrep*length(rateMove));
		out$statsTable<-matrix(out$statsTable,byrow=FALSE,ncol=Nrep)
		good <- which(!is.nan(out$statsTable[, 1]))
		out$statsTable <- out$statsTable[good, ]
		#statsTable now includes 7 stats
		#LDsynLikSpat doesn't account for this, make statsTable only include first two stats for backward compatibility
		#out$statsTable<-out$statsTable[1:(2*(length(genIntervals)-1)), ]
		infestH <- out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1
		out$indexInfest <- infestH
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]
		out$timeI <- timeH

		return(out)
		# return(list(infestH, timeI))
	}

	getPosteriorMapsLD<-function(Fit,Data,repByTheta=1){
		Data$Nrep<-repByTheta
		meanMap<-0*rep(0,dim(Data$maps)[1])
		thetas<-as.matrix(Fit$Posterior2)
		nbThetas<-dim(thetas)[1]
		for(numTheta in 1:nbThetas){
			theta<-thetas[numTheta,]
			ModelOutBest<-Model(theta,Data,postDraw=TRUE)
			meanMap<-meanMap+ModelOutBest$yhat
		}
		
		meanMap<-meanMap/nbThetas
		attributes(meanMap)$nbThetas<-nbThetas
		attributes(meanMap)$repByTheta<-repByTheta

		return(meanMap)
	}

	generate_prob_mat_C <- function(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, dist_mat, blockIndex, L=sqrt(length(dist_mat)), cumul=FALSE )
	{
		prob_mat <- mat.or.vec(L, L)
		out <- .C("generateProbMat",
			 halfDistJ = as.numeric(halfDistJ),
			 halfDistH = as.numeric(halfDistH), 
			 useDelta = as.integer(useDelta), 
			 delta = as.numeric(delta), 
			 rateHopInMove = as.numeric(rateHopInMove), 
			 rateSkipInMove = as.numeric(rateSkipInMove), 
			 rateJumpInMove = as.numeric(rateJumpInMove), 
			 dist_mat = as.numeric(dist_mat),
			 prob_mat = as.numeric(prob_mat),
			 blockIndex = as.integer(blockIndex),
			 cumul = as.integer(cumul),
			 L = as.integer(L))
		if(cumul){
			byrow<-FALSE
		}else{
			byrow<-TRUE
		}

		return(matrix(out$prob_mat, L, L, byrow = byrow))
	}

	generate_prob_mat_Grant_C <- function(limitHopSkips,  weightSkip, weightJump, dist_mat, blockIndex, L=sqrt(length(dist_mat)), cumul=TRUE){
		prob_mat <- mat.or.vec(L, L)
		out <- .C("generateProbMatGrant",DUP=FALSE,NAOK=TRUE,
			 limitHopSkips = as.numeric(limitHopSkips),
			 weightSkip = as.numeric(weightSkip), 
			 weightJump = as.numeric(weightJump), 
			 dist_mat = as.numeric(dist_mat),
			 prob_mat = as.numeric(prob_mat),
			 blockIndex = as.integer(blockIndex),
			 L = as.integer(L),
			 cumul = as.integer(cumul)
			 )
		if(cumul){
			byrow<-FALSE
		}else{
			byrow<-TRUE
		}

		return(matrix(out$prob_mat, L, L, byrow = byrow))
	}
}



# take the out from gillespie and transform it in normal maps
# (create a binary list INFEST where 1 = infested, 0 = uninfested)
infestSerieToMaps<-function(outGillespie, sp, mapTime=outGillespie$time){
	
	infestH <- unlist(outGillespie[1])
	timeH <- unlist(outGillespie[2])

	# reset times above mapTime
	indLastInf<-max(which(timeH<=mapTime))
	timeH<-timeH[1:indLastInf]
	infestH<-infestH[1:indLastInf]

	# make the map
	map<-as.data.frame(sp)
	names(map)<-c("X","Y")
	map$ages <- map$infest<-rep(0, length(sp[,1]))

	map$infest[infestH] <- 1

	map$ages[infestH]<-mapTime-timeH

	return(map)
}

updatePredict<-function(Fit,Data,infestHints,repByTheta=1000){
  # make a lot of simulations and select the ones with same house
  # initmap 
  # or same block (to add: ,sameBlockApprox=FALSE)
      Data$Nrep<-repByTheta
      meanMap<-0*rep(0,dim(Data$maps)[1])
      thetas<-as.matrix(Fit$Posterior2)
      nbThetas<-dim(thetas)[1]
      for(numTheta in 1:nbThetas){
	theta<-thetas[numTheta,]
	ModelOutBest<-Model(theta,Data,postDraw=TRUE,infestHints)
	meanMap<-meanMap+ModelOutBest$yhat
      }
      meanMap<-meanMap/nbThetas
      attributes(meanMap)$nbThetas<-nbThetas
      attributes(meanMap)$repByTheta<-repByTheta

  outBase <- multiGilStat(cumulProbMat=cumulProbMat, blockIndex, infestH, timeH=rep(-1,length(infestH)), endTime = nbit, rateMove, Nrep, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE)

      return(meanMap)
    }

simulObserved<-function(infested,openRate,detectRate){
  # initial state
  observed<-infested
  nbHouses<-length(infested)
  cat("Total:", nbHouses,"; Infested:",sum(infested));

  # only openRate of houses observed 
  nbHousesNonObserved<-rpois(n=1,lambda=nbHouses*(1-openRate))
  observed[sample(1:nbHouses,nbHousesNonObserved)]<-0
  cat("; Openned Inf:", sum(observed));

  # only XX% of infested are observed infested
  nbNonDetected<-rpois(n=1,lambda=sum(observed)*(1-detectRate))
  nonDetected<-sample(which(observed==1),nbNonDetected)
  observed[nonDetected]<-0
  cat("; Obs Inf:", sum(observed),"\n");

  return(observed)
}

getBestPredict<-function(maps,cumulProbMat,init,nbit=52*2,Nrep=10000){
  infestH<-which(init==1)
  outBase <- multiGilStat(cumulProbMat=cumulProbMat, blockIndex=maps$blockIndex, infestH, timeH=rep(-1,length(infestH)), endTime = nbit, rateMove, Nrep, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE,getStats=FALSE)
  return(outBase$infestedDens)
}

getPosteriorMaps<-function(maps,thetas,cumulProbMat=NULL,initInf,nbit,repByTheta=10){

  thetas<-as.matrix(thetas)
  nbThetas<-dim(thetas)[1]

  infestH<-which(initInf==1)
  infestedDens <- multiThetaMultiGilStat(cumulProbMat=cumulProbMat, blockIndex=maps$blockIndex, infestH, timeH=rep(-1,length(infestH)), endTime = nbit, thetas, repByTheta, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE,getStats=FALSE)$infestedDens

  return(infestedDens)
}


#===========================
# simple model: regression principles
#===========================

# use times 1 year and 2 years to 
# - set beta of the link 
# (distance to nearest infested at t)-> (infestation  at t+1 years)
# (- beta for number of streets)
# - make a prediction for 3 years 

# in addition of the problems to extrapolates this thing at different years
# we expect to show that this is not a good predictor of the possible 
# => set != seeds at t=3years, and evaluate the quality of prediction of the possibles
#    then do the same for the wood approach and hopefully will be better
#    the problem is the criterium, so look at: mean quality, sd quality


#============================
# simplest model: glm with closest infested neighbor at t-1
#============================
# get the nearest neighbor of coordsOrig in coordsNeigh
getNearestNeighDist<-function(coordsOrig,coordsNeigh=coordsOrig,tr=NULL){
	if(is.null(tr)){ # set tr to the maximum possible distance
		xrange<-range(c(coordsOrig$X,coordsNeigh$X))
		yrange<-range(c(coordsOrig$Y,coordsNeigh$Y))
		xextent<-xrange[2]-xrange[1]
		yextent<-yrange[2]-yrange[1]
		tr<-sqrt(xextent^2+yextent^2)
	}
	distToInf<-nearest.dist(coordsOrig,coordsNeigh,delta=tr)
	distToInf<-as.matrix(distToInf)
	nearestNeigh<-apply(distToInf,1,min)
	return(nearestNeigh)
}

# get the basic model
getBasicModel<-function(maps,infest1Name,infest2Name){
  maps$nearestInf<-getNearestNeighDist(maps[,c("X","Y")],maps[which(maps[,infest1Name]==1),c("X","Y")])

  basicModel<-glm(maps[,infest2Name]~nearestInf,data=maps,family=binomial())
}
# make prediction with basic model
predictBasicModel<-function(maps,infestInitName,model=basicModel,timePredictOverTimeTraining=1){
  maps$nearestInf<-getNearestNeighDist(maps[,c("X","Y")],maps[which(maps[,infestInitName]==1),c("X","Y")])/timePredictOverTimeTraining
  predicted<-predict(model,newdata=maps,type="response")
  return(predicted)
}

# Tests
# test_file("test-functions_migration.R")

