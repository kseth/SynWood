# functions for the migration/spatial stats
# data for the tests
library(geoR)
library(testthat)
# source("RanalysisFunctions.R")

# compile the .c if needed
compilLoad<-function(sourcef,options=""){
	try(file.remove(gsub(".c$",".o",sourcef)),silent=TRUE)
	if(file.exists(sourcef)){
		exitCode<-system(paste("R CMD SHLIB",options,sourcef))
		if(exitCode!=0){
			stop("Compilation of ",sourcef," failed")
		}else{
			importOk<-try(dyn.load(gsub(".c$",".so",sourcef)),silent=TRUE)
		}
	}else{
		importOk<-paste(sourcef,"missing")
		class(importOk)<-"try-error"

	}
	return(importOk)
}

#======================
# Basic functions
#======================
# dist function
makeDistMat<-function(xs,ys){
	expect_equal(length(xs),length(ys))
	dists<-mat.or.vec(length(xs),length(xs))
	out<-.C("makeDistMat",
	   xc = as.numeric(xs),
	   L = as.integer(length(xs)),
	   yc = as.numeric(ys),
	   dists = as.numeric(dists)
	   )
	dists<-matrix(out$dists,length(xs))

	return(dists)
}

# create an infestH/indexInfest vector:
# from a vector with the number of positive unit per location
# for each positive unit gives the index of the location 
# assumes not specific order for infestation
MakeIndexFromNinfest <- function(nInfs){
  indInf <- rep(0,sum(nInfs))
  b<-1
  for(i in 1:length(nInfs)){
    if(nInfs[i]>0){
      e<- b+nInfs[i]-1
      indInf[b:e]<-i
      b<-e+1
    }
  }
  return(indInf)
}

MakeNinfestFromIndex <- function(ind,nMacro){
  nInfs<-rep(0,nMacro)
  nInfs[as.numeric(names(table(ind)))] <- table(ind)
  return(nInfs)
}



#======================
# Generate hop/skip/jump matrixes
#======================

##
# returns matrixes hopMat, skipMat, jumpMat
# with 1 if can hop/skip/jump there
# only need to pass the @colindices and @rowpointers to the c methods
# coords - X, Y coords
# limitHopSkip - the limit for local movement
# limitJump - the upper limit for long distance movement
# lowerLimitJump - the lower limit for long distance movement
# blockIndex - should blocks be used to make skips (upper distance class), pass blockIndices (same length as coords)
# lowerLimitSkip - if no blocks but still want skips, pass a lower limit for skips
generate_stratified_mat <- function(coords, limitHopSkip, limitJump, lowerLimitJump=0, blockIndex=NULL, lowerLimitSkip=NULL,autLocalHops=FALSE)
{

	###=====================
       	## make the jump  matrix
	###=====================
	# upper limit for jumps
	dist_mat_jump <- nearest.dist(coords, y=NULL, method = "euclidian", delta = limitJump, upper = NULL)
	dist_mat_jump <- cleanup(dist_mat_jump)

	# if we also have a lower limit for the jumps
	if(lowerLimitJump > 0)
	{
		dist_mat_jump_low <- nearest.dist(coords, y=NULL, method = "euclidian", delta = lowerLimitJump, upper = NULL)
		dist_mat_jump_low <- cleanup(dist_mat_jump_low)

		dist_mat_jump <- dist_mat_jump - dist_mat_jump_low #subtract lower values from total
		dist_mat_jump <- as.spam(dist_mat_jump) #respam the matrix
	}
	
	#define jumpMat
	jumpMat <- dist_mat_jump
	jumpMat@entries <- rep(1, length(jumpMat@entries))
	jumpMat <- as.spam(jumpMat)

	###=====================
	## make the hop/skip matrices
	###=====================
	# upper limit for hops + skips
	dist_mat_hop_skip <- nearest.dist(coords, y=NULL, method = "euclidian", delta = limitHopSkip, upper = NULL)
	dist_mat_hop_skip <- cleanup(dist_mat_hop_skip)
	
	#define hopMat (definition modified if blockIndices or lowerLimitSkip passed)
	hopMat <- dist_mat_hop_skip
	hopMat@entries <- rep(1, length(hopMat@entries))
	hopMat <- as.spam(hopMat)

	if(!is.null(blockIndex)){
		# if blockIndices are passed
		# make SameBlock (SB) matrix of households in same block
		# define skipMat as across blocks
		# hopMat as within blocks

		SB <- nearest.dist(x=cbind(blockIndex,rep(0,length(blockIndex))), method="euclidian", upper=NULL, delta=0.1)
		SB@entries <- rep(1,length(SB@entries))
		SB<-as.spam(SB)

		skipMat <- dist_mat_hop_skip
		skipMat@entries <- rep(1, length(skipMat@entries))

		hopMat <- hopMat * SB
		hopMat <- as.spam(hopMat)

		skipMat <- skipMat - hopMat
		skipMat <- as.spam(skipMat)

	}else if(!is.null(lowerLimitSkip)){
			# if lowerLimitSkip is passed
			# skip from lowerLimitSkip to limitHopSkip
			# hop from 0 to lowerLimitSkip

			dist_mat_hopskip_low <- nearest.dist(coords, y=NULL, method = "euclidian", delta = lowerLimitSkip, upper = NULL)
			dist_mat_hopskip_low <- cleanup(dist_mat_hopskip_low)

			hopMat <- dist_mat_hopskip_low
			hopMat@entries <- rep(1, length(hopMat@entries))
			hopMat <- as.spam(hopMat)

			skipMat <- dist_mat_hop_skip
			skipMat@entries <- rep(1, length(skipMat@entries))
			skipMat <- skipMat - hopMat
			skipMat <- as.spam(skipMat)
		}
	if(autLocalHops){
	  diag(hopMat) <- 1
	}

	if(!exists("skipMat"))
		return(list(hopMat = hopMat, jumpMat = jumpMat))
	else
		return(list(hopMat = hopMat, skipMat = skipMat, jumpMat = jumpMat))

}

# visualisation of 
VisuHSJ<-function(maps,stratHopSkipJump,house){
  plot(maps$X,maps$Y,asp=1,pch=".") #all
  with(maps[house,],points(X,Y,col="blue")) # house
  with(maps[which(stratHopSkipJump$hopMat[,house]==1),],
       points(X,Y,col="red")) # hops
  with(maps[which(stratHopSkipJump$skipMat[,house]==1),],
       points(X,Y,col="green")) # skip
  with(maps[which(stratHopSkipJump$jumpMat[,house]==1),],
       points(X,Y,col="purple")) # jump
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

	if(cumul){
		#redefine probMat as the cumulative sum on each line
		for(col in 1:L){
			probMat[,col]<-cumsum(probMat[,col])
		}
	}

	return(probMat)
}


#=============================
# Simple grid definition
#=============================
# num.row, num.col - number rows, number cols
# row.dist, col.dist - distance between adjacent rows, adjacent columns
# house 1 is bottom left, go up vertically, then go to the next column
# house num.rows * num.cols is at top right
makeGrid <- function(num.rows, num.cols, row.dist, col.dist = row.dist){

	rows <- (1:(num.rows*num.cols) - 1)
	rows <- floor(rows/num.rows)
	rows <- rows * row.dist
	cols <- (1:(num.rows*num.cols) - 1)
	cols <- floor(cols%%num.rows)
	cols <- cols*col.dist
	
	maps <- data.frame(X = rows, Y = cols)

	return(maps)	
}

#==========================
# Partition of map (with grid or kmeans)
#==========================
## divides the map according to the partitionSizes and the type of divide under consideration
divideMap <- function(maps, partitionSizes, typeDivide = "grid", language = "R"){

	map.partitions <- list()
	length(map.partitions) <- length(partitionSizes) #number different grid partitions will be used
	
	if(typeDivide == "grid"){ ## if grid

		for(part in 1:length(partitionSizes))
			map.partitions[[part]] <- partitionMap(maps$X, maps$Y, partitionSizes[part]) 
	
	}else if(typeDivide == "kmeans"){ ## if kmeans

		for(part in 1:length(partitionSizes))
			map.partitions[[part]] <- partitionKMeans(maps$X, maps$Y, partitionSizes[part], language = language)

	}

	return(map.partitions)
}	

# partition.rows == number of rows in final grid
# partition.cols == number of cols in final grid
# if want 1 by 1 final grid pass 1, 1
partitionMap <- function(X, Y, partition.rows, partition.cols = partition.rows){

	maxX <- max(X)
	minX <- min(X)
	maxY <- max(Y)
	minY <- min(Y)

	yBreaks <- seq(from = minY, to = maxY, length.out = partition.rows+1)

	uniqY <- unique(Y)
	outY <- unlist(lapply(uniqY, function(insert, breaks){
								breaks <- (breaks > insert)
								comparison <- xor(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)])
								if(any(comparison))
									return(which(comparison))
								else
									return(length(breaks)-1)
							}, breaks = yBreaks))

	matchY <- match(Y, uniqY)
	indexY <- outY[matchY]
	
	xBreaks <- seq(from = minX, to = maxX, length.out = partition.cols+1)
	uniqX <- unique(X)
	outX <- unlist(lapply(uniqX, function(insert, breaks){
								breaks <- (breaks > insert)
								comparison <- xor(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)])
								if(any(comparison))
									return(which(comparison))
								else
									return(length(breaks)-1)
							}, breaks = xBreaks))

	matchX <- match(X, uniqX)
	indexX <- outX[matchX]

	## want all unique pairs of indexX, indexY -> use cantor's pairing function
	combXY <- 1/2*(indexX+indexY)*(indexX+indexY+1) + indexY
	uniqXY <- unique(combXY)
	
	## the grid index of each of the final outputs
	indexXY	<- match(combXY, uniqXY)

	## plot the output
	## colors <- rainbow(length(uniqXY))
	## plot(X, Y, col = colors[indexXY])

	## calculate the number of houses in each cell
	minIndex <- min(indexXY)
	maxIndex <- max(indexXY)
	housesPerCell <- unlist(lapply(minIndex:maxIndex, function(x, indexes){ return(length(which(indexes == x))) }, indexes = indexXY))

	return(list(index = indexXY, num_cells = length(uniqXY),  housesPerCell = housesPerCell))

}

importkmeans_Ok <- try(dyn.load("kmeans.so"), silent=TRUE)
if(class(importkmeans_Ok)=="try-error"){
	importkmeans_Ok <-compilLoad("kmeans.c")
	
}
## X the x coords of the points
## Y the y coords of the points
## num_clusts - the number of clusters to make
## language (whether to do it in R or in C) 
##	 C better for smaller clusters, R more stable for larger clusters
partitionKMeans <- function(X, Y, num_clusts, language = "R"){
	
	if(language == "C" && class(importkmeans_Ok) != "try-error"){

		cat("using C kmeans partitioning\n")

		# XY positions converted to c format
		x <- as.vector((matrix(c(X, Y), ncol = 2))) 
		# choose random centers
		# converted to c format
		clust_centers <- runif(num_clusts, 1, length(X)) # choose random cluster centers
		d <- as.vector((matrix(c(X[clust_centers], Y[clust_centers]), ncol = 2))) 	

		# stores final deviation per cluster
		dev <- rep(0, num_clusts)

		# stores index of each XY at the end
		b <- rep(0, length(X))

		# workspace
		f <- rep(0, length(X))

		# stores observations per cluster
		e <- rep(0, num_clusts)

		# number of observations 
		i <- length(X)

		# number of dimensions
		j <- 2

		# minimum size of each cluster
		nz <- floor(i/num_clusts)*0.85
		out <- .C("clustr",
	   	  		 x=as.numeric(x), 
		  		 d=as.numeric(d),
		  		 dev=as.numeric(dev),
		  		 b=as.integer(b),
		  		 f=as.numeric(f),
		  		 e=as.integer(e),
		  		 i=as.integer(i),
		  		 j=as.integer(j),
		  		 n=as.integer(num_clusts),
		  		 nz=as.integer(nz),
		  		 k=as.integer(num_clusts))

		housesPerCell <- out$e

		while(any(housesPerCell == 0)){ ## if the partitioning messed up, redo!

			out <- .C("clustr",
	   	  			 x=as.numeric(x), 
		  			 d=as.numeric(d),
			  		 dev=as.numeric(dev),
			  		 b=as.integer(b),
			  		 f=as.numeric(f),
			  		 e=as.integer(e),
			  		 i=as.integer(i),
			  		 j=as.integer(j),
			  		 n=as.integer(num_clusts),
			  		 nz=as.integer(nz*.85),
			  		 k=as.integer(num_clusts))
     
			housesPerCell <- out$e
		}
     
		index <- out$b
		num_cells <- num_clusts

		## plot the output
		## colors <- rainbow(num_clusts)
		## plot(X, Y, col = colors[index])

		return(list(index=index, num_cells=num_cells, housesPerCell=housesPerCell))

	} else {
		# cat("using R kmeans partitioning\n")

		points <- matrix(c(X, Y), ncol = 2)
		out <- kmeans(points, centers = num_clusts, iter = 100, nstart = 2)

		index <- out$cluster
		num_cells <- num_clusts
		housesPerCell <- out$size
				
		## plot the output
		## colors <- rainbow(num_clusts)
		## plot(X, Y, col = colors[index])

		return(list(index=index, num_cells=num_cells, housesPerCell=housesPerCell))
	}		
}


#=============================
# Concentric Circles
#=============================

conc.circles <- function(X, Y, distClasses, initInfested){

	out <- makeDistClasses(X, Y, distClasses)

	# since out$CClassIndex is triangular and the diagonals are 0
	all.indexes <- matrix(out$CClassIndex, byrow = TRUE, nrow = length(X)) + matrix(out$CClassIndex, byrow = FALSE, nrow = length(X))


	if(distClasses[1] != 0) #if dont want to include self in concentric circles
		diag(all.indexes) <- -1

	# return only the indexes of the initInfested
	all.indexes <- all.indexes[initInfested, ]

	# count how many values exist per distance class per house to get prevalence
	prev.counts <- mat.or.vec(length(initInfested), length(distClasses)-1)
	for(house in 1:length(initInfested))
		for(class in 0:(length(distClasses)-2)){

			howmany <- length(which(all.indexes[house, ] == class))
			prev.counts[house, class+1] <- howmany
		} 

	return(list(circleIndexes = all.indexes, counts = prev.counts))

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

#==========================
## Functions defined in C
## all the main loop + statistical functions functions
#==========================

# import C functions if possible
importOk<-try(dyn.load("functions_migration.so"), silent=FALSE)
if(class(importOk)=="try-error"){
	importOk <-compilLoad("functions_migration.c","-lgsl -lgslcblas -lm samlmu.f")
	
}
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
	multiGilStat<-function(cumulProbMat, blockIndex, infestH, timeH, endTime, rateMove, Nrep, coords, breaksGenVar, seed=1, simul=TRUE, getStats=TRUE, halfDistJ = -1, halfDistH = -1, useDelta = -1, delta = -1, rateHopInMove = -1, rateSkipInMove = -1, rateJumpInMove = -1, dist_out = NULL){
		
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
	
		# if cumulProbMat is not passed, create blank cumulProbMat for C computation
		# set useProbMat to FALSE	
		if(is.null(cumulProbMat)){
			cumulProbMat <- mat.or.vec(L, L)
			useProbMat <- FALSE
		}else{ # else pass dummy dist_mat so as not to take up memory space	
			dist_mat <- 0
            		useProbMat <- TRUE
		}

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

			dist_indices <- dist_out$CClassIndex
			cbin <- dist_out$classSize
			cbinas <- dist_out$classSizeAS
			cbinsb <- dist_out$classSizeSB
	
		}else{
 	      		sizeVvar<-0
			nbStats<-0
			statsTable<-0

			dist_indices <- 0
			cbin <- 0
			cbinas <- 0
			cbinsb <- 0
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
		
		keep<-intersect(notNAN)
		out$statsTable<-out$statsTable[keep, ]

		infestH <- out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1
		out$indexInfest <- infestH
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]
		out$timeI <- timeH

		return(out)
	}

	#' @title Distance stratified dispersion simulator 
	#'
	#' @description #to be added
	#' @param stratHopSkipJump hop skip jump move possibilities (list of hoppable, skippable, jumpable locations given as spam matrices, result of call to generate_stratified_mat)
	#' @param blockIndex the block index (or ID) for each household, pass NULL if want to use model w/o streets
	#' @param infestH the infested households (index numbers)
	#' @param timeH the infestation time of the households
	#' @param endTime duration of simulation (when to cutoff gillespie)
	#' @param rateMove the frequency of movement (number of infestations per unit of time per infested house, in same units as endTime)
	#' @param rateHopInMove the percentage of hops
	#' @param rateSkipInMove the percentage of skips
	#' @param rateJumpInMove the percentage of jumps
	#' @param Nrep the number of simulations
	#' @param coords the x and y coordinates of all households in the map
	#' @param seed the starting seed for the simulations (defaults to 1)
	#' @param simul should simulations be run (defaults to TRUE), pass FALSE if want only to calculate statistics on infestH
	#' @param maxInfest the maximum number of infestations for each of the households (defaults to 1 per household)
	#' @param getStats should statistics be obtained from the simulations (defaults to TRUE), pass FALSE if only want simulation output
	#' @param dist_out a semivariance statistics object (defaults to NULL, see makeDistClassesWithStreets or makeDistClasses)
	#' @param map.partitions a partition statistics object (defaults to NULL, see partitionMap)
	#' @param conc.circs a circles statistics object (defaults to NULL, see conc.circles)
	#' @param atRisk.trs thresholds for atRisk statistics (defaults to NULL, breaks between houses, should be between 0 and infinity)
	#' @param atRisk.ncoefs number of coefficients to be calculated (defaults to NULL)
	#' @param typeStat which statistics to calculate (defaults to "semivariance" for backward compatibility), should be any or multiple of "semivariance", "grid", "circles", "atRisk", "num_inf" - appropriate distance objects should also be passed, see above
	#' @param whichPairwise which pairwise (semivariance) to calculate (defaults to c("semivariance", "moran", "geary", "ripley")), should be any or multiple of these
	#' @param detectRate whether to withold data (defaults to 1, detectRate = 0.7, 30% of data randomly withheld when generating statistics + counts)
       	#' @param rateIntro the rate of introductions (new infestations per unit time, same as endTime, defaults to 0)	
	#' @param nPartCoefs the number of partition regression coefficients to have per partition (defaults to 0)
	#' @param iPartLMoments the indices of the partition L-moments to have per partition (defaults to 1, 2, 3 or Variance, Skewness, Kurtosis)
	#' @return a list with:\itemize{
	#'         \item{ToBeAdded}{ The output least needs cleaning} 
	#'         }
	#' @details #to be added
	#' @author Corentin M. Barbu, Karthik Sethuraman
	#' @seealso #to be added
	#' @examples # to be added
	#'
	#' @export noKernelMultiGilStat
	noKernelMultiGilStat <- function(
	    	stratHopSkipJump, 
		blockIndex, 
		infestH, 
		timeH, 
		endTime, 
		rateMove, rateHopInMove, rateSkipInMove, rateJumpInMove, 
		Nrep, 
		coords, 
		seed=1, 
		simul=TRUE,
		maxInfest=rep(1,dim(coords)[1]), 
		getStats=TRUE, 
		dist_out = NULL, map.partitions = NULL, conc.circs = NULL, atRisk.trs = NULL, atRisk.ncoefs = NULL,
		typeStat = "semivariance",
	        whichPairwise = c("semivariance", "moran", "geary", "ripley"),	
		detectRate = 1, rateIntro = 0,
		nPartCoefs = 0, iPartLMoments = c(1, 2)){

		# do we have blocks?
		haveBlocks <- (!missing(blockIndex) && !is.null(blockIndex))

		#=====================
	        # initialize all necessary variables
		#=====================

	  	L<-dim(coords)[1]
		indexInfest <- rep(-1, sum(maxInfest)) 
		timeI <- rep(-1, sum(maxInfest))
		infested <- rep(0, L)
		infestedDens<-rep(0,L)
		endIndex <- as.integer(length(infestH) - 1)
		# cat("tot maxInfest:",sum(maxInfest),"\n");
		if(length(infestH)>0){
		  indexInfest[1:length(infestH)] <- infestH - 1
		  timeI[1:length(timeH)] <- timeH
		  infested <- MakeNinfestFromIndex(infestH,L)
		}else{
		  infestH <- -1 # avoid to pass a NULL to .C
		}

		# if stratHopSkipJump, throw error 	
		if(is.null(stratHopSkipJump)){
			stop("need to pass a stratified hop/skip/jump matrix; see generate_stratifed_mat")
		}

		# implemented stats
		implStats <- c("semivariance", "grid", "circles", "atRisk", "num_inf")

		# initialize all the statistics to 0
		matchStats <- 0

		# if getStats and specific statistics are used, then change their value
		# semivariance statistics
		nbins <- 0
		dist_indices <- 0
		cbin <- 0
		cbinas <- 0
		cbinsb <- 0
		semivar.nbStats <- 0
		semivar.statsTable <- 0

		# grid/partition statistics
		numDiffGrids <- 0
		gridIndexes <- 0
		gridNumCells <- 0
		gridCountCells <- 0
		grid.nbStats <- 0
		grid.numCoeffs <- 0 
		grid.numLmoments <- 0
		grid.statsTable <- 0

		#circle statistics
		numDiffCircles <- 0
		numDiffCenters <- 0
		circleIndexes <- 0
		circleCounts <- 0
		circle.nbStats <- 0
		circle.statsTable <- 0

		# num_inf statistics
		inf.nbStats <- 0
		inf.statsTable <- 0

		# at_risk statistics
		atRisk.nbCoefs<-0
		atRisk.statsTable<-0

		if(getStats){

			## pass the corresponding matched numbers to C instead of the name of the statistics
			matchStats <- match(typeStat, implStats)

			if(any(is.na(matchStats))){ #throw an error regarding stats not yet implemented
				stop(paste0(typeStat[is.na(matchStats)], " not implemented! Only implemented ", implStats))
			}

			if("num_inf" %in% typeStat){
				# 1 if no blocks, 3 if blocks
				# stats selection
				##===============
				# Number Locations Infested (location or macrounits > 0)
				# Number Units infested (total microunits +)
				# Number Units Infested / Number Locations Infested
				# If haveBlocks:
				#	Number Blocks Infested
				#	Number Infested / Number Blocks Infested
				inf.nbStats <- 3 + haveBlocks*2
				inf.statsTable <- mat.or.vec(inf.nbStats, Nrep)
			}	

			# if want to calculate semivariance stats
			if("semivariance" %in% typeStat){
				if(is.null(dist_out)){

					stop("no semivariance dist_out object passed, cannot execute!")
				}

				dist_indices <- dist_out$CClassIndex
				cbin <- dist_out$classSize
				nbins <- length(cbin) + 1 #number of breaks not bins

				# stats selection
				###===================================
				## CURRENT STATS:
				## General Semivariance (new - new)
				## Moran's I
				## Geary's C
				## Ripley's L
				##	= 4 * length(cbin)
				## if haveBlocks
				##	By block Semivariance (same block - across streets)
				##	By block Semivariance Std. Dev
				## 	= 6 * length(cbin)
				## DEPRECATED STATS:
				## General Semivariance Std. Dev. (new - new)
				## General Semivariance (old - new)
				## General Semivariance Std. Dev (old - new)
				###===================================
				if(!haveBlocks)
					semivar.nbStats <- 4*length(cbin)
				else{
					semivar.nbStats <- 6*length(cbin)
					cbinas <- dist_out$classSizeAS
					cbinsb <- dist_out$classSizeSB
				}

				semivar.statsTable<-mat.or.vec(semivar.nbStats,Nrep)
			}

			if("grid" %in% typeStat){

				if(is.null(map.partitions)){ # if an indexing of the map hasn't yet been passed, throw error
					stop("map.partitions passed as null, cannot execute!")
				}
				
				## the number of different indexings present in gridIndexes
				numDiffGrids <- length(map.partitions)

				## unlist the first list
				map.partitions <- unlist(map.partitions, recursive = FALSE)
		
				## define the indexing systems
				## every house should have numDiffGrids indexes
				gridIndexes <- map.partitions[which(names(map.partitions) %in% "index")]
				gridIndexes <- unlist(gridIndexes)
				names(gridIndexes) <- NULL

				## the length of each indexing - the number of cells per grid/index system
				gridNumCells <- map.partitions[which(names(map.partitions) %in% "num_cells")] 
				gridNumCells <- unlist(gridNumCells)
				names(gridNumCells) <- NULL

				# countCells will keep track of total number of houses per cell
				gridCountCells <- map.partitions[which(names(map.partitions) %in% "housesPerCell")]
				gridCountCells <- unlist(gridCountCells)
				names(gridCountCells) <- NULL

				# stats selection
				###===================================
				## CURRENT STATS 
				## (by grid system):
				## Variance of % positive per cell
				## Number Cells with at least 1 positive 
				## Fit quantile distribution to polynomial
				## a + bx + cx^2 + dx^3 + ... (grid.numCoeffs stats)
			        ##      = 2*numDiffGrids + 2*sum(grid.numCoeffs)
			       	## L-moment statistics (taken from quantile distribution)
			        ## (2nd, 3rd, 4th L-moments, L-scale, L-skewness, L-kurtosis)
				## L-mean should be ~ to median (also to num_inf)	
				###===================================
				grid.numCoeffs <- rep(max(1, nPartCoefs), length(gridNumCells)) #take max of nPartCoefs, 1 - need to pass in dummy object
				grid.numLmoments <- max(1, max(iPartLMoments)) #take the max of the indices 
				grid.nbStats <- 2*numDiffGrids + sum(grid.numCoeffs) + grid.numLmoments*numDiffGrids	
				grid.statsTable <- mat.or.vec(grid.nbStats, Nrep)
			}

			if("circles" %in% typeStat){
				if(is.null(conc.circs)){ # if an indexing of concentric circles hasn't yet been passed, throw error
					stop("conc.circles passed as null, cannot execute!")
				}

				numDiffCenters <- dim(conc.circs$counts)[1]
				numDiffCircles <- dim(conc.circs$counts)[2]
				circleIndexes <- t(conc.circs$circleIndexes) 
				circleCounts <- t(conc.circs$counts) 

				# stats selection
				###===================================
				## CURRENT STATS 
				## (by numDiffCircles):
				## Variance of % positive (across initInfested)
				## Mean of % positive (across initInfested) 
				##	= 2 * numDiffCircles
				###===================================
				circle.nbStats <- 2*numDiffCircles
				circle.statsTable <- mat.or.vec(circle.nbStats, Nrep)
			}
			
			if("atRisk" %in% typeStat){
				atRisk.nbCoefs<-atRisk.ncoefs
				atRisk.nbStats<-length(atRisk.trs)+atRisk.ncoefs
				atRisk.statsTable<-mat.or.vec(Nrep,atRisk.nbStats)
			}
		}

		#=============================
		# Simulate in C
		#=============================	

		skipColIndex <- if(!is.null(stratHopSkipJump$skipMat)){
		  as.integer(stratHopSkipJump$skipMat@colindices-1L)
		}else{
		  as.integer(-1L)
		}
		skipRowPointer <-if(!is.null(stratHopSkipJump$skipMat)){
		  as.integer(stratHopSkipJump$skipMat@rowpointers-1L)}else{
		    as.integer(-1L)}

		out<- .C("noKernelMultiGilStat",
			 # simulation parameters
			 hopColIndex = as.integer(stratHopSkipJump$hopMat@colindices-1L),
			 hopRowPointer = as.integer(stratHopSkipJump$hopMat@rowpointers-1L), 
			 skipColIndex = skipColIndex, # if no skips, pass dummy skip
			 skipRowPointer = skipRowPointer, # if no skips, pass dummy skip
			 jumpColIndex = as.integer(stratHopSkipJump$jumpMat@colindices-1L),
			 jumpRowPointer = as.integer(stratHopSkipJump$jumpMat@rowpointers-1L),
			 rateHopInMove = as.numeric(rateHopInMove), 
			 rateSkipInMove = as.numeric(rateSkipInMove), 
			 rateJumpInMove = as.numeric(rateJumpInMove), 
			 blockIndex = if(haveBlocks){as.integer(blockIndex)}else{as.integer(0)},
			 simul = as.integer(simul),
			 infested = as.integer(infested),
			 maxInfest = as.integer(maxInfest),
			 infestedDens = as.numeric(infestedDens),
			 endIndex = endIndex,
			 L = as.integer(L),
			 endTime = as.numeric(endTime),
			 indexInfest = as.integer(indexInfest),
			 timeI = as.numeric(timeI),
			 rateMove = as.numeric(rateMove),
			 rateIntro = as.numeric(rateIntro),
			 seed = as.integer(seed),
			 Nrep = as.integer(Nrep),
			 # stats
			 getStats=as.integer(getStats),
			 matchStats=as.integer(matchStats),
			 lengthStats=as.integer(length(matchStats)),
			 nbins = as.integer(nbins),
			 cbin = as.integer(cbin),
			 cbinas = as.integer(cbinas),
			 cbinsb = as.integer(cbinsb),
			 indices = as.integer(dist_indices),
			 semivar.statsTable = as.numeric(semivar.statsTable),
			 semivar.nbStats = as.integer(semivar.nbStats),
			 haveBlocks = as.integer(haveBlocks),
			 numDiffGrids = as.integer(numDiffGrids),
			 gridIndexes = as.integer(gridIndexes-1), #subtract 1 because C is 0 indexed
			 gridNumCells = as.integer(gridNumCells),
			 gridCountCells = as.integer(gridCountCells),
			 grid.nbStats = as.integer(grid.nbStats),
			 grid.numCoeffs = as.integer(grid.numCoeffs),
			 grid.numLmoments = as.integer(grid.numLmoments),
			 grid.statsTable = as.numeric(grid.statsTable),
			 numDiffCircles = as.integer(numDiffCircles),
			 numDiffCenters = as.integer(numDiffCenters),
			 circleIndexes = as.integer(circleIndexes), #already C indexed (was computed in C)
			 circleCounts = as.integer(circleCounts),
			 circle.nbStats = as.integer(circle.nbStats),
			 circle.statsTable = as.numeric(circle.statsTable),
			 inf.nbStats = as.integer(inf.nbStats),
			 inf.statsTable = as.numeric(inf.statsTable), 
			 atRisk.trs = as.numeric(atRisk.trs),
			 atRisk.ntrs = as.integer(length(atRisk.trs)),
			 atRisk.statsTable = as.numeric(atRisk.statsTable),
			 atRisk.nbCoefs = as.integer(atRisk.nbCoefs),
			 xs = as.numeric(coords$X),
			 ys = as.numeric(coords$Y),
			 detectRate = as.numeric(detectRate) 
			 )

		out$infestedDens<-out$infestedDens/Nrep;

		# make matrix out of semivar.statsTable
		out$semivar.statsTable<-matrix(out$semivar.statsTable,byrow=FALSE,ncol=Nrep)

		# need to remove the ones that are NAN
		notNAN <- which(!is.nan(out$semivar.statsTable[, 1]))

		# which of the pairwise to keep in final statistics
		orderPairwise <- c("semivariance", "moran", "geary", "ripley")
		whichkeep <- match(whichPairwise, orderPairwise) - 1 #positions are 0 indexed
		
		if(any(is.na(whichkeep))){
			warning("some pairwise supplied that are NA")
	      		whichkeep <- whichkeep[!is.na[whichkeep]]
		}

		keepable <- unlist(lapply(length(cbin)*whichkeep, "+", 1:length(cbin)))

		out$semivar.statsTable <- out$semivar.statsTable[intersect(keepable, notNAN), ]
	
		# make matrix out of grid.statsTable
		out$grid.statsTable <- matrix(out$grid.statsTable,byrow=FALSE,ncol=Nrep)

		nGridUniqueStats <- 2 + max(nPartCoefs, 1) + grid.numLmoments
		if(nPartCoefs > 0)
			keepGrid <- 2 + c(1:nPartCoefs, nPartCoefs+iPartLMoments)
		else
			keepGrid <- 2 + 1 + iPartLMoments
		keepGrid = keepGrid %% nGridUniqueStats
		
		keepIndices <- which((1:dim(out$grid.statsTable)[1] %% nGridUniqueStats) %in% keepGrid)
		out$grid.statsTable <- out$grid.statsTable[keepIndices, ] 

		# make matrix out of circle.statsTable
		out$circle.statsTable <- matrix(out$circle.statsTable, byrow=FALSE, ncol=Nrep)

		# make matrix out of inf.statsTable
		out$inf.statsTable <- matrix(out$inf.statsTable, byrow = FALSE, ncol = Nrep)
		if(sum(maxInfest) == length(maxInfest)) # throw away macro/micro unit statistics
			out$inf.statsTable <- out$inf.statsTable[-(2:3), ]

		# make matrix out of atRisk.statsTable
		out$atRisk.statsTable <- matrix(out$atRisk.statsTable, byrow = FALSE, ncol = Nrep)

		statsTable <- 0
		degenerateStats <- integer(0)

		if(getStats){ ## if want to get statistics, need to make the statsTable
	
			# put all the stats into one list for making statsTable
			allStats <- list(out$semivar.statsTable, out$grid.statsTable, out$circle.statsTable, out$atRisk.statsTable, out$inf.statsTable)
			if(Nrep==1){
				## if only one repetition, stats have to be handled as vectors
				for(statsWant in matchStats)
			 		statsTable <- c(statsTable, allStats[[statsWant]])

				statsTable <- statsTable[-1]
			}else{

				statsTable <- matrix(0, 1, Nrep)
				for(statsWant in matchStats)
					statsTable <- rbind(statsTable, allStats[[statsWant]])
				
				statsTable <- statsTable[-1, ]
				
				#figure out which stats are degenerate (important to do prestats removal!)
				vars <- apply(statsTable, 1, var)
				degenerateStats <- which(vars == 0)
			}
		}

		infestH <- out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1
		out$indexInfest <- infestH
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]
		out$timeI <- timeH

		## add the compiled statsTable and degenerateStats to out
		oldnames <- names(out)
		length(out) <- length(out) + 2 
		names(out) <- c(oldnames, "statsTable", "degenerateStats")
		out$statsTable <- statsTable
		out$degenerateStats <- degenerateStats
			
		return(out)
	}

	ToCstratHopSkipJump <- function(stratHopSkipJump){
	  for(i in 1:length(stratHopSkipJump)){
	    stratHopSkipJump@colindices <- as.integer(stratHopSkipJump@colindices-1L)
	  }

	  return(stratHopSkipJump)
	}

	noKernelMultiGilStatTweak <- function(
	    stratHopSkipJump, 
		blockIndex, 
		infestH, 
		timeH, 
		endTime, 
		rateMove,
		rateHopInMove, 
		rateSkipInMove, 
		rateJumpInMove, 
		Nrep, 
		coords, 
		seed=1, 
		simul=TRUE,
		maxInfest=rep(1,dim(coords)[1]), 
		getStats=TRUE, 
		dist_out = NULL, 
		map.partitions = NULL, 
		conc.circs = NULL, 
		atRisk.trs = NULL,
		atRisk.ncoefs = NULL,
		typeStat = "semivariance",
	        whichPairwise = c("semivariance", "moran", "geary", "ripley"),	
		detectRate = 1, 
		rateIntro = 0){

		# do we have blocks?
		haveBlocks <- !is.null(blockIndex)

	  	L<-dim(coords)[1]
		indexInfest <- rep(-1, sum(maxInfest)) 
		timeI <- rep(-1, sum(maxInfest))
		infested <- rep(0, L)
		infestedDens<-rep(0,L)
		endIndex <- as.integer(length(infestH) - 1)
		cat("tot maxInfest:",sum(maxInfest),"\n");
		if(length(infestH)>0){
		  indexInfest[1:length(infestH)] <- infestH - 1
		  timeI[1:length(timeH)] <- timeH
		  infested <- MakeNinfestFromIndex(infestH,L)
		}else{
		  infestH <- -1 # avoid to pass a NULL to .C
		}

		# if stratHopSkipJump, throw error 	
		if(is.null(stratHopSkipJump)){
			stop("need to pass a stratified hop/skip/jump matrix; see generate_stratifed_mat")
		}

		# implemented stats
		implStats <- c("semivariance", "grid", "circles", "atRisk", "num_inf")

		# initialize all the statistics to 0
		matchStats <- 0

		# if getStats and specific statistics are used, then change their value
		# semivariance statistics
		nbins <- 0
		dist_indices <- 0
		cbin <- 0
		cbinas <- 0
		cbinsb <- 0
		semivar.nbStats <- 0
		semivar.statsTable <- 0

		# grid/partition statistics
		numDiffGrids <- 0
		gridIndexes <- 0
		gridNumCells <- 0
		gridCountCells <- 0
		grid.nbStats <- 0
		grid.numCoeffs <- 0 
		grid.numLmoments <- 0
		grid.statsTable <- 0

		#circle statistics
		numDiffCircles <- 0
		numDiffCenters <- 0
		circleIndexes <- 0
		circleCounts <- 0
		circle.nbStats <- 0
		circle.statsTable <- 0

		# num_inf statistics
		inf.nbStats <- 0
		inf.statsTable <- 0

		# at_risk statistics
		atRisk.nbCoefs<-0
		atRisk.statsTable<-0

		if(getStats){

			## pass the corresponding matched numbers to C instead of the name of the statistics
			matchStats <- match(typeStat, implStats)

			if(any(is.na(matchStats))){ #throw an error regarding stats not yet implemented
				stop(paste0(typeStat[is.na(matchStats)], " not implemented! Only implemented ", implStats))
			}

			if("num_inf" %in% typeStat){
				# 1 if no blocks, 3 if blocks
				# stats selection
				##===============
				# Number Locations Infested (location or macrounits > 0)
				# Number Units infested (total microunits +)
				# Number Units Infested / Number Locations Infested
				# If haveBlocks:
				#	Number Blocks Infested
				#	Number Infested / Number Blocks Infested
				inf.nbStats <- 3 + haveBlocks*2
				inf.statsTable <- mat.or.vec(inf.nbStats, Nrep)
			}	

			# if want to calculate semivariance stats
			if("semivariance" %in% typeStat){
				if(is.null(dist_out)){

					stop("no semivariance dist_out object passed, cannot execute!")
				}

				dist_indices <- dist_out$CClassIndex
				cbin <- dist_out$classSize
				nbins <- length(cbin) + 1 #number of breaks not bins

				# stats selection
				###===================================
				## CURRENT STATS:
				## General Semivariance (new - new)
				## Moran's I
				## Geary's C
				## Ripley's L
				##	= 4 * length(cbin)
				## if haveBlocks
				##	By block Semivariance (same block - across streets)
				##	By block Semivariance Std. Dev
				## 	= 6 * length(cbin)
				## DEPRECATED STATS:
				## General Semivariance Std. Dev. (new - new)
				## General Semivariance (old - new)
				## General Semivariance Std. Dev (old - new)
				###===================================
				if(!haveBlocks)
					semivar.nbStats <- 4*length(cbin)
				else{
					semivar.nbStats <- 6*length(cbin)
					cbinas <- dist_out$classSizeAS
					cbinsb <- dist_out$classSizeSB
				}

				semivar.statsTable<-mat.or.vec(semivar.nbStats,Nrep)
			}

			if("grid" %in% typeStat){

				if(is.null(map.partitions)){ # if an indexing of the map hasn't yet been passed, throw error
					stop("map.partitions passed as null, cannot execute!")
				}
				
				## the number of different indexings present in gridIndexes
				numDiffGrids <- length(map.partitions)

				## unlist the first list
				map.partitions <- unlist(map.partitions, recursive = FALSE)
		
				## define the indexing systems
				## every house should have numDiffGrids indexes
				gridIndexes <- map.partitions[which(names(map.partitions) %in% "index")]
				gridIndexes <- unlist(gridIndexes)
				names(gridIndexes) <- NULL

				## the length of each indexing - the number of cells per grid/index system
				gridNumCells <- map.partitions[which(names(map.partitions) %in% "num_cells")] 
				gridNumCells <- unlist(gridNumCells)
				names(gridNumCells) <- NULL

				# countCells will keep track of total number of houses per cell
				gridCountCells <- map.partitions[which(names(map.partitions) %in% "housesPerCell")]
				gridCountCells <- unlist(gridCountCells)
				names(gridCountCells) <- NULL

				# stats selection
				###===================================
				## CURRENT STATS 
				## (by grid system):
				## Variance of % positive per cell
				## Number Cells with at least 1 positive 
				## Fit quantile distribution to polynomial
				## a + bx + cx^2 + dx^3 + ... (grid.numCoeffs stats)
			        ##      = 2*numDiffGrids + 2*sum(grid.numCoeffs)
			       	## L-moment statistics (taken from quantile distribution)
			        ## (2nd, 3rd, 4th L-moments, L-scale, L-skewness, L-kurtosis)
				## L-mean should be ~ to median (also to num_inf)	
				###===================================
				grid.numCoeffs <- rep(1, length(gridNumCells)) #should normally be 4, 1 for quickness #grid.numCoeffs <- c(2, 4, 6, 6, 4, 2)
				grid.numLmoments <- 3 #should be 3
				grid.nbStats <- 2*numDiffGrids + sum(grid.numCoeffs) + grid.numLmoments*numDiffGrids	
				grid.statsTable <- mat.or.vec(grid.nbStats, Nrep)
			}

			if("circles" %in% typeStat){
				if(is.null(conc.circs)){ # if an indexing of concentric circles hasn't yet been passed, throw error
					stop("conc.circles passed as null, cannot execute!")
				}

				numDiffCenters <- dim(conc.circs$counts)[1]
				numDiffCircles <- dim(conc.circs$counts)[2]
				circleIndexes <- t(conc.circs$circleIndexes) 
				circleCounts <- t(conc.circs$counts) 

				# stats selection
				###===================================
				## CURRENT STATS 
				## (by numDiffCircles):
				## Variance of % positive (across initInfested)
				## Mean of % positive (across initInfested) 
				##	= 2 * numDiffCircles
				###===================================
				circle.nbStats <- 2*numDiffCircles
				circle.statsTable <- mat.or.vec(circle.nbStats, Nrep)
			}
			
			if("atRisk" %in% typeStat){
				atRisk.nbCoefs<-atRisk.ncoefs
				atRisk.nbStats<-length(atRisk.trs)+atRisk.ncoefs
				atRisk.statsTable<-mat.or.vec(Nrep,atRisk.nbStats)
			}
		}

		# if don't have blocks, pass 0 as the value for skipColIndex + skipRowPointer (note, these values should never be used since rateSkipInMove == 0)
		# also pass 0 as the value to blockIndex

		skipColIndex <- if(!is.null(stratHopSkipJump$skipMat)){
		  as.integer(stratHopSkipJump$skipMat@colindices-1L)
		}else{
		  as.integer(-1L)
		}
		skipRowPointer <-if(!is.null(stratHopSkipJump$skipMat)){
		  as.integer(stratHopSkipJump$skipMat@rowpointers-1L)}else{
		    as.integer(-1L)}

		out<- .C("noKernelMultiGilStat",
			 # simulation parameters
			 hopColIndex = as.integer(stratHopSkipJump$hopMat@colindices-1L),
			 hopRowPointer = as.integer(stratHopSkipJump$hopMat@rowpointers-1L), 
			 skipColIndex = skipColIndex, # if no skips, pass dummy skip
			 skipRowPointer = skipRowPointer, # if no skips, pass dummy skip
			 jumpColIndex = as.integer(stratHopSkipJump$jumpMat@colindices-1L),
			 jumpRowPointer = as.integer(stratHopSkipJump$jumpMat@rowpointers-1L),
			 rateHopInMove = as.numeric(rateHopInMove), 
			 rateSkipInMove = as.numeric(rateSkipInMove), 
			 rateJumpInMove = as.numeric(rateJumpInMove), 
			 blockIndex = if(haveBlocks){as.integer(blockIndex)}else{as.integer(0)},
			 simul = as.integer(simul),
			 infested = as.integer(infested),
			 maxInfest = as.integer(maxInfest),
			 infestedDens = as.numeric(infestedDens),
			 endIndex = endIndex,
			 L = as.integer(L),
			 endTime = as.numeric(endTime),
			 indexInfest = as.integer(indexInfest),
			 timeI = as.numeric(timeI),
			 rateMove = as.numeric(rateMove),
			 rateIntro = as.numeric(rateIntro),
			 seed = as.integer(seed),
			 Nrep = as.integer(Nrep),
			 # stats
			 getStats=as.integer(getStats),
			 matchStats=as.integer(matchStats),
			 lengthStats=as.integer(length(matchStats)),
			 nbins = as.integer(nbins),
			 cbin = as.integer(cbin),
			 cbinas = as.integer(cbinas),
			 cbinsb = as.integer(cbinsb),
			 indices = as.integer(dist_indices),
			 semivar.statsTable = as.numeric(semivar.statsTable),
			 semivar.nbStats = as.integer(semivar.nbStats),
			 haveBlocks = as.integer(haveBlocks),
			 numDiffGrids = as.integer(numDiffGrids),
			 gridIndexes = as.integer(gridIndexes-1), #subtract 1 because C is 0 indexed
			 gridNumCells = as.integer(gridNumCells),
			 gridCountCells = as.integer(gridCountCells),
			 grid.nbStats = as.integer(grid.nbStats),
			 grid.numCoeffs = as.integer(grid.numCoeffs),
			 grid.numLmoments = as.integer(grid.numLmoments),
			 grid.statsTable = as.numeric(grid.statsTable),
			 numDiffCircles = as.integer(numDiffCircles),
			 numDiffCenters = as.integer(numDiffCenters),
			 circleIndexes = as.integer(circleIndexes), #already C indexed (was computed in C)
			 circleCounts = as.integer(circleCounts),
			 circle.nbStats = as.integer(circle.nbStats),
			 circle.statsTable = as.numeric(circle.statsTable),
			 inf.nbStats = as.integer(inf.nbStats),
			 inf.statsTable = as.numeric(inf.statsTable), 
			 atRisk.trs = as.numeric(atRisk.trs),
			 atRisk.ntrs = as.integer(length(atRisk.trs)),
			 atRisk.statsTable = as.numeric(atRisk.statsTable),
			 atRisk.nbCoefs = as.integer(atRisk.nbCoefs),
			 xs = as.numeric(coords$X),
			 ys = as.numeric(coords$Y),
			 detectRate = as.numeric(detectRate) 
			 )

		out$infestedDens<-out$infestedDens/Nrep;

		# make matrix out of semivar.statsTable
		out$semivar.statsTable<-matrix(out$semivar.statsTable,byrow=FALSE,ncol=Nrep)

		# need to remove the ones that are NAN
		notNAN <- which(!is.nan(out$semivar.statsTable[, 1]))

		# which of the pairwise to keep in final statistics
		orderPairwise <- c("semivariance", "moran", "geary", "ripley")
		whichkeep <- match(whichPairwise, orderPairwise) - 1 #positions are 0 indexed
		
		if(any(is.na(whichkeep))){
			warning("some pairwise supplied that are NA")
	      		whichkeep <- whichkeep[!is.na[whichkeep]]
		}

		keepable <- unlist(lapply(length(cbin)*whichkeep, "+", 1:length(cbin)))

		out$semivar.statsTable <- out$semivar.statsTable[intersect(keepable, notNAN), ]
	
		# make matrix out of grid.statsTable
		out$grid.statsTable <- matrix(out$grid.statsTable,byrow=FALSE,ncol=Nrep)

		## only keep L-moments (if only want to keep 3rd, 4th, change to 5, 0) 
		keepLmoments <- which((1:dim(out$grid.statsTable)[1] %% 6) %in% c(4, 5, 0))
		out$grid.statsTable <- out$grid.statsTable[keepLmoments, ] 

		## only keep regression coefficients
		## keepCoeffs <- which((1:dim(out$grid.statsTable)[1] %% 7) %in% c(3, 4, 5, 6))
		## out$grid.statsTable <- out$grid.statsTable[keepCoeffs, ]

		# make matrix out of circle.statsTable
		out$circle.statsTable <- matrix(out$circle.statsTable, byrow=FALSE, ncol=Nrep)

		# make matrix out of inf.statsTable
		out$inf.statsTable <- matrix(out$inf.statsTable, byrow = FALSE, ncol = Nrep)

		# make matrix out of atRisk.statsTable
		out$atRisk.statsTable <- matrix(out$atRisk.statsTable, byrow = FALSE, ncol = Nrep)

		statsTable <- 0
		degenerateStats <- integer(0)

		if(getStats){ ## if want to get statistics, need to make the statsTable
	
			# put all the stats into one list for making statsTable
			allStats <- list(out$semivar.statsTable, out$grid.statsTable, out$circle.statsTable, out$atRisk.statsTable, out$inf.statsTable)
			if(Nrep==1){
				## if only one repetition, stats have to be handled as vectors
				for(statsWant in matchStats)
			 		statsTable <- c(statsTable, allStats[[statsWant]])

				statsTable <- statsTable[-1]
			}else{

				statsTable <- matrix(0, 1, Nrep)
				for(statsWant in matchStats)
					statsTable <- rbind(statsTable, allStats[[statsWant]])
				
				statsTable <- statsTable[-1, ]
				
				#figure out which stats are degenerate (important to do prestats removal!)
				vars <- apply(statsTable, 1, var)
				degenerateStats <- which(vars == 0)
			}
		}

		infestH <- out$indexInfest
		infestH <- infestH[which(infestH != -1)] + 1
		out$indexInfest <- infestH
		timeH <- out$timeI
		timeH <- timeH[which(infestH != -1)]
		out$timeI <- timeH

		## add the compiled statsTable and degenerateStats to out
		oldnames <- names(out)
		length(out) <- length(out) + 2 
		names(out) <- c(oldnames, "statsTable", "degenerateStats")
		out$statsTable <- statsTable
		out$degenerateStats <- degenerateStats
			
		return(out)
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
	
	# group the nodes that are within tr of each others
	# Nodes: 
	percolation_circle<-function(dists,tr){
	  n<-dim(dists)[1]
	  Groups<-rep(-1,n)
	  out<-.C("percolation_circle",
		  Nodes=as.integer(Groups),
		  n=as.integer(n),
		  dists=as.numeric(dists),
		  tr=as.numeric(tr))
	  return(out$Nodes)
	}
	
	# flag at risk nodes according to a vector of threshold distances
	# and a matrix of distances
	get_at_risk_indicator<-function(posnodes,dists,trs){
	  # posnodes: nodes generating the risk
	  # dists: square distance matrix for all the nodes
	  # trs: increasing threshold distances
	  n <- dim(dists)[1]
	  at_risk<-rep(0,length(trs)*n)
	  posnodes<- posnodes-1 # to correspond to C numering
	
	  out<-.C("get_at_risk_indicator",
		  at_risk = as.integer(at_risk),
		  n = as.integer(n),
		  posnodes = as.integer(posnodes),
		  nPosnodes = as.integer(length(posnodes)),
		  dists = as.numeric(dists),
		  trs = as.numeric(trs),
		  nTr = as.integer(length(trs))
		  )
	  return(matrix(out$at_risk,ncol=length(trs),byrow=TRUE))
	}
	
	# get indicator at_Risk and transform it in raw stat
	get_at_risk_stat<-function(pos,dists,trs){
		stat<- rep(0,length(trs))
		L<-dim(dists)[1]
		pos<- pos-1 # to correspond to C numering
	
		out<-.C("get_at_risk_stat",
			at_risk_current = as.numeric(stat),
			L = as.integer(L),
			pos = as.integer(pos),
			npos = as.integer(length(pos)),
			dists = as.numeric(dists),
			trs_at_risk = as.numeric(trs),
			ntr_at_risk = as.integer(length(trs))
			)
		rawstat<-out$at_risk_current
		return(rawstat)
	}
	
	# fit the raw stat and formalize in big matrix
	get_stats_at_risk<-function(numRep,pos,dists,trs,currentAtRiskStat,ncoefs){
		pos<- pos-1 # to correspond to C numering
		out<-.C("get_stats_at_risk",
			rep = as.integer(numRep),
			L = as.integer(dim(dists)[1]),
			pos = as.integer(pos),
			npos = as.integer(length(pos)),
			dists = as.numeric(dists),
			trs_at_risk = as.numeric(trs),
			ntr_at_risk = as.integer(length(trs)),
			at_risk_stat = as.double(currentAtRiskStat),
			ncoefsAtRisk = as.integer(ncoefs)
			)
	
		at_risk_stats<-matrix(out$at_risk_stat,ncol=dim(currentAtRiskStat)[2],byrow=TRUE)
		return(at_risk_stats)
	}

	get_stats_grid <- function(infested, maxInfest, map.partitions, nPartCoefs=0, iPartLMoments=c(1, 2, 3)){
				
		## the number of different indexings present in gridIndexes
		numDiffGrids <- length(map.partitions)

		## unlist the first list
		map.partitions <- unlist(map.partitions, recursive = FALSE)
	
		## define the indexing systems
		## every house should have numDiffGrids indexes
		gridIndexes <- map.partitions[which(names(map.partitions) %in% "index")]
		gridIndexes <- unlist(gridIndexes)
		names(gridIndexes) <- NULL

		## the length of each indexing - the number of cells per grid/index system
		gridNumCells <- map.partitions[which(names(map.partitions) %in% "num_cells")] 
		gridNumCells <- unlist(gridNumCells)
		names(gridNumCells) <- NULL

		# countCells will keep track of total number of houses per cell
		gridCountCells <- map.partitions[which(names(map.partitions) %in% "housesPerCell")]
		gridCountCells <- unlist(gridCountCells)
		names(gridCountCells) <- NULL

		# stats selection
		###===================================
		## CURRENT STATS 
		## (by grid system):
		## Variance of % positive per cell
		## Number Cells with at least 1 positive 
		## Fit quantile distribution to polynomial
		## a + bx + cx^2 + dx^3 + ... (grid.numCoeffs stats)
	        ##      = 2*numDiffGrids + 2*sum(grid.numCoeffs)
	       	## L-moment statistics (taken from quantile distribution)
	        ## (2nd, 3rd, 4th L-moments, L-scale, L-skewness, L-kurtosis)
		## L-mean should be ~ to median (also to num_inf)	
		###===================================
		grid.numCoeffs <- rep(max(1, nPartCoefs), length(gridNumCells)) #take max of nPartCoefs, 1 - need to pass in dummy object
		grid.numLmoments <- max(1, max(iPartLMoments)) #take the max of the indices 
		grid.nbStats <- 2*numDiffGrids + sum(grid.numCoeffs) + grid.numLmoments*numDiffGrids	
		grid.statsTable <- mat.or.vec(grid.nbStats, 1)
		out <- .C(	"get_stats_grid",
				rep = as.integer(0),
				L = as.integer(length(infested)),
				infestedInit = as.integer(infested),
				maxInfest = as.integer(maxInfest),
				gridnbStats = as.integer(grid.nbStats),
				numDiffGrids = as.integer(numDiffGrids),
				gridIndexes = as.integer(gridIndexes-1),
				gridNumCells = as.integer(gridNumCells),
				gridCountCells = as.integer(gridCountCells),
				numCoeffs = as.integer(grid.numCoeffs),
				numLmoments = as.integer(grid.numLmoments),
				gridstats = as.numeric(grid.statsTable))
		
		nGridUniqueStats <- 2 + max(nPartCoefs, 1) + grid.numLmoments
		if(nPartCoefs > 0){
			keepGrid <- 2 + c(1:nPartCoefs, nPartCoefs+iPartLMoments)
		}else{
			keepGrid <- 2 + 1 + iPartLMoments
		}
		keepGrid = keepGrid %% nGridUniqueStats
		
		keepIndices <- which((1:length(out$gridstats) %% nGridUniqueStats) %in% keepGrid)
		out$gridstats_mod <- out$gridstats[keepIndices] 


		return(list(gridstats = out$gridstats, gridstats_mod = out$gridstats_mod))
		
	}



}

#=====================
## Analysis of outputs of multiGils
## Way to simulate messiness in data
#=====================

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
  cat("; Opened Inf:", sum(observed));

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
