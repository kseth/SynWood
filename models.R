source("functions_migration.R")

#==================================
## Declare a global variable to hold minimum Likelihood
## 	use min LL in bad cases
#==================================
### "Pirat" data (global because need to be updated)
### If run into a case where Wood LD method doesn't work, use this LL
### Wood LD only doesn't work in extreme + unlikely cases
minLLever=-10e6

#==================================
## declare noKernel (sampling over pre-normalized weights) Model for LD
#==================================
# List of data to pass to model + sampler
#==================================
### Data <- list(y=statsData,
###	     trans=NULL,
###	     stratHopSkipJump = stratHopSkipJump,
###	     blockIndex=blockIndex,
###	     dist_out = makeDistClassesWithStreets(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals, blockIndex),
###	     map.partitions = map.partitions,
###	     conc.circs = NULL,
###	     useStats = useStats, 
###	     infestH=infestH,
###	     timeH=timeH,
###	     endTime=nbit,
###	     maps=maps,
###	     nbit=nbit,
###	     Nrep=Nrep,
###	     priorMeans=priorMeans,
###	     priorSd=priorSd,
###	     priorType=priorType,
###	     priorIntervals=priorIntervals,
###	     defaultDR=defaultDR,
###	     genIntervals=genIntervals,
###	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

noKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parm.names
	
	# set seed for c simulations
	seed <-runif(1,min=0,(2^16-1)) 
	
	# simulations
	start <- Sys.time()
	
	if(postDraw){
		getStats<-FALSE
		# then in wrapper use infested if Nrep =1 and infested$Dens if Nrep=a few
		# to show a posterior predictive map
	}else{
		getStats<-TRUE
	}
	
	#pass all variables correctly!
	out <- noKernelMultiGilStat(stratHopSkipJump = Data$stratHopSkipJump, blockIndex = Data$blockIndex, 
			    infestH = Data$infestH, timeH = Data$timeH, endTime = Data$nbit, 
			    rateMove = theta["rateMove"],
			    weightHopInMove = 1,
			    weightSkipInMove = theta["weightSkipInMove"],
			    weightJumpInMove = theta["weightJumpInMove"], 
			    Nrep = Data$Nrep, 
			    coords = Data$maps[, c("X", "Y")], 
			    breaksGenVar = Data$genIntervals,
			    seed=seed,
			    getStats=getStats,
			    dist_out=Data$dist_out, map.partitions=Data$map.partitions, conc.circs=Data$conc.circs, typeStat=Data$useStats,
			    detectRate=ifelse("detectRate" %in% Data$parm.names, theta["detectRate"], defaultDR))

	end <- Sys.time()

	## cat("t multiGil: ")
	## print(end-start)

	start <- Sys.time()
	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
	}else{
		yhat<-out$statsTable[,1]
		# synthetic likelihood
		epsilon <- 1/(Data$Nrep+1)
		degen <- out$degenerateStats # the degenerate stats all have dirac distributions
		print(degen)	
		if(length(degen) > 0 && length(degen) < length(Data$y)){ #if some but not all stats degenerate

			degenY <- Data$y[degen]
			degenTheta <- out$statsTable[degen, 1]

			# still try() b/c could have cov issues (cov==0 if stats are proportional)
			success<-try(ll<-synLik(out$statsTable[-degen,],Data$y[-degen],Data$trans))

			if(class(success)=="try-error")
				  ll<-minLLever
			else{
				  minLLever<<-min(ll,minLLever)	
				  numGood <- length(which(degenY - degenTheta == 0))

				  #put in epsilon for degeneracy
				  ll <- ll + numGood*log(1-epsilon) + (length(degen)-numGood)*log(epsilon) 
			}		
		}else if(length(degen) == length(Data$y)){ #if all stats degenerate

			degenY <- Data$y[degen]
			degenTheta <- out$statsTable[degen, 1]
			numGood <- length(which(degenY - degenTheta == 0))

			ll <- minLLever + numGood*log(1-epsilon) + (length(degen)-numGood)*log(epsilon)

		}else{
			# still try() b/c could have cov issues (cov==0 if stats are proportional)
			success<-try(ll<-synLik(out$statsTable,Data$y,Data$trans))

			if(class(success)=="try-error")
				  ll<-minLLever
			else
				  minLLever<<-min(ll,minLLever)	
		}

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL

		LP <- LL

		for(name in Data$parm.names){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
			if(Data$priorType[name] == "lnorm")
				priorLL <- dlnorm(theta[name], meanlog = log(Data$priorMeans[name]), sdlog = Data$priorSd[name], log = TRUE)
			else if(Data$priorType[name] == "norm")
				priorLL <- dnorm(theta[name], mean = Data$priorMeans[name], sd = Data$priorSd[name], log = TRUE)
			else if(Data$priorType[name] == "boundednorm")
				priorLL <- dtnorm(theta[name], mean = Data$priorMeans[name], sd = Data$priorSd[name], lower = Data$priorIntervals[[name]][1], upper = Data$priorIntervals[[name]][2], log = TRUE)
			else if(Data$priorType[name] == "boundedlnorm")
				priorLL <- {if(theta[name] > Data$priorIntervals[[name]][1] && theta[name] < Data$priorIntervals[[name]][2]) dlnorm(theta[name], meanlog = log(Data$priorMeans[name]), sdlog = Data$priorSd[name], log = TRUE) else -Inf}
			else priorLL <- 0

			attributes(priorLL)<-NULL
			LP <- LP + priorLL
		}

		# LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	end <- Sys.time()
	## cat("t synLik: ")
	## print(end-start)
	
	Modelout <- list(LP=LP, # joint posterior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP, theta), # to be monitored/ploted, carefull to change namesMonitor if modified
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}

#==================================
## declare noKernel Binomial Likelihood (sampling over pre-normalized weights) Model for LD
#==================================
# List of data to pass to model + sampler
#==================================
### Data <- list(y=binomEndInfested,
###	     trans=NULL,
###	     stratHopSkipJump = stratHopSkipJump,
###	     blockIndex=blockIndex,
###	     dist_out = NULL, 
###	     map.partitions = NULL,
###	     conc.circs = NULL,
###	     useStats = useStats,
###	     infestH=infestH,
###	     timeH=timeH,
###	     endTime=nbit,
###	     maps=maps,
###	     nbit=nbit,
###	     Nrep=Nrep,
###	     priorMeans=priorMeans,
###	     priorSd=priorSd,
###	     priorType=priorType,
###	     priorIntervals=priorIntervals,
###	     defaultDR=defaultDR,
###	     genIntervals=genIntervals,
###	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

binomNoKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parm.names

	# set seed for c simulations
	seed <-runif(1,min=0,(2^16-1)) 
	
	# simulations
	start <- Sys.time()
	
	# we don't need the summary statistics
	getStats <- FALSE
	
	#pass all variables correctly!
	out <- noKernelMultiGilStat(stratHopSkipJump = Data$stratHopSkipJump, blockIndex = Data$blockIndex, 
			    infestH = Data$infestH, timeH = Data$timeH, endTime = Data$nbit, 
			    rateMove = theta["rateMove"],
			    weightHopInMove = 1,
			    weightSkipInMove = theta["weightSkipInMove"],
			    weightJumpInMove = theta["weightJumpInMove"], 
			    Nrep = Data$Nrep, 
			    coords = Data$maps[, c("X", "Y")], 
			    breaksGenVar = Data$genIntervals,
			    seed=seed,
			    getStats=getStats,
			    dist_out = Data$dist_out, map.partitions = Data$map.partitions, conc.circs = Data$conc.circs, typeStat = Data$useStats,
			    detectRate=ifelse("detectRate" %in% Data$parm.names, theta["detectRate"], defaultDR))

	end <- Sys.time()
	# cat("t multiGil:",end-start,"\n")

	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
	}else{
		yhat<-out$infestedDens

		# binomial likelihood
		# take the product of all the infestedDens that are infested at the end time point y
		# then multiply by (1-infestedDens) for all those that are not infested at end time point y
		inf <- which(Data$y == 1)
		uninf <- which(Data$y == 0)
		pred1 <- which(yhat == 1)
		pred0 <- which(yhat == 0)
		predMid <- which(yhat > 0 & yhat < 1)

		# ff == fudge factor, the epsilon that takes into account the error from a limited number of simulations
		# ff <- 1/(length(Data$Nrep)+1)
		ff <- 10e-6

		ll <- sum(log(yhat[intersect(inf, union(predMid, pred1))]))
		ll <- ll + sum(log(1-yhat[intersect(uninf, union(predMid, pred0))]))
		ll <- ll + sum(log(rep(ff, length(union(intersect(inf, pred0), intersect(uninf, pred1))))))

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL

		LP <- LL

		for(name in Data$parm.names){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
			if(Data$priorType[name] == "lnorm")
				priorLL <- dlnorm(theta[name], meanlog = log(Data$priorMeans[name]), sdlog = Data$priorSd[name], log = TRUE)
			else if(Data$priorType[name] == "norm")
				priorLL <- dnorm(theta[name], mean = Data$priorMeans[name], sd = Data$priorSd[name], log = TRUE)
			else if(Data$priorType[name] == "boundednorm")
				priorLL <- dtnorm(theta[name], mean = Data$priorMeans[name], sd = Data$priorSd[name], lower = Data$priorIntervals[[name]][1], upper = Data$priorIntervals[[name]][2], log = TRUE)
			else if(Data$priorType[name] == "boundedlnorm")
				priorLL <- {if(theta[name] > Data$priorIntervals[[name]][1] && theta[name] < Data$priorIntervals[[name]][2]) dlnorm(theta[name], meanlog = log(Data$priorMeans[name]), sdlog = Data$priorSd[name], log = TRUE) else -Inf}
			else priorLL <- 0

			attributes(priorLL)<-NULL
			LP <- LP + priorLL
		}

		# LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	Modelout <- list(LP=LP, # joint posterior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP, theta), # to be monitored/ploted, carefull to change namesMonitor if modified
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}



#==================================
## declare kernel Model for LD
## this model is outdated and has not been used since November 2012
## samples over possibly rateMove, rateJumpInMove, rateSkipInMove, halfDistH, halfDistJ
#==================================
## Data must contain:

### If cumulProbMat set to NULL multiGilStat will use theta to compute it
### Must calculate cumulProbMat in C each time
### Data <- list(y=statsData,
###	     trans=NULL,
###	     cumulProbMat=NULL,
###	     blockIndex=blockIndex,
###	     dist_out = makeDistClassesWithStreets(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals, blockIndex), 
###	     infestH=infestH,
###	     timeH=timeH,
###	     endTime=nbit,
###	     maps=maps,
###	     nbit=nbit,
###	     Nrep=Nrep,
###	     priorMeans=priorMeans,
###	     priorSd=priorSd,
###	     priorType=priorType,
###	     priorIntervals=priorIntervals,
###	     genIntervals=genIntervals,
###	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

kernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parm.names

	# set seed for c simulations
	seed <-runif(1,min=0,(2^16-1)) 
	
	# simulations
	start <- Sys.time()
	
	if(postDraw){
		getStats<-FALSE
		# then in wrapper use infested if Nrep =1 and infested$Dens if Nrep=a few
		# to show a posterior predictive map
	}else{
		getStats<-TRUE
	}
	
	# if only sampled parameter is rateMove
	if(length(theta) == 1 & all.equal(Data$parm.names[1], "rateMove")){
		out <- multiGilStat(cumulProbMat = Data$cumulProbMat, blockIndex = Data$blockIndex, 
			    infestH = Data$infestH, timeH = Data$timeH, endTime = Data$nbit, 
			    rateMove = theta["rateMove"], 
			    Nrep = Data$Nrep, 
			    coords = Data$maps[, c("X", "Y")], 
			    breaksGenVar = Data$genIntervals,
			    getStats=getStats,
			    seed=seed,
			    dist_out = Data$dist_out)
	}else{
		#pass all variables correctly!
		#need to make probMat everytime, this is more time consuming.
		#note rateHopInMove = 1-rateJumpInMove
		out <- multiGilStat(cumulProbMat = Data$cumulProbMat, blockIndex = Data$blockIndex, 
			    infestH = Data$infestH, timeH = Data$timeH, endTime = Data$nbit, 
			    rateMove = theta["rateMove"], 
			    Nrep = Data$Nrep, 
			    coords = Data$maps[, c("X", "Y")], 
			    breaksGenVar = Data$genIntervals,
			    getStats=getStats,
			    seed=seed,
			    halfDistJ = ifelse(is.na(theta["halfDistJ"]), halfDistJ, theta["halfDistJ"]), # default is halfDistJ if not in theta
			    halfDistH = ifelse(is.na(theta["halfDistH"]), halfDistH, theta["halfDistH"]), 
			    useDelta = ifelse(is.na(theta["useDelta"]), useDelta, theta["useDelta"]), 
			    delta = ifelse(is.na(theta["delta"]), delta, theta["delta"]), 
			    rateHopInMove = ifelse(is.na(theta["rateJumpInMove"]), 1-rateJumpInMove, 1-theta["rateJumpInMove"]), 
			    rateSkipInMove = ifelse(is.na(theta["rateSkipInMove"]), rateSkipInMove, theta["rateSkipInMove"]), 
			    rateJumpInMove = ifelse(is.na(theta["rateJumpInMove"]), rateJumpInMove, theta["rateJumpInMove"]),
			    dist_out = Data$dist_out)
	}

	end <- Sys.time()
	# cat("t multiGil:",end-start,"\n")

	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
	}else{
		yhat<-out$statsTable[,1]

		# synthetic likelihood
		success<-try(ll<-synLik(out$statsTable,Data$y,Data$trans))
		if(class(success)=="try-error"){
		  ll<-minLLever
		}else{
		  minLLever<<-min(ll,minLLever)
		}

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL

		LP <- LL
		for(name in Data$parm.names){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
			if(Data$priorType[name] == "lnorm")
				priorLL <- dlnorm(theta[name], meanlog = log(Data$priorMeans[name]), sdlog = Data$priorSd[name], log = TRUE)
			else if(Data$priorType[name] == "norm")
				priorLL <- dnorm(theta[name], mean = Data$priorMeans[name], sd = Data$priorSd[name], log = TRUE)
			else if(Data$priorType[name] == "boundednorm")
				priorLL <- dtnorm(theta[name], mean = Data$priorMeans[name], sd = Data$priorSd[name], lower = Data$priorIntervals[[name]][1], upper = Data$priorIntervals[[name]][2], log = TRUE)
			else if(Data$priorType[name] == "boundedlnorm")
				priorLL <- {if(theta[name] > Data$priorIntervals[[name]][1] && theta[name] < Data$priorIntervals[[name]][2]) dlnorm(theta[name], meanlog = log(Data$priorMeans[name]), sdlog = Data$priorSd[name], log = TRUE) else -Inf}
			else priorLL <- 0

			attributes(priorLL)<-NULL
			LP <- LP + priorLL
		}


		# LP<-LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}
	
	Modelout <- list(LP=LP, # joint posterior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP, theta), # to be monitored/ploted, carefull to change namesMonitor if modified
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}

