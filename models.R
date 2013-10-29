source("functions_migration.R")
library(sn)

#==================================
# Binomial likelihood function
#==================================
# yhat - the mean observation over all iterations (the infestation density of the field)
# y - the observed infestation
# ff - fudge factor, the epsilon that takes into account the error from a limited number of simulations

binLik <- function(yhat, y, ff=10e-6, log=TRUE){

	probs <- 1 - abs(yhat - y)
	probs[which(probs == 0)] <- ff

	if(log){
	        ll <- sum(log(probs))
	}else{
		ll <- prod(probs)
	}

	return(ll)
}

#==================================
## Declare a global variable to hold minimum Likelihood
## 	use min LL in bad cases
#==================================
### "Pirat" data (global because need to be updated)
### If run into a case where Wood LD method doesn't work, use this LL
### Wood LD only doesn't work in extreme + unlikely cases
minLLever=-10e6

#==================================
## declare skewNoKernelModel (sampling over pre-normalized weights) Model for LD
#==================================
# List of data to pass to model + sampler
#==================================
### Data <- list(y=statsData,
###	     trans=NULL,
###	     stratHopSkipJump = stratHopSkipJump,
###	     blockIndex=blockIndex,
###	     dist_out = dist_out, #semivariance object
###	     map.partitions = map.partitions, #partitions object
###	     conc.circs = conc.circs, #circles object
###	     useStats = useStats, 
###	     whichPairwise = whichPairwise,
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
###	     initValues=initValues,
###	     default=default,
###	     monNames=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parmNames=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

skewNoKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parmNames
	
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
			    rateMove = ifelse("rateMove" %in% Data$parmNames, theta["rateMove"], Data$default["rateMove"]),
			    weightSkipInMove = ifelse("weightSkipInMove" %in% Data$parmNames, theta["weightSkipInMove"], Data$default["weightSkipInMove"]),
			    weightJumpInMove = ifelse("weightJumpInMove" %in% Data$parmNames,theta["weightJumpInMove"], Data$default["weightJumpInMove"]),
			    Nrep = Data$Nrep, 
			    coords = Data$maps[, c("X", "Y")], 
			    seed=seed,
			    getStats=getStats,
			    dist_out=Data$dist_out, map.partitions=Data$map.partitions, conc.circs=Data$conc.circs, typeStat=Data$useStats,
			    detectRate=ifelse("detectRate" %in% Data$parmNames, theta["detectRate"], Data$default["detectRate"]),
			    rateIntro=ifelse("rateIntro" %in% Data$parmNames, theta["rateIntro"], Data$default["rateIntro"]), 
			    whichPairwise=Data$whichPairwise)

	end <- Sys.time()

	statsTable <- t(out$statsTable)

	## cat("t multiGil: ")
	## print(end-start)

	start <- Sys.time()
	if(postDraw){
		yhat<-statsTable[,1]
		LL<-NA
		LP<-NA
	}else{
		yhat<-statsTable[1,]
		# synthetic likelihood
		degen <- out$degenerateStats # the degenerate stats all have dirac distributions	
		if(length(degen) > 0){ #if some stats degenerate
			if(length(degen) == length(Data$y)) #if all stats are degenerate
				ll <- -Inf
			else{
				degenY <- Data$y[degen]
				degenTheta <- statsTable[1, degen]

				numGood <- length(which(degenY - degenTheta == 0)) # stats exactly correct
				
				if(numGood == length(degenY)){ #if all degenerate stats match their actual values
					# still try() b/c could have cov issues (cov==0 if stats are proportional)
					success<-try(skew.fit<-msn.mle(y=statsTable[, -degen]))

					if(class(success)=="try-error"){
						  ll<-minLLever
					}else{
						  ll<-dmsn(Data$y[-degen],dp=skew.fit$dp,log=TRUE)
						  minLLever<<-min(ll,minLLever)	
					}
				}else{ #if even one degenerate stat does not match actual value
						ll <- -Inf
				}
			}
		}else{
			# still try() b/c could have cov issues (cor==1 if stats are proportional)
			success<-try(skew.fit<-msn.mle(y=statsTable))

			if(class(success)=="try-error"){
				  ll<-minLLever
			}else{
				ll<-dmsn(Data$y,dp=skew.fit$dp,log=TRUE)
				minLLever<<-min(ll,minLLever)	
			}
		}

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL

		LP <- LL
		Lprioronly <- 0 #keep track of just the prior sum

		for(name in Data$parmNames){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
			if(is.function(Data$priorType[name]))
				priorLL <- Data$priorType[name](theta[name])
			else if(Data$priorType[name] == "lnorm")
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
			Lprioronly <- Lprioronly + priorLL
		}

		# LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	end <- Sys.time()
	## cat("t synLik: ")
	## print(end-start)
	
	Modelout <- list(LP=LP, # joint posterior
			 LL=LL, # likelihood of stats w/o prior
			 Lprioronly=Lprioronly, # likelihood of parmeters according only to the prior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP,theta), # to be monitored/ploted, carefull to change namesMonitor if modified
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}

#==================================
## declare noKernel (sampling over pre-normalized weights) Model for LD
#==================================
# List of data to pass to model + sampler
#==================================
### Data <- list(y=statsData,
###	     trans=NULL,
###	     stratHopSkipJump = stratHopSkipJump,
###	     blockIndex=blockIndex,
###	     dist_out = dist_out, #semivariance object
###	     map.partitions = map.partitions, #partitions object
###	     conc.circs = conc.circs, #circles object
###	     useStats = useStats, 
###	     whichPairwise = whichPairwise,
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
###	     initValues=initValues,
###	     default=default,
###	     monNames=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parmNames=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

noKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parmNames
	
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
	out <- noKernelMultiGilStat(
		stratHopSkipJump = Data$stratHopSkipJump, 
		blockIndex = Data$blockIndex, 
		infestH = Data$infestH, 
		timeH = Data$timeH, 
		endTime = Data$nbit, 
		rateMove = ifelse("rateMove" %in% Data$parmNames, theta["rateMove"], Data$default["rateMove"]),
		weightSkipInMove = ifelse("weightSkipInMove" %in% Data$parmNames, theta["weightSkipInMove"], Data$default["weightSkipInMove"]),
		weightJumpInMove = ifelse("weightJumpInMove" %in% Data$parmNames,theta["weightJumpInMove"], Data$default["weightJumpInMove"]),
		Nrep = Data$Nrep, 
		coords = Data$maps[, c("X", "Y")], 
		seed=seed,
		simul=TRUE,
		maxInfest=Data$maxInfest, 
		getStats=getStats,
		dist_out=Data$dist_out, 
		map.partitions=Data$map.partitions, 
		conc.circs=Data$conc.circs, 
		typeStat=Data$useStats,
		whichPairwise=Data$whichPairwise,
		detectRate=ifelse("detectRate" %in% Data$parmNames, theta["detectRate"], Data$default["detectRate"]),
		rateIntro=ifelse("rateIntro" %in% Data$parmNames, theta["rateIntro"], Data$default["rateIntro"]))

	end <- Sys.time()

	cat("t multiGil: ")
	print(end-start)

	start <- Sys.time()
	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
		Lprioronly <- 0 #keep track of just the prior sum
	}else{
		yhat<-out$statsTable[,1]
		# synthetic likelihood
		degen <- out$degenerateStats # the degenerate stats all have dirac distributions	
		if(length(degen) > 0){ #if some stats degenerate
			if(length(degen) == length(Data$y))
				ll <- -Inf
			else{
				degenY <- Data$y[degen]
				degenTheta <- out$statsTable[degen, 1]

				numGood <- length(which(degenY - degenTheta == 0)) # stats exactly correct
				
				if(numGood == length(degenY)){ #if all degenerate stats match their actual values
					# still try() b/c could have cov issues (cov==0 if stats are proportional)
					success<-try(ll<-synLik(out$statsTable[-degen,],Data$y[-degen],Data$trans))
					if(class(success)=="try-error")
						  ll<-minLLever
					else
						  minLLever<<-min(ll,minLLever)	
				}else{ #if even one degenerate stat does not match actual value
						ll <- -Inf
				}
			}
		}else{
			# still try() b/c could have cov issues (cor==1 if stats are proportional)
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
		Lprioronly <- 0 #keep track of just the prior sum

		for(name in Data$parmNames){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
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
			Lprioronly <- Lprioronly + priorLL
		}

		# LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	end <- Sys.time()
	cat("t synLik: ")
	print(end-start)
	
	Modelout <- list(LP=LP, # joint posterior
			 LL=LL, # likelihood of stats w/o prior
			 Lprioronly=Lprioronly, # likelihood of parmeters according only to the prior
			 Dev=-2*LL, # deviance, probably not to be changed
			 monitor=c(LL,LP,theta), # to be monitored/ploted, carefull to change namesMonitor if modified
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
###	     whichPairwise = whichPairwise,
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
###	     initValues=initValues,
###	     default=default,
###	     monNames=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parmNames=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

binomNoKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parmNames

	# set seed for c simulations
	seed <-runif(1,min=0,(2^16-1)) 
	
	# simulations
	start <- Sys.time()
	
	# we don't need the summary statistics
	getStats <- FALSE
	
	#pass all variables correctly!
	out <- noKernelMultiGilStat(
		stratHopSkipJump = Data$stratHopSkipJump, blockIndex = Data$blockIndex, 
		infestH = Data$infestH, timeH = Data$timeH, endTime = Data$nbit, 
		rateMove = ifelse("rateMove" %in% Data$parmNames, theta["rateMove"], Data$default["rateMove"]),
		weightSkipInMove = ifelse("weightSkipInMove" %in% Data$parmNames, theta["weightSkipInMove"], Data$default["weightSkipInMove"]),
		weightJumpInMove = ifelse("weightJumpInMove" %in% Data$parmNames,theta["weightJumpInMove"], Data$default["weightJumpInMove"]),
		Nrep = Data$Nrep, 
		coords = Data$maps[, c("X", "Y")], 
		seed=seed,
		simul=TRUE,
		maxInfest=Data$maxInfest,
		getStats=getStats,
		dist_out=Data$dist_out, 
		map.partitions=Data$map.partitions, 
		conc.circs=Data$conc.circs, 
		typeStat=Data$useStats,
		detectRate=ifelse("detectRate" %in% Data$parmNames, theta["detectRate"], Data$default["detectRate"]),
		rateIntro=ifelse("rateIntro" %in% Data$parmNames, theta["rateIntro"], Data$default["rateIntro"]))

	end <- Sys.time()
	# cat("t multiGil:",end-start,"\n")

	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
	}else{
		yhat<-out$infestedDens
		ll<-binLik(yhat=yhat, y=Data$y, ff=10e-6)

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL

		LP <- LL
		Lprioronly <- 0 #keep track of prior likelihood

		for(name in Data$parmNames){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
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
			Lprioronly <- Lprioronly + priorLL
		}

		# LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	Modelout <- list(LP=LP, # joint posterior
			 LL=LL, # likelihood of stats w/o prior
			 Lprioronly=Lprioronly, #likelihood of parameters according to prior
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
###	     initValues=initValues,
###	     genIntervals=genIntervals,
###	     monNames=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parmNames=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

kernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parmNames

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
	if(length(theta) == 1 & all.equal(Data$parmNames[1], "rateMove")){
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
		Lprioronly <- 0

		for(name in Data$parmNames){ #factor the priors in (if don't want to use priors, just pass priorType not listed)
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
			Lprioronly <- Lprioronly + priorLL
		}


		# LP<-LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}
	
	Modelout <- list(LP=LP, # joint posterior
			 LL=LL, # likelihood of stats w/o prior
			 Lprioronly=Lprioronly, #likelihood of parameters according to prior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP, theta), # to be monitored/ploted, carefull to change namesMonitor if modified
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}

