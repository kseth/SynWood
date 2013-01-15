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
###	     infestH=infestH,
###	     timeH=timeH,
###	     endTime=nbit,
###	     maps=maps,
###	     nbit=nbit,
###	     Nrep=Nrep,
###	     priorMeans=priorMeans,
###	     priorSdlog=priorSdlog,
###	     genIntervals=genIntervals,
###	     paramInf=paramInf,
###	     paramSup=paramSup,
###	     PGF=PGF,
###	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

noKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parm.names
	
	# coerce theta, a priori all positive
	# set intervals for the variables
	# take care of this problem by doing try catch on wood syn likelihood
	# for(param in names(theta)){
        # theta[param]<-interval(theta[param],a=Data$paramInf[param],b=Data$paramSup[param])
    	# }
	# cat("theta:",theta, "\n")

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
			    dist_out = Data$dist_out)

	end <- Sys.time()
	# cat("t multiGil:",end-start,"\n")

	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
	}else{
		yhat<-out$statsTable[,1]
       	
		# print(out$statsTable)	
		# # see for the following multiGilStat
		# statsToKeep<-1:length(yhat) # c(1:16,39:42)
		# statsTable<-out$statsTable[statsToKeep,]
		# y<-Data$y[statsToKeep]
		# statsTable<-out$statsTable[!duplicated(yhat),]
	
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
		LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	# return
	
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
###	     infestH=infestH,
###	     timeH=timeH,
###	     endTime=nbit,
###	     maps=maps,
###	     nbit=nbit,
###	     Nrep=Nrep,
###	     priorMeans=priorMeans,
###	     priorSdlog=priorSdlog,
###	     genIntervals=genIntervals,
###	     paramInf=paramInf,
###	     paramSup=paramSup,
###	     PGF=PGF,
###	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

binomNoKernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parm.names
	
	# coerce theta, a priori all positive
	# set intervals for the variables
	# take care of this problem by doing try catch on wood syn likelihood
	# for(param in names(theta)){
        # theta[param]<-interval(theta[param],a=Data$paramInf[param],b=Data$paramSup[param])
    	# }
	# cat("theta:",theta, "\n")

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
			    dist_out = Data$dist_out)

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
		# then multiply by (1-infestedDense) for all those that are not infested at end time point y
		ll <- prod(yhat[which(Data$y == 1)])
		ll <- ll * prod((1-yhat[which(Data$y == 0)]))		

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL
		LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	# return
	
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
###	     priorSdlog=priorSdlog,
###	     genIntervals=genIntervals,
###	     paramInf=paramInf,
###	     paramSup=paramSup,
###	     PGF=PGF,
###	     mon.names=c("LL","LP", names(priorMeans)), # monitored variables (like in Model)
###	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
###	     sampling=sampling # method of sampling parameters
###		)

kernelModel <- function(theta,Data,postDraw=FALSE){

	theta<-theta
	names(theta)<-Data$parm.names
	
	# coerce theta, a priori all positive
	# set intervals for the variables
	# take care of this problem by doing try catch on wood syn likelihood
	# for(param in names(theta)){
        # theta[param]<-interval(theta[param],a=Data$paramInf[param],b=Data$paramSup[param])
    	# }
	# cat("theta:",theta, "\n")

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
       	
		# print(out$statsTable)	
		# # see for the following multiGilStat
		# statsToKeep<-1:length(yhat) # c(1:16,39:42)
		# statsTable<-out$statsTable[statsToKeep,]
		# y<-Data$y[statsToKeep]
		# statsTable<-out$statsTable[!duplicated(yhat),]
	
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
		LP <- LL + sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	# return
	
	Modelout <- list(LP=LP, # joint posterior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP, theta), # to be monitored/ploted, carefull to change namesMonitor if modified
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}

 

 

