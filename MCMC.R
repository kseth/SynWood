source("functions_sampling.R")

#===========================
## The function containing the adaptive metropolis hastings
##	MyDataFullSample - has all the data for the Model
##		see models.R for specifications on what is expected in each MyDataFullSample (aka Data)	 
##	Model - the sampling model
##	functionSample - the sampling function (defaults to omniSample in functions_sampling.R)
##	nbsimul - the starting number of simulations (also the max number if useAutoStop = FALSE)
##	upFreq, saveFreq - how frequently to print updates to the stdout, save updates to save files 
##	sdprop - the initial standard deviations of proposals
##	adaptOK - need to adapt sdprop?
##	if(!adaptOK), sdprop should allow acceptance between low- and highAcceptRate
##	useAutoStop - use geweke + raftery to test for convergence + autostop?
##	checkAutoStop - initial #iterations after finishing adapting sampling variance when to check autostop
#===========================

MCMC <- function(MyDataFullSample, Model, functionSample = omniSample, nbsimul = 600, upFreq = 1, saveFreq = 20, sdprop = 0.4, adaptOK = FALSE, checkAdapt = 20, lowAcceptRate = 0.15, highAcceptRate = 0.40, useAutoStop = TRUE, checkAutoStop = 100, monitor.file = "thetasamples_all.txt"){

#===========================
# Init values 
#===========================

	theta <- MyDataFullSample$initValues
	nparams <- length(theta)
	
	if(length(sdprop) != nparams)
		sdprop <- rep(sdprop[1], nparams)

	names(sdprop) <- MyDataFullSample$parmNames

	beginEstimate <- -1 #value containing position of when adaptation of sampling variance is complete

	# init of theta attributes and saving scheme
	outModel<-Model(theta,MyDataFullSample)
	Monitor<-mat.or.vec(nbsimul+1,length(MyDataFullSample$monNames))
	Monitor[1,]<-outModel$Monitor
	attributes(theta)$outModel<-outModel

	accepts<-as.data.frame(matrix(rep(0,nparams*nbsimul),ncol=nparams))
	names(accepts)<-MyDataFullSample$parmNames

	# write headers to table
	write.table(t(MyDataFullSample$monNames), monitor.file, sep="\t",append=FALSE,col.names=FALSE,row.names=FALSE)
	write.table(t(MyDataFullSample$parmNames), paste("acceptsamples", monitor.file, sep = ""), sep="\t",append=FALSE,col.names=FALSE,row.names=FALSE)


#============================
# Launch sampler
#============================

	#Rprof()
	
	numit <- 1

	while(numit <= nbsimul){		

		## display at every update frequency
		if(upFreq!=0 && numit%%upFreq==0){
			  cat("it:",numit,"of", nbsimul, "current theta:", theta,"\n");
		}

		##save at every save frequency
		if(saveFreq!=0 && numit%%saveFreq==0){
			write.table(Monitor[(numit-saveFreq+1):numit, ], monitor.file, sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			write.table(accepts[(numit-saveFreq):(numit-1), ], paste("acceptsamples", monitor.file, sep = ""), sep ="\t",append=TRUE,col.names=FALSE,row.names=FALSE) ## very first save index can be zero, it's okay
			cat("numit: ",numit,"\nnbsimul: ",nbsimul,"\nadaptOK :",adaptOK,"\ncheckAdapt: ",checkAdapt,"\nsdprop: ", sdprop,"\nbeginEstimate: ",beginEstimate,"\nuseAutoStop: ",useAutoStop,"\ncheckAutoStop: ",checkAutoStop,"\nsaveFreq: ",saveFreq, "\n\n", file =  paste("MCMCvalues", monitor.file, sep = ""), append = TRUE)
}
	
		## adapt the sampling variance	
		if(!adaptOK && numit%%checkAdapt==0){
			adaptOK <- TRUE
			for(paramName in MyDataFullSample$parmNames){
	      			logSDprop <- adaptSDProp(sdprop[paramName], accepts[1:numit, paramName], lowAcceptRate, highAcceptRate, tailLength = 20)
				adaptOK <- adaptOK && attributes(logSDprop)$noupdate
				sdprop[paramName] <- logSDprop
	    	 	}
			
			if(adaptOK){ #if we've just adapted

				cat("adaption of sampling variance complete: beginning final run from numit: ", numit, "\n")
				beginEstimate <- numit
				nbsimul <- beginEstimate + nbsimul
				Monitor <- resized(Monitor, nr=nbsimul+1) #resize the Monitor
				accepts <- resized(accepts, nr=nbsimul) #resize the accepts
			}

			## if the variances haven't yet been adapted and running out of simulations
			## double the simulation size and resize
			if(!adaptOK && numit+checkAdapt >= nbsimul)
			{
				
				cat("sampling variance adaptation not complete: numit: ", numit, "doubling number simulations\n")
				nbsimul <- nbsimul * 2
				Monitor<-resized(Monitor,nr=nbsimul+1) #resize the Monitor
				accepts<-resized(accepts,nr=nbsimul) #resize the accepts

			}
		}
		
		## sample the variables
		for(paramName in MyDataFullSample$parmNames){
	      		theta<-functionSample(Model,MyDataFullSample,theta,paramName,sdprop[paramName])
	      		accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    	 }

	    	Monitor[numit+1,]<-attributes(theta)$outModel$Monitor

		## auto stopping
		if(useAutoStop && adaptOK && numit == beginEstimate + checkAutoStop){

			cb<-cb.diag(Monitor[(1+beginEstimate):numit, ],logfile="convergence_tests.txt")
			checkAutoStop<-min(cb$newNbIt,(numit-beginEstimate)*3)
			
			if(!cb$ok){
					cat("checking auto stop: numit: ", numit, "next check: ", checkAutoStop)
					
					if(nbsimul <= beginEstimate + checkAutoStop)
					{
						nbsimul <- beginEstimate + checkAutoStop
						Monitor<-resized(Monitor,nr=nbsimul+1) #resize the Monitor
						accepts<-resized(accepts,nr=nbsimul) #resize the accepts
					}
			}
			else{
				cat("auto stop okay! terminating chain\n")
				break
			}
		}	
	
		# increase the iteration count
		numit <- numit + 1
	}

	if(!exists("cb")){
		cat("didn't check autostop, no automatic burnIn, throwing away first 500 iterations\n")
		cb <- list(burnIn = 500)
	}

	write.table(Monitor[(beginEstimate+cb$burnIn):min(numit, nbsimul), ], paste0("complete", monitor.file), sep="\t",append=FALSE,col.names=TRUE,row.names=FALSE)
	write.table(accepts[(beginEstimate+cb$burnIn):min(numit, nbsimul), ], paste0("completeaccepts", monitor.file), sep ="\t",append=FALSE,col.names=TRUE,row.names=FALSE)
	
	# Rprof(NULL)
}
