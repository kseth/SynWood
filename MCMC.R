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

MCMC <- function(MyDataFullSample, Model, functionSample = omniSample, nbsimul = 600, upFreq = 1, saveFreq = 20, sdprop = rep(0.4, nparams), adaptOK = FALSE, checkAdapt = 20, lowAcceptRate = 0.15, highAcceptRate = 0.40, useAutoStop = TRUE, checkAutoStop = 100, monitor.file = "thetasamples_all.txt"){

#===========================
# Init values 
#===========================

	theta <- MyDataFullSample$priorMeans
	nparams = length(theta)
	names(sdprop) <- MyDataFullSample$parm.names

	beginEstimate <- -1 #value containing position of when adaptation of sampling variance is complete

	# init of theta attributes and saving scheme
	outModel<-Model(theta,MyDataFullSample)
	Monitor<-mat.or.vec(nbsimul+1,length(MyDataFullSample$mon.names))
	Monitor[1,]<-outModel$Monitor
	attributes(theta)$outModel<-outModel

	accepts<-as.data.frame(matrix(rep(0,nparams*nbsimul),ncol=nparams))
	names(accepts)<-MyDataFullSample$parm.names

	# write headers to table
	write.table(t(MyDataFullSample$mon.names), monitor.file, sep="\t",append=FALSE,col.names=FALSE,row.names=FALSE)
	write.table(t(MyDataFullSample$parm.names), paste("acceptsamples", monitor.file, sep = ""), sep="\t",append=FALSE,col.names=FALSE,row.names=FALSE)


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
			write.table(Monitor[(numit-saveFreq):numit, ], monitor.file, sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			write.table(accepts[(numit-saveFreq):numit, ], paste("acceptsamples", monitor.file, sep = ""), sep ="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			cat("numit: ",numit,"\nnbsimul: ",nbsimul,"\nadaptOK :",adaptOK,"\ncheckAdapt: ",checkAdapt,"\nsdprop: ", sdprop,"\nbeginEstimate: ",beginEstimate,"\nuseAutoStop: ",useAutoStop,"\ncheckAutoStop: ",checkAutoStop,"\nsaveFreq: ",saveFreq, "\n", file =  paste("MCMCvalues", monitor.file, sep = ""), append = TRUE)
}
	
		## adapt the sampling variance	
		if(!adaptOK && numit%%checkAdapt==0){
			adaptOK <- TRUE
			for(paramName in MyDataFullSample$parm.names){
	      			logSDprop <- adaptSDProp(sdprop[paramName], accepts[1:numit, paramName], lowAcceptRate, highAcceptRate, tailLength = 20)
				adaptOK <- adaptOK && attributes(logSDprop)$noupdate
				sdprop[paramName] <- logSDprop
	    	 	}
			
			if(adaptOK){

				cat("adaption of sampling variance complete: beginning final run from numit: ", numit, "\n")
				beginEstimate <- numit
				nbsimul <- beginEstimate + nbsimul
			}

			## if the variances haven't yet been adapted and running out of simulations
			## double the simulation size and resize
			if(!adaptOK && numit+checkAdapt >= nbsimul)
			{
				
				cat("sampling variance adaptation not complete: numit: ", numit, "doubling number simulations\n")
				nbsimul <- nbsimul * 2
				Monitor<-resized(Monitor,nr=nbsimul+1)
				accepts<-resized(accepts,nr=nbsimul)

			}
		}
		
		## sample the variables
		for(paramName in MyDataFullSample$parm.names){
	      		theta<-functionSample(Model,MyDataFullSample,theta,paramName,sdprop[paramName])
	      		accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    	 }

	    	Monitor[numit+1,]<-attributes(theta)$outModel$Monitor

		## auto stopping
		if(useAutoStop && adaptOK && numit == beginEstimate + checkAutoStop){

			cb<-cb.diag(Monitor[(1+beginEstimate):numit, ],logfile="convergence_tests.txt")
			checkAutoStop<-min(cb$newNbIt,numit*3)
			
			if(!cb$ok){
					cat("checking auto stop: numit: ", numit, "next check: ", numit + checkAutoStop)
					
					if(nbsimul <= beginEstimate + checkAutoStop)
					{
						nbsimul <- beginEstimate + checkAutoStop
						Monitor<-resized(Monitor,nr=nbsimul+1)
						accepts<-resized(accepts,nr=nbsimul)
					}
			}
			else{
				finalTestResults <- cb
				cat("auto stop okay! terminating chain")
				break
			}
		}	
	
		# increase the iteration count
		numit <- numit + 1
	}

	write.table(Monitor[1:min(numit, nbsimul), ], "thetasamplescomplete.txt", sep="\t",append=FALSE,col.names=TRUE,row.names=FALSE)
	write.table(accepts[1:min(numit, nbsimul), ], "acceptsamplescomplete.txt", sep ="\t",append=FALSE,col.names=TRUE,row.names=FALSE)
	
	# Rprof(NULL)

	stop("MCMC finished")
}
