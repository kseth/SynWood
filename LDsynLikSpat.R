library(LaplacesDemon)
library(testthat)
library(verification)
source("RanalysisFunctions.R")
source("logLik.r") # override some sl functions and add synLik
source("functions_sampling.R")
source("convergence_test.r")
source("param.r")

## load data 
## use simulated data as in base_regression.R
source("maps_basic_regression.R") # gives maps
source("functions_migration.R")

#==================
#  Params, if modifying stats: look for multiGilStat
#==================
Nrep=100
set.seed(1)
genIntervals <- c(seq(10, 100, 15), seq(130, 250, 30)) # distance classes for the general variogram

#===================
# Prep for simulations
# declare Data for LD
#===================
threshold <- 2000
dist_mat <- nearest.dist(x=sp, y=NULL, method="euclidian", delta=threshold, upper=NULL);          
dist_mat <- as.matrix(dist_mat)

# cumulProbMat <- fast_prob_mat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, threshold, sp, dist_mat,SB,cumul=TRUE)
cumulProbMat <- generate_prob_mat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, threshold, sp, dist_mat,SB,cumul=TRUE)

# test cumulProbMat
test_that("cumulProbMat is cumulative, ending in ones",{
expect_equal(sum(cumulProbMat[dim(cumulProbMat)[1],]),dim(cumulProbMat)[1])
expect_true(all.equal(cumulProbMat[dim(cumulProbMat)[1],],rep(1,dim(cumulProbMat)[1])))
})
# all.equal(cumulProbMatRef,cumulProbMat)

### the vec of stats for the data
infestH<-which(maps$infest3==1)
stats <- multiGilStat(cumulProbMat=cumulProbMat, blockIndex = blockIndex, infestH = infestH, timeH=maps$ages[infestH], endTime = 1, rateMove = rateMove, Nrep = Nrep, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=FALSE)
if(!is.vector(stats$statsTable)){
statsData <- stats$statsTable[, 1]
}else{
statsData <- stats$statsTable
}

### the standard call to multiGilStat
snapTimes<-attributes(maps)$snapTimes
nbit<-snapTimes[3]-snapTimes[2]
infestH<-which(maps$infest2==1)
start<- Sys.time()
outBase <- multiGilStat(cumulProbMat=cumulProbMat, blockIndex = blockIndex, infestH = infestH, timeH=rep(-1,length(infestH)), endTime = nbit, rateMove = rateMove, Nrep = Nrep, coords = maps[, c("X", "Y")], breaksGenVar = genIntervals, simul=TRUE)
cat(Sys.time()-start)

# # test commentes as following would only be true with enough Nrep
# test_that("multiGilStat with true params -> stats ~= statsData",{
# expect_true(count(apply(outBase$statsTable,1,max)<statsData)==0)
# expect_true(count(apply(outBase$statsTable,1,min)>statsData)==0)
# })

### starting point for simulations
infestH <- which(maps$infest2 > 0)
timeH <- maps$ages[infestH]

### Priors (also the place to change the parameters)
priorMeans<-c(0.03, 0.04, 0.35, 7, 40)
priorSdlog <- c(1, 1, 1, 1, 1)
realMeans<-c(rateMove, rateJumpInMove, delta, halfDistH, halfDistJ)
sampling <- c("lnorm", "lnorm", "lnorm", "lnorm", "lnorm")
names(priorMeans)<-c("rateMove", "rateJumpInMove", "delta", "halfDistH", "halfDistJ")
names(sampling)<-names(priorMeans)
names(realMeans)<-names(priorMeans)
names(priorSdlog)<-names(priorMeans)

### Intervals of definition for the parameters
paramInf<-c(0.002,0.0001,0.0001, 0.5, 25)
paramSup<-c(0.30,0.30,0.7, 20, 100)
names(paramInf)<-names(priorMeans)
names(paramSup)<-names(priorMeans)

### LD formalism for Data (no setting should be made past this declaration)
PGF<-function(Data){ # parameters generating functions (for init etc...)
	
	priorMeans<-Data$priorMeans
	priorSdlog<-Data$priorSdlog
	
	values<-rlnorm(length(priorMeans),mean=log(priorMeans), sd=priorSdlog)
	
	return(values)
}

### If cumulProbMat set to NULL multiGilStat will use theta to compute it
### Must calculate cumulProbMat in C each time
MyDataFullSample <- list(y=statsData,
	     trans=NULL,
	     cumulProbMat=NULL,
	     blockIndex=blockIndex,
	     dist_out = makeDistClassesWithStreets(X = as.vector(maps[, "X"]), Y = as.vector(maps[, "Y"]), genIntervals, blockIndex), 
	     infestH=infestH,
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
	     mon.names=c("ll","lp", names(priorMeans)), # monitored variables (like in Model)
	     parm.names=names(priorMeans), # parameters names (like in Model and Initial.Values)
	     sampling=sampling # method of sampling parameters
		)

rm(cumulProbMat,dist_mat)

### "Pirat" data (global because need to be updated)
minLLever=-10e6

#==================================
## declare Model for LD
#==================================
Model <- function(theta,Data,postDraw=FALSE){
	# coerce theta, a priori all positive
	theta<-theta
	names(theta)<-Data$parm.names
	
	# set intervals for the variables
    for(param in names(theta)){
        theta[param]<-interval(theta[param],a=Data$paramInf[param],b=Data$paramSup[param])
    }
	
	# cat("theta:",theta, "\n")

	seed<-runif(1,min=0,(2^16-1)) # NB: for repeatability, R seed 
					#     Should be set

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
		LP<-LL+sum(dlnorm(theta, meanlog=log(Data$priorMeans), sdlog=Data$priorSdlog, log = TRUE))
		# LP<-LL+sum(dnorm(log(theta),mean=log(Data$priorMeans),sd=1, log = TRUE))
		# cat("LL:",LL,"LP:",LP,"\n")
	}

	# return
	Modelout <- list(LP=LP, # joint posterior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP, theta), # to be monitored/ploted
			 yhat=yhat, # data generated for that set of parameter
			 	    # will be used for posterior check
			 parm=theta # the parameters, possibly constrained by the model
			 )

	return(Modelout)
}

 
#===========================
# Testing Model/Data
#===========================
start<-Sys.time()
ModelOutGood<-Model(priorMeans,MyDataFullSample)
cat(Sys.time()-start)
ModelOutBest<-Model(realMeans,MyDataFullSample)
expect_true(ModelOutGood$Dev>ModelOutBest$Dev-4)

# theta<-realMeans
# theta["rateMove"]<-log(0.5)
# ModelOutBest<-Model(theta,MyDataFullSample)

# testRateMove<-seq(0.025,0.05,0.0005)
# testLL<-0*testRateMove
# for(i in 1:length(testRateMove)){
#	rM<-testRateMove[i]
#	cat("RateMove:",rM," ")
#	theta["rateMove"]<-log(rM)
#	ModelOutGood<-Model(theta=theta["rateMove"],MyDataFullSample)
#	testLL[i]<-ModelOutGood$LP
# }

signedLog<-function(signedBigNums){
	signNum<-sign(signedBigNums)
	signedLog<-log(abs(signedBigNums))*signNum
	return(signedLog)
}

# stop()

# par(mfrow=c(1,2))
# plot(log(testRateMove),signedLog(testLL),type="l",ylab="Log(Log(likelihood))",xlab="Log(rateMove)",main="Synthetic likelihood profile for rateMove \n blue is true");
# abline(v=log(rateMove),col="blue")
# int<-0.5
# sel<-which(log(testRateMove) > log(rateMove)-int & log(testRateMove) < log(rateMove)+int)
# sel<-which(testLL > -0)
# plot(testRateMove[sel],exp(testLL[sel]),type="l",ylab="Log(likelihood)",xlab="Log(rateMove)",main="Synthetic likelihood profile for rateMove \n blue is true (Zoom)");
# abline(v=rateMove,col="blue")
# dev.print(device=pdf,"SynLikProfileRateMove.pdf")
# save.image("profilSynLLH_rateMove.img")



# weibull order plotting to check if stats are normal
# for(i in 1:dim(outBase$statsTable)[1])
# {
#  
# order <- order(outBase$statsTable[i, ])
# plot((1:dim(outBase$statsTable)[2])/(dim(outBase$statsTable)[2]+1), outBase$statsTable[i, order], main = i)
# Sys.sleep(0.5)
# 
# }

stop()

#===========================
# Init values 
#===========================

# Initial.Values <- GIV(Model, MyData, PGF=TRUE,n=100) #GIV: generate initial values

Initial.Values <- priorMeans
theta <- Initial.Values
nparams <-length(theta)

nbsimul <- 20000 #starting simulation size

upFreq <- 1 #update frequency
saveFreq <- 1000 #save frequency

sdprop <- rep(0.4, nparams)
names(sdprop) <- MyDataFullSample$parm.names
adaptOK <- FALSE
checkAdapt <- 20
beginEstimate <- -1
lowAcceptRate <- 0.15
highAcceptRate <- 0.4

useAutoStop <- TRUE
checkAutoStop <- 100 #initial length of time after which to check the autostop

# init of theta attributes and saving scheme
outModel<-Model(theta,MyDataFullSample)
Monitor<-mat.or.vec(nbsimul+1,length(MyDataFullSample$mon.names))
Monitor[1,]<-outModel$Monitor
attributes(theta)$outModel<-outModel

accepts<-as.data.frame(matrix(rep(0,nparams*nbsimul),ncol=nparams))
names(accepts)<-MyDataFullSample$parm.names


#============================
# Launch sampler from functions sampling
#============================

	Rprof()
	
	numit <- 1

	while(numit < nbsimul){
		
		## display at every update frequency
		if(upFreq!=0 && numit%%upFreq==0){
			  cat("it:",numit,"of", nbsimul, "current theta:", theta,"\n");
		}

		##save at every save frequency
		if(saveFreq!=0 && numit%%saveFreq==0){
			write.table(Monitor[(numit-saveFreq):numit, ], "thetasamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			write.table(accepts[numit-saveFreq):numit, ], "acceptsamples.txt", sep ="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			cat("numit: ",numit,"\nnbsimul: ",nbsimul,"\nadaptOK :",adaptOK,"\ncheckAdapt: ",checkAdapt,"\nsdprop: ", sdprop,"\nbeginEstimate: ",beginEstimate,"\nuseAutoStop: ",useAutoStop,"\ncheckAutoStop: ",checkAutoStop,"\nsaveFreq: ",saveFreq, "\n", file = "mcmc_values.txt", append = TRUE)
		
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
			if(!adaptOK && numit+checkAdapt > nbsimul)
			{
				
				cat("sampling variance adaptation not complete: numit: ", numit, "doubling number simulations\n")
				nbsimul <- nbsimul * 2
				Monitor<-resized(Monitor,nr=nbsimul+1)
				accepts<-resized(accepts,nr=nbsimul)

			}
		}
		
		## sample the variables
		for(paramName in MyDataFullSample$parm.names){
	      		theta<-omniSample(Model,MyDataFullSample,theta,paramName,sdprop[paramName])
	      		accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    	 }

	    	Monitor[numit+1,]<-attributes(theta)$outModel$Monitor

		## auto stopping
		if(useAutoStop && adaptOK && numit == beginEstimate + checkAutoStop){

			cb<-cb.diag(Monitor[(1+beginEstimate):numit, ],logfile="convergence_tests.txt")
			checkAutoStop<-min(cb$newNbIt,numit*3)
			
			if(!cb$ok){
					cat("checking auto stop: numit: ", numit, "next check: ", numit + checkAutoStop)
					nbsimul <- beginEstimate + checkAutoStop
					Monitor<-resized(Monitor,nr=nbsimul+1)
					accepts<-resized(accepts,nr=nbsimul)
			}
		}	
	
		# increase the iteration count
		numit <- numit + 1
	}

	write.table(Monitor, "thetasamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
	write.table(accepts, "acceptsamples.txt", sep ="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
	
	Rprof(NULL)

#===========================
## Calculate estimates from MCMC
#===========================

# post treatment

Monitor<-as.data.frame(Monitor)
names(Monitor)<-MyDataFullSample$mon.names
burn.in<-ceiling(n.mc/5)

estMeans <- apply(Monitor[-(1:burn.in), -(1:2)], 2, mean)
names(estMeans) <- MyDataFullSample$parm.names
yMean<-mean(MyDataFullSample$y)
ySd<-sd(MyDataFullSample$y)


stop()

#===========================
# Below is older LD framework for MCMC
#	No longer used
#===========================

#===========================
## launch LD
#===========================

start1<-Sys.time()
Fit1 <- LaplacesDemon(Model, 
		     Data=MyData, 
		     Initial.Values, 
		     Covar=NULL, # proposal covariance matrix
		     Iterations=n.mc, # total number of iterations
		     Status=10, # how often status plotted 
		     Thinning=1,   # allow to save memory...
		     Algorithm="",  # see details for this and specs in ?LaplacesDemon
		     Specs = list(SIV = NULL, n1 = 4, at = 6, aw = 1.5) # standard recommended specifications for twalk
		     )

#===========================
## results
#===========================
print(Fit1)
Consort(Fit1)

# PosteriorChecks(Fit1)
MyData$infestH<-which(maps$infest3==1)
postMap<-getPosteriorMapsLD(Fit1,MyData,repByTheta=10)

# Brier scores
predReg<-maps$predict4
predReg[which(maps$infest3==1)]<-1

verWood<-verify(maps$infest4,maps$postMap)
# b<-verify(mapsholds<-exp(seq(log(0.01),log(1),0.2))
verReg<-verify(maps$infest4,maps$predict4)

dev.new()
probthresholds<-exp(seq(log(0.01),log(1),0.2))
par(mfrow=c(2,3))
plot(verReg)
rocReg<-roc.plot(maps$infest4,maps$predict4,main="predict4",CI=TRUE,n.boot=100)
valueReg<-value(maps$infest4,maps$predict4,thresholds=probthresholds,all=TRUE)

plot(verWood)
rocWood<-roc.plot(maps$infest4,maps$postMap,main="postMap",CI=TRUE,n.boot=100)
valueReg<-value(maps$infest4,maps$postMap,all=TRUE,thresholds=probthresholds)

dev.new()
discrimination.plot(maps$infest4,maps$predict4,main="Reg")
dev.new()
discrimination.plot(maps$infest4,maps$postMap,main="Wood")

dev.new()
par(mfrow=c(3,3))
plot_reel(maps$X,maps$Y,maps$infest3,base=0)
plot_reel(maps$X,maps$Y,maps$infest3,base=0,main="Init")
plot_reel(maps$X,maps$Y,maps$infest3,base=0)
plot_reel(maps$X,maps$Y,maps$postMap,base=0,main=paste("Posterior Map \nBrier:",verWood$bs,"Skill:",verWood$ss))
plot_reel(maps$X,maps$Y,maps$infest4,base=0,main="Observed")
plot_reel(maps$X,maps$Y,maps$predict4,base=0,main=paste("Prediction regression\nBrier:",verReg$bs,"Skill:",verReg$ss))
plot_reel(maps$X,maps$Y,maps$postMap,base=0,top=0.1)
plot_reel(maps$X,maps$Y,maps$infest4,base=0,top=0.1,main="idem topped at 0.1")
plot_reel(maps$X,maps$Y,maps$predict4,base=0,top=0.1)

caterpillar.plot(Fit1, Parms="rateMove")
BurnIn <- Fit1$Rec.BurnIn.Thinned
plot(Fit1, BurnIn, MyData, PDF=FALSE)
hist(exp(Fit1$Posterior1))
Pred <- predict(Fit1, Model, MyData)
summary(Pred, Discrep="Chi-Square")
plot(Pred, Style="Covariates", Data=MyData)
plot(Pred, Style="Density", Rows=1:9)
plot(Pred, Style="ECDF")
plot(Pred, Style="Fitted")
plot(Pred, Style="Jarque-Bera")
plot(Pred, Style="Predictive Quantiles")
plot(Pred, Style="Residual Density")
plot(Pred, Style="Residuals")
Levene.Test(Pred)
Importance(Fit1, Model, MyData, Discrep="Chi-Square")

# prediction
