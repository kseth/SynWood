library(sl)
library(LaplacesDemon)
library(testthat)
library(verification)
source("RanalysisFunctions.R")
source("logLik.r") # override some sl functions and add synLik

source("param.r")

## load data 
## use simulated data as in base_regression.R
source("maps_basic_regression.R") # gives maps
source("functions_migration.R")

#==================
#  Params, if modifying stats: look for multiGilStat
#==================
Nrep=300;
n.mc=30000
set.seed(1)
genIntervals <- seq(0, 250, 30) # distance classes for the general variogram

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

### the vector of stats for the data
infestH<-which(maps$infest3==1)
stats <- multiGilStat(cumulProbMat=cumulProbMat, blockIndex, infestH, timeH=maps$ages[infestH], endTime = 1, rateMove, Nrep, coords = maps[, c("X", "Y")], breaks = genIntervals, simul=FALSE)
statsData <- stats$statsTable[, 1]

### the standard call to multiGilStat
snapTimes<-attributes(maps)$snapTimes
nbit<-snapTimes[3]-snapTimes[2]
infestH<-which(maps$infest2==1)
outBase <- multiGilStat(cumulProbMat=cumulProbMat, blockIndex, infestH, timeH=rep(-1,length(infestH)), endTime = nbit, rateMove, Nrep, coords = maps[, c("X", "Y")], breaks = genIntervals, simul=TRUE)
# test
test_that("multiGilStat with true params -> stats ~= statsData",{
expect_true(count(apply(outBase$statsTable,1,max)<statsData)==0)
expect_true(count(apply(outBase$statsTable,1,min)>statsData)==0)
})

### starting point for simulations
infestH <- which(maps$infest2 > 0)
timeH <- maps$ages[infestH]

### Priors (also the place to change the parameters)
priorMeans<-log(c(rateMove))
realMeans<-log(c(rateMove))
names(priorMeans)<-c("rateMove")
names(realMeans)<-names(priorMeans)

### Intervals of definition for the parameters
paramInf<-c(0.002)
paramSup<-c(0.24)
names(paramInf)<-names(priorMeans)
names(paramSup)<-names(priorMeans)

### LD formalism for Data (no setting should be made past this declaration)
PGF<-function(Data){ # parameters generating functions (for init etc...)
	priorMeans<-Data$priorMeans
  values<-rnorm(length(priorMeans),mean=priorMeans,sd=0.15)
  return(values)
}

MyData <- list(y=statsData,
	     trans=NULL,
	     cumulProbMat=cumulProbMat,
	     blockIndex=blockIndex,
	     infestH=infestH,
	     timeH=timeH,
	     endTime=nbit,
	     maps=maps,
	     nbit=nbit,
	     Nrep=Nrep,
	     priorMeans=priorMeans,
	     genIntervals=genIntervals,
	     paramInf=paramInf,
	     paramSup=paramSup,
	     PGF=PGF,
	     mon.names=c("ll","lp"), # monitored variables (like in Model)
	     parm.names=names(priorMeans)# parameters names (like in Model and Initial.Values)
	     )

rm(cumulProbMat,dist_mat)

### "Pirat" data (global because need to be updated)
minLLever=-10e6

#==================================
## declare Model for LD
#==================================
Model<-function(theta,Data,postDraw=FALSE,infestHints=NULL){
	# coerce theta, a priori all positive
	theta<-exp(theta)
	names(theta)<-Data$parm.names
	# theta["rateMove"]<-interval(theta["rateMove"],a=0.002,b=0.24)
	theta["rateMove"]<-interval(theta["rateMove"],a=Data$paramInf["rateMove"],b=Data$paramSup["rateMove"])
	names(theta)<-Data$parm.names
	cat("theta:",theta)

	# simulations
	start <- Sys.time()
	seed<-runif(1,min=0,(2^16-1)) # NB: for repeatability, R seed 
	if(postDraw){
		getStats<-FALSE
					#     Should be set
		
		# then in wrapper use infested if Nrep =1 and infested$Dens if Nrep=a few
		# to show a posterior predictive map
	}else{
		getStats<-TRUE
	}
	out <- multiGilStat(Data$cumulProbMat, Data$blockIndex, 
			    Data$infestH, Data$timeH, endTime = Data$nbit, 
			    theta["rateMove"], 
			    Data$Nrep, 
			    coords = Data$maps[, c("X", "Y")], 
			    breaks = Data$genIntervals,
			    getStats=getStats,
			    seed=seed
			    )
	end <- Sys.time()
	cat("t multiGil:",end-start,"\n")

	if(postDraw){
		yhat<-out$infestedDens
		LL<-NA
		LP<-NA
	}else{
		yhat<-out$statsTable[,1]

		statsToKeep<-c(1:16,39:42)
		statsTable<-out$statsTable[statsToKeep,]
		y<-Data$y[statsToKeep]

		# statsTable<-out$statsTable[!duplicated(yhat),]
	
		# synthetic likelihood
		success<-try(ll<-synLik(statsTable,y,Data$trans))
		if(class(success)=="try-error"){
		  ll<-minLLever
		}else{
		  minLLever<<-min(ll,minLLever)
		}

		# get likelihood with priors
		LL<-ll
		attributes(LL)<-NULL
		LP<-LL+sum(dnorm(log(theta),mean=Data$priorMeans,sd=1))
		# LP<-LL+sum(dlnorm(theta,meanlog=log(Data$priorMeans),sdlog=1))
		cat("LL:",LL,"LP:",LP,"\n")
	}

	# return
	Modelout <- list(LP=LP, # joint posterior
			 Dev=-2*LL, # deviance, probably not to be changed
			 Monitor=c(LL,LP), # to be monitored/ploted
			 yhat=yhat, # data generated for that set of parameter
			 # will be used for posterior check
			 parm=log(theta) # the parameters, possibly constrained by the model
			 )
} 
#===========================
# Testing Model/Data
#===========================
ModelOutGood<-Model(priorMeans,MyData)
ModelOutBest<-Model(realMeans,MyData)
expect_true(ModelOutGood$Dev>ModelOutBest$Dev-4)

theta<-realMeans
theta["rateMove"]<-log(0.5)
ModelOutBest<-Model(theta,MyData)

testRateMove<-seq(0.025,0.05,0.0005)
testLL<-0*testRateMove
for(i in 1:length(testRateMove)){
	rM<-testRateMove[i]
	cat("RateMove:",rM," ")
	theta["rateMove"]<-log(rM)
	ModelOutGood<-Model(theta,MyData)
	testLL[i]<-ModelOutGood$LP
}
signedLog<-function(signedBigNums){
	signNum<-sign(signedBigNums)
	signedLog<-log(abs(signedBigNums))*signNum
	return(signedLog)
}

par(mfrow=c(1,2))
plot(log(testRateMove),signedLog(testLL),type="l",ylab="Log(Log(likelihood))",xlab="Log(rateMove)",main="Synthetic likelihood profile for rateMove \n blue is true");
abline(v=log(rateMove),col="blue")
int<-0.5
# sel<-which(log(testRateMove) > log(rateMove)-int & log(testRateMove) < log(rateMove)+int)
sel<-which(testLL > -0)
plot(testRateMove[sel],exp(testLL[sel]),type="l",ylab="Log(likelihood)",xlab="Log(rateMove)",main="Synthetic likelihood profile for rateMove \n blue is true (Zoom)");
abline(v=rateMove,col="blue")
# dev.print(device=pdf,"SynLikProfileRateMove.pdf")
# save.image("profilSynLLH_rateMove.img")
stop()

#===========================
# Init values 
#===========================

# Initial.Values <- GIV(Model, MyData, PGF=TRUE,n=100) #GIV: generate initial values
Initial.Values <- priorMeans

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
stop1<-Sys.time()
cat("chain time:",stop1-start1)

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
