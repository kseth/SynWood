# aim is to find the worse case scenario for the regression prediction
# compared to best Predictions
# and for a bunch of different simulations

source("param.r")
source("RanalysisFunctions.R")
source("functions_migration.R")
library(verification)

snapTimes<-c(0,52*c(3,5,7))
nbSimuls<-3
EvalOnDifferentSimul<-FALSE 
     # if FALSE:
     	# training is done on 2 and 3rd snapTimes
     	# prediction is done on 3 to 4th snapTimes
     # if TRUE
     	# training is done on 2 and 3rd snapTimes
     	# prediction is done on 2 to 3rd snapTimes of other simul

#===========================
## Set up
#===========================
scale <- 1/rateMove
threshold <- 2000
dist_mat <-nearest.dist(x=sp, y=NULL, method="euclidian", delta=threshold, upper=NULL);          
dist_mat <- as.matrix(dist_mat)
cumulProbMat<- generate_prob_mat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, threshold, sp, dist_mat,SB,cumul=TRUE)
infestH <- which(infestHouse>0) # infestH - locations of all infested houses
timeH <- c(-1)# investation a t-1, easier to manage after
multiSeeds<-list()

#===========================
## simulations
#===========================

for(numSimul in 1:nbSimuls){
  cat("numsimul:",numSimul,"/",nbSimuls,"\n")
  seed<-numSimul
set.seed(seed)
out <- gillespie(probMat=NULL, cumulProbMat, infestH, timeH,endTime=max(snapTimes), scale,seed=seed)

#==============================
# make maps at different times
#==============================
nsnap<-length(snapTimes)
maps<-infestSerieToMaps(out,sp,snapTimes[nsnap])
maps<-changeNameCol(maps,"infest",paste("infest",nsnap,sep=""))

while(nsnap>0){
	colName<-paste("infest",nsnap,sep="")
	maps[,colName]<-infestSerieToMaps(out,sp,snapTimes[nsnap])$infest
	nsnap=nsnap-1
}

#===========================
# simulation of missing data
#===========================
maps$blockIndex<-blockIndex
detectRate<-90/100
openRate<-70/100
maps$obs2<-simulObserved(maps$infest2,openRate,detectRate)
maps$obs3<-simulObserved(maps$infest3,openRate,detectRate)
maps$obs4<-simulObserved(maps$infest4,openRate,detectRate)

#===========================
# basic model
#===========================
basicModelImp<-getBasicModel(maps,"obs2","obs3") # training

#===========================
# Data to use for prediction/check
#===========================
if(EvalOnDifferentSimul){
  # new simul
  out <- gillespie(probMat=NULL, cumulProbMat, infestH, timeH,endTime=max(snapTimes), scale,seed=seed+1) 
  maps$initInf<-infestSerieToMaps(out,sp,snapTimes[2])$infest
  maps$initObs<-simulObserved(maps$initInf,openRate,detectRate)
  maps$toPred<-infestSerieToMaps(out,sp,snapTimes[3])$infest
}else{
  maps$initInf<-maps$infest3
  maps$initObs<-maps$obs3
  maps$toPred<-maps$infest4
}

#============================
# Predictions
#============================
maps$predictGlm<-predictBasicModel(maps,"initInf",model=basicModelImp)
maps$predictGlmImp<-predictBasicModel(maps,"initObs",model=basicModelImp)
maps$predictBest<-getBestPredict(maps,cumulProbMat,maps$initInf)
maps$predictBestImp<-getBestPredict(maps,cumulProbMat,maps$initObs)
# Perfect predictions over know to be imperfect data
# estInitProb<-IntelliSmooth(maps,maps$obs3)
# maps$bestPredictImpC4<-getBestPredict(maps,cumulProbMat,estInitProb)

#============================
# Assessment
#============================
bsglm<-verify(maps$toPred,maps$predictGlm,show=FALSE)$bs
bsglmImp<-verify(maps$toPred,maps$predictGlmImp,show=FALSE)$bs
bsbest<-verify(maps$toPred,maps$predictBest,show=FALSE)$bs
bsbestImp<-verify(maps$toPred,maps$predictBestImp,show=FALSE)$bs

percWorse<- (bsglm-bsbest)/bsbest
percWorseImp<- (bsglmImp-bsbestImp)/bsbestImp
cat("% glm worse compared to perfect",percWorse," perfect from imp. obs",percWorseImp,"\n")

dev.new()
par(mfrow=c(2,3))
plot_reel(maps$X,maps$Y,maps$toPred,base=0,main="To predict")
plot_reel(maps$X,maps$Y,maps$predictGlmImp,base=0,
	  main=paste("Predict forward\nBSinf:",
		     bsglmImp,
	  "BSobs:",verify(maps$obs4,maps$predictGlmImp,show=FALSE)$bs))


plot_reel(maps$X,maps$Y,maps$predictBest,base=0,
	  main=paste("Predict forward\nBSinf:",
		     bsbest,
	  "BSobs:",verify(maps$obs4,maps$predictBest,show=FALSE)$bs))

plot_reel(maps$X,maps$Y,maps$toPred,base=0,main="To predict",top=0.1)
plot_reel(maps$X,maps$Y,maps$predictGlmImp,base=0,top=0.1,
	  main=paste("Predict forward\nBSinf:",
		     bsglmImp,
	  "BSobs:",verify(maps$obs4,maps$predictGlmImp,show=FALSE)$bs))


plot_reel(maps$X,maps$Y,maps$predictBest,base=0,top=0.1,
	  main=paste("Predict forward\nBSinf:",
		     bsbest,
	  "BSobs:",verify(maps$obs4,maps$predictBest,show=FALSE)$bs))



#============================
# Saving
#============================
multiSeeds$seed[numSimul]<-seed
multiSeeds$bsglmImp[numSimul]<-bsglmImp
multiSeeds$bsbest[numSimul]<-bsbest
multiSeeds$bsbestImp[numSimul]<-bsbestImp
multiSeeds$percWorse[numSimul]<-percWorse
multiSeeds$percWorseImp[numSimul]<-percWorseImp
}
dev.new()
par(mfrow=c(2,3))
plot(1,1,type="n")
plot(multiSeeds$percWorse,ylim=c(-1,1))
plot(multiSeeds$percWorseImp,ylim=c(-1,1))
plot(multiSeeds$bsglmImp,ylim=c(0,1))
plot(multiSeeds$bsbestImp,ylim=c(0,1))
plot(multiSeeds$bsbest,ylim=c(0,1))

cat("mean percWorse:",mean(multiSeeds$percWorse),
    "\nmean bs glm :",mean(multiSeeds$bsglmImp),"sd:",sd(multiSeeds$bsglmImp),
    "\nmean bs best:",mean(multiSeeds$bsbest),"sd:",sd(multiSeeds$bsbest),"\n")
cat("mean percWorseImp:",mean(multiSeeds$percWorseImp),
    "\nmean bs glm     :",mean(multiSeeds$bsglmImp),"sd:",sd(multiSeeds$bsglmImp),
    "\nmean bs best Imp:",mean(multiSeeds$bsbestImp),"sd:",sd(multiSeeds$bsbestImp),"\n")

