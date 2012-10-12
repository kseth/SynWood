# aim is to find the worse case scenario for the regression prediction
# compared to best Predictions
# and for a bunch of different simulations

source("param.r")
source("RanalysisFunctions.R")
source("functions_migration.R")
library(verification)

#===========================
## Parameters
#===========================
trainTimes<-c(52*c(3,5))
predTimes<-c(52*c(3,5))
EvalOnSameSimul<-FALSE 
displayGraphs<-FALSE
waitForDisplay<-FALSE # if TRUE wait for entry on command line before continuing
compToSL<-TRUE # if TRUE use thetas.R to generate postMap
	       # if FALSE use the real model to generate the postMap

nbSimuls<-30
openRate<-70/100
detectRate<-90/100

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

thetas<-rnorm(100,mean=rateMove,sd=0.003) # according to grid sampling
thetas<-thetas[which(thetas>0)]

#===========================
## simulations
#===========================

for(numSimul in 1:nbSimuls){
  cat("numsimul:",numSimul,"/",nbSimuls,"\n")
  seed<-numSimul
  set.seed(seed)
  #=============================
  # training dataset
  #=============================
  outTrain <- gillespie(probMat=NULL, cumulProbMat, infestH, timeH,endTime=max(trainTimes,predTimes), scale,seed=seed)

  #==============================
  # make maps of training
  #==============================
  # simulations infestation
  maps<-infestSerieToMaps(outTrain,sp,trainTimes[1])
  maps<-changeNameCol(maps,"infest","initTrainInf")

  maps$endTrainInf<-infestSerieToMaps(outTrain,sp,trainTimes[2])$infest

  maps$blockIndex<-blockIndex

  # simulation of missing data
  maps$initTrainObs<-simulObserved(maps$initTrainInf,openRate,detectRate)
  maps$endTrainObs<-simulObserved(maps$endTrainInf,openRate,detectRate)

  #===========================
  # train glm
  #===========================
  basicModel<-getBasicModel(maps,"initTrainInf","endTrainInf") # training
  basicModelImp<-getBasicModel(maps,"initTrainObs","endTrainObs") # training

  #===========================
  # Data to use for prediction/check
  #===========================
  if(EvalOnSameSimul){
    # use training simulation
    outPred<-outTrain
  }else{
    # new simul
    outPred <- gillespie(probMat=NULL, cumulProbMat, infestH, timeH,endTime=max(predTimes), scale,seed=seed+1) 
  }
  maps$initInf<-infestSerieToMaps(outPred,sp,predTimes[1])$infest
  maps$initObs<-simulObserved(maps$initInf,openRate,detectRate)
  maps$toPred<-infestSerieToMaps(outPred,sp,predTimes[2])$infest

  #============================
  # Predictions
  #============================
  maps$predictGlm<-predictBasicModel(maps,"initInf",model=basicModel)
  maps$predictGlmImp<-predictBasicModel(maps,"initObs",model=basicModelImp)
  if(compToSL){
	  # cumulProbMat <- generate_prob_mat_C(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, dist_mat, blockIndex, cumul=TRUE)

	  maps$predictSL<-getPosteriorMaps(maps,thetas,
			cumulProbMat=cumulProbMat,
			maps$initInf,nbit=predTimes[2]-predTimes[1],
			repByTheta=100)

	  maps$predictSLImp<-getPosteriorMaps(maps,thetas,
			cumulProbMat=cumulProbMat,
			maps$initObs,nbit=predTimes[2]-predTimes[1],
			repByTheta=100)
  }else{
	  maps$predictSL<-NA
	  maps$predictSLImp<-NA
  }
  maps$predictBest<-getBestPredict(maps,cumulProbMat,maps$initInf,nbit=predTimes[2]-predTimes[1])
  maps$predictBestImp<-getBestPredict(maps,cumulProbMat,maps$initObs,nbit=predTimes[2]-predTimes[1])
  # # Perfect predictions over know to be imperfect data (to do)
  # estInitProb<-IntelliSmooth(maps,maps$initObs)
  # maps$predictBestImpC<-getBestPredict(maps,cumulProbMat,estInitProb)

  #============================
  # Assessment
  #============================
  bsglm<-verify(maps$toPred,maps$predictGlm,show=FALSE)$bs
  bsglmImp<-verify(maps$toPred,maps$predictGlmImp,show=FALSE)$bs
  bsSL<-verify(maps$toPred,maps$predictSL,show=FALSE)$bs
  bsSLImp<-verify(maps$toPred,maps$predictSLImp,show=FALSE)$bs
  bsbest<-verify(maps$toPred,maps$predictBest,show=FALSE)$bs
  bsbestImp<-verify(maps$toPred,maps$predictBestImp,show=FALSE)$bs

  percglmSL<- (bsglm-bsSL)/bsSL
  percglmSLImp<- (bsglmImp-bsSLImp)/bsSLImp
  percWorse<- (bsglm-bsbest)/bsbest
  percWorseImp<- (bsglmImp-bsbestImp)/bsbestImp
  cat("% glm worse compared to perfect",percWorse," perfect from imp. obs",percWorseImp,"comp to SL on inf",percglmSL,"to SL on Obs",percglmSLImp,"\n")

  if(displayGraphs){
    # dev.new()
    graphics.off()
    # training
    par(mfcol=c(2,4))
    plot_reel(maps$X,maps$Y,maps$initTrainInf,base=0,main="initTrainInf")
    plot_reel(maps$X,maps$Y,maps$endTrainInf,base=0,main="endTrainInf")
    plot_reel(maps$X,maps$Y,maps$initTrainObs,base=0,main="initTrainObs")
    plot_reel(maps$X,maps$Y,maps$endTrainObs,base=0,main="endTrainObs")

    # prediction
    plot_reel(maps$X,maps$Y,maps$initInf,base=0,main="initInf")
    plot_reel(maps$X,maps$Y,maps$toPred,base=0,main="To predict")
    plot_reel(maps$X,maps$Y,maps$initObs,base=0,main="initObs")

    dev.new()
    if(compToSL) ncolDisp<-6 else ncolDisp<-4

    par(mfrow=c(2,ncolDisp))
    for(top in c(1,0.1)){
      plot_reel(maps$X,maps$Y,maps$predictGlm,base=0,top=top,
		main=paste("glm predict on Inf\nBSinf:",
			   bsglm,
			   "BSobs:",verify(maps$toPred,maps$predictGlm,show=FALSE)$bs))
      plot_reel(maps$X,maps$Y,maps$predictGlmImp,base=0,top=top,
		main=paste("glm predict on Obs\nBSinf:",
			   bsglmImp,
			   "BSobs:",verify(maps$toPred,maps$predictGlmImp,show=FALSE)$bs))
      plot_reel(maps$X,maps$Y,maps$predictBest,base=0,top=top,
		main=paste("Perfect predict on Inf\nBSinf:",
			   bsbest,
			   "BSobs:",verify(maps$toPred,maps$predictBest,show=FALSE)$bs))
      plot_reel(maps$X,maps$Y,maps$predictBestImp,base=0,top=top,
		main=paste("Perfect predict on Obs\nBSinf:",
			   bsbestImp,
			   "BSobs:",verify(maps$toPred,maps$predictBestImp,show=FALSE)$bs))
      if(compToSL) {
	plot_reel(maps$X,maps$Y,maps$predictSL,base=0,top=top,
		  main=paste("SL predict on Inf\nBSinf:",
			     bsSL,
			     "BSobs:",verify(maps$toPred,maps$predictSL,show=FALSE)$bs))
	plot_reel(maps$X,maps$Y,maps$predictSLImp,base=0,top=top,
		  main=paste("SL predict on Obs\nBSinf:",
			     bsSLImp,
			     "BSobs:",verify(maps$toPred,maps$predictSLImp,show=FALSE)$bs))

      }
    }

    if(waitForDisplay) readline()
  }

  #============================
  # Saving
  #============================
  multiSeeds$seed[numSimul]<-seed
  multiSeeds$bsglm[numSimul]<-bsglm
  multiSeeds$bsglmImp[numSimul]<-bsglmImp
  multiSeeds$bsbest[numSimul]<-bsbest
  multiSeeds$bsbestImp[numSimul]<-bsbestImp
  multiSeeds$bsSL[numSimul]<-bsSL
  multiSeeds$bsSLImp[numSimul]<-bsSLImp
  multiSeeds$percglmSL[numSimul]<-percglmSL
  multiSeeds$percglmSLImp[numSimul]<-percglmSLImp
  multiSeeds$percWorse[numSimul]<-percWorse
  multiSeeds$percWorseImp[numSimul]<-percWorseImp
}
dev.new()
par(mfcol=c(2,2))
plot(multiSeeds$percWorse,ylim=c(-1,1),main="glm/best")
plot(multiSeeds$percWorseImp,ylim=c(-1,1))
plot(multiSeeds$percglmSL,ylim=c(-1,1),main="glm/SL")
plot(multiSeeds$percglmSLImp,ylim=c(-1,1))

par(mfcol=c(2,3))
plot(multiSeeds$bsglm,ylim=c(0,1),main="glm")
plot(multiSeeds$bsglmImp,ylim=c(0,1))
plot(multiSeeds$bsbest,ylim=c(0,1),main="best")
plot(multiSeeds$bsbestImp,ylim=c(0,1))
plot(multiSeeds$bsSL,ylim=c(0,1),main="SL")
plot(multiSeeds$bsSLImp,ylim=c(0,1))

cat("mean percWorse:",mean(multiSeeds$percWorse),
    "\nmean bs glm :",mean(multiSeeds$bsglm),"sd:",sd(multiSeeds$bsglm),
    "\nmean bs best:",mean(multiSeeds$bsbest),"sd:",sd(multiSeeds$bsbest),"\n")

cat("mean percWorseImp:",mean(multiSeeds$percWorseImp),
    "\nmean bs glm Imp    :",mean(multiSeeds$bsglmImp),"sd:",sd(multiSeeds$bsglmImp),
    "\nmean bs best Imp:",mean(multiSeeds$bsbestImp),"sd:",sd(multiSeeds$bsbestImp),"\n")

cat("mean percglmSL:",mean(multiSeeds$percglmSL),
    "\nmean bs SL:",mean(multiSeeds$bsSL),"sd:",sd(multiSeeds$bsSL),"\n")

cat("mean percglmSLImp:",mean(multiSeeds$percglmSLImp),
    "\nmean bs SL Imp:",mean(multiSeeds$bsSLImp),"sd:",sd(multiSeeds$bsSLImp),"\n")

percImproveBestGlm<-(multiSeeds$bsglm-multiSeeds$bsbest)/multiSeeds$bsglm
percImproveBestGlmImp<-(multiSeeds$bsglmImp-multiSeeds$bsbestImp)/multiSeeds$bsglmImp
percImproveSLGlm<-(multiSeeds$bsglm-multiSeeds$bsSL)/multiSeeds$bsglm
percImproveSLGlmImp<-(multiSeeds$bsglmImp-multiSeeds$bsglmImp)/multiSeeds$bsSLImp
cat("mean percImprove Best/glm:",mean(percImproveBestGlm))
cat("mean percImprove Best/glm Imp:",mean(percImproveBestGlmImp))
cat("mean percImprove SL/glm:",mean(percImproveSLGlm))
cat("mean percImprove SL/glm Imp:",mean(percImproveSLGlmImp))

## tests of difference of mean
par(mfrow=c(1,2))
plot(density(multiSeeds$bsSL),col="blue",main="Perfect Obs\nBlue SL; black glm")
lines(density(multiSeeds$bsglm))
plot(density(multiSeeds$bsSLImp),col="blue",main="Imperfect Obs\nBlue SL; black glm")
lines(density(multiSeeds$bsglmImp))
# => seems reasonably bell shaped, can use t-test (Welch t-test)
print(t.test(multiSeeds$bsglm,multiSeeds$bsSL))
print(t.test(multiSeeds$bsglmImp,multiSeeds$bsSLImp))

# save.image("basic_reg_at_faults_diffSit_sameHighTime.img")

# # For final plots/graphs
# # keep worse case for glm
# worseCase<-order(percImproveSLGlmImp)[length(percImproveSLGlmImp)]
# # => relaunch the loop only on that Case
# 
# 
# # using plot_reel
# dev.new()
# par(mfrow=c(1,3))
# plot_reel(maps$X,maps$Y,maps$initInf,base=0,main="initInf")
# plot_reel(maps$X,maps$Y,maps$toPred,base=0,main="To predict")
# plot_reel(maps$X,maps$Y,maps$initObs,base=0,main="initObs")
# 
# # dev.print(device=pdf,"image_data.pdf")
# # dev.print(device=png,"image_data.png",width=1200)
# 
# dev.new()
# par(mfrow=c(2,2))
# top=1
# plot_reel(maps$X,maps$Y,maps$predictGlm,base=0,top=top,
# 	  main=paste("glm predict on Inf\nBSinf:",
# 		     bsglm,
# 		     "BSobs:",verify(maps$toPred,maps$predictGlm,show=FALSE)$bs))
# plot_reel(maps$X,maps$Y,maps$predictGlmImp,base=0,top=top,
# 	  main=paste("glm predict on Obs\nBSinf:",
# 		     bsglmImp,
# 		     "BSobs:",verify(maps$toPred,maps$predictGlmImp,show=FALSE)$bs))
# plot_reel(maps$X,maps$Y,maps$predictSL,base=0,top=top,
# 	  main=paste("SL predict on Inf\nBSinf:",
# 		     bsSL,
# 		     "BSobs:",verify(maps$toPred,maps$predictSL,show=FALSE)$bs))
# plot_reel(maps$X,maps$Y,maps$predictSLImp,base=0,top=top,
# 	  main=paste("SL predict on Obs\nBSinf:",
# 		     bsSLImp,
# 		     "BSobs:",verify(maps$toPred,maps$predictSLImp,show=FALSE)$bs))
# # dev.print(device=pdf,"plot_reel_predictions.pdf")
# # dev.print(device=png,"plot_reel_predictions.png",width=1200)
# 
# # using image
# par(mfrow=c(1,3))
# nCote<-sqrt(length(maps$initInf))
# colors<-yellow.colors(100)
# image(matrix(maps$initInf,ncol=nCote),zlim=c(0,1),col=colors,asp=1,main="Init")
# image(matrix(maps$toPred,ncol=nCote),zlim=c(0,1),col=colors,asp=1,main="To Pred")
# image(matrix(maps$initObs,ncol=nCote),zlim=c(0,1),col=colors,asp=1,main="Init Obs")
# # dev.print(device=pdf,"image_data_yellow.pdf")
# # dev.print(device=png,"image_data_yellow.png",width=1200)
# 
# par(mfrow=c(2,2))
# image(matrix(maps$predictGlm,ncol=sqrt(length(maps$predictGlm))),col=colors,zlim=c(0,1),asp=1,main="Glm")
# image(matrix(maps$predictSL,ncol=sqrt(length(maps$predictSL))),zlim=c(0,1),col=colors,asp=1,main="SL")
# image(matrix(maps$predictGlmImp,ncol=sqrt(length(maps$predictSLImp))),zlim=c(0,1),col=colors,asp=1,main="Glm Imp")
# image(matrix(maps$predictSLImp,ncol=sqrt(length(maps$predictSLImp))),zlim=c(0,1),asp=1,col=colors,main="SL Imp")
# # dev.print(device=pdf,"image_predictions_yellow.pdf")
# # dev.print(device=png,"image_predictions_yellow.png",width=1200)


