source("param.r")
source("RanalysisFunctions.R")
source("functions_migration.R")
# source("morans_functions.r")
library(geoR)

snapTimes<-c(0,52*c(1,3,5))
seed<-1

set.seed(seed)
scale <- 1/rateMove
threshold <- 2000
dist_mat <-nearest.dist(x=sp, y=NULL, method="euclidian", delta=threshold, upper=NULL);          
dist_mat <- as.matrix(dist_mat)
cumulProbMat<- generate_prob_mat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, threshold, sp, dist_mat,SB,cumul=TRUE)

#===========================
## simulation
#===========================
infestH <- which(infestHouse>0) # infestH - locations of all infested houses
timeH <- c(-1)# investation a t-1, easier to manage after
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
## visualisation simulation
#===========================
dev.new()
par(mfcol=c(2,length(snapTimes)))
for(itime in 1:length(snapTimes)){
	mapTime<-snapTimes[itime]
	mapsTemp<-infestSerieToMaps(out,sp,mapTime)
	plot_reel(mapsTemp$X,mapsTemp$Y,mapsTemp$infest,base=0,main=paste(mapTime,"time unit (weeks)"))
	plot_reel(mapsTemp$X,mapsTemp$Y,log(mapsTemp$ages+1),base=0,top=max(log(mapsTemp$ages+1)))
}
#===========================
# simulation of missing data
#===========================
detectRate<-90/100
openRate<-70/100
maps$obs2<-simulObserved(maps$infest2,openRate,detectRate)
maps$obs3<-simulObserved(maps$infest3,openRate,detectRate)
maps$obs4<-simulObserved(maps$infest4,openRate,detectRate)

#============================
# On perfect data
#============================
# train on switch 2->3
basicModel<-getBasicModel(maps,"infest2","infest3")

# training prediction
maps$predict3<-predictBasicModel(maps,"infest2")

# forward prediction
maps$predict4<-predictBasicModel(maps,"infest3")

dev.new()
par(mfrow=c(2,3))
plot_reel(maps$X,maps$Y,maps$infest2,base=0,main="Init observed")
plot_reel(maps$X,maps$Y,maps$infest3,base=0,main="2nd observed")
plot_reel(maps$X,maps$Y,maps$infest4,base=0,main="To predict")

plot_reel(maps$X,maps$Y,maps$infest2,base=0,main="Init")
plot_reel(maps$X,maps$Y,maps$predict3,base=0,main="Predict training")
plot_reel(maps$X,maps$Y,maps$predict4,base=0,main="Predict forward")

#============================
# On non perfect data
#============================
# train on switch 2->3
basicModelImp<-getBasicModel(maps,"obs2","obs3")

# training prediction
maps$predictImp3<-predictBasicModel(maps,"obs2",model=basicModelImp)

# forward prediction
maps$predictImp4<-predictBasicModel(maps,"obs3",model=basicModelImp)

#============================
# Perfect predictions using directly the model
#============================
maps$blockIndex<-blockIndex
maps$bestPredict4<-getBestPredict(maps,cumulProbMat,maps$infest3)
maps$bestPredict3<-getBestPredict(maps,cumulProbMat,maps$infest2)

#============================
# Plot Quality of predictions (+ Brier Score)
#============================
library(verification)

dev.new()
par(mfcol=c(3,5))
plot_reel(maps$X,maps$Y,maps$infest2,base=0,main="Init Infested")
plot_reel(maps$X,maps$Y,maps$infest3,base=0,main="2nd Infested")
plot_reel(maps$X,maps$Y,maps$infest4,base=0,main="To predict")

plot_reel(maps$X,maps$Y,maps$infest2,base=0,main="Init perf Obs")
plot_reel(maps$X,maps$Y,maps$predict3,base=0,main=paste("Predict training (perf obs)\nBS:",verify(maps$infest3,maps$predict3)$bs))
plot_reel(maps$X,maps$Y,maps$predict4,base=0,main=paste("Predict forward (perf Obs)\nBS:",verify(maps$infest4,maps$predict4)$bs))

plot_reel(maps$X,maps$Y,maps$obs2,base=0,main="Init Observed")
plot_reel(maps$X,maps$Y,maps$obs3,base=0,main="2nd Infested Obs")
plot_reel(maps$X,maps$Y,maps$obs4,base=0,main="Obs To predict")

plot_reel(maps$X,maps$Y,maps$obs2,base=0,main="Init Observed")
plot_reel(maps$X,maps$Y,maps$predictImp3,base=0,
	  main=paste("Predict training\nBSinf:",
		     verify(maps$infest3,maps$predictImp3)$bs,
	  "BSobs:",verify(maps$obs3,maps$predictImp3)$bs))
plot_reel(maps$X,maps$Y,maps$predictImp4,base=0,
	  main=paste("Predict forward\nBSinf:",
		     verify(maps$infest4,maps$predictImp4)$bs,
	  "BSobs:",verify(maps$obs4,maps$predictImp4)$bs))

plot_reel(maps$X,maps$Y,maps$infest2,base=0,main="Init perf Obs")
plot_reel(maps$X,maps$Y,maps$bestPredict3,base=0,
	  main=paste("Best Predict training\nBSinf:",
		     verify(maps$infest3,maps$bestPredict3)$bs,
	  "BSobs:",verify(maps$obs3,maps$bestPredict3)$bs))
plot_reel(maps$X,maps$Y,maps$bestPredict4,base=0,
	  main=paste("Predict forward\nBSinf:",
		     verify(maps$infest4,maps$bestPredict4)$bs,
	  "BSobs:",verify(maps$obs4,maps$bestPredict4)$bs))

bsbest<-verify(maps$infest4,maps$bestPredict4)$bs
bsglm<-verify(maps$infest4,maps$predictImp4)$bs

percWorse<- (bsglm-bsbest)/bsbest
cat("% degradation of Brier score by glm compared to perfect",percWorse,"\n")


dev.print(device=pdf,"basic_regression_prediction.pdf")
attributes(maps)$snapTimes<-snapTimes
dump("maps",file="maps_basic_regression.R")

