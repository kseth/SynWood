# Makes a simple plot of all the values in the MCMCs to see 
# the densities

outFolders <- "."

# all? moran's only files for 0.04,0.8
outFolders <- list.files("..",pattern="*_Moran_0.04.*",full.names=TRUE)

outfiles <- list.files(outFolders,pattern="^thetasamples_all.*.txt",full.names=TRUE)

allRuns <- read.table(file = outfiles[1], header = TRUE)
allLengths <- dim(allRuns)[1]

for(num in 1:length(outfiles)){
	allRuns <- rbind(allRuns, read.table(file = outfiles[num], header = TRUE))
	allLengths <- c(allLengths, dim(allRuns)[1])
}

# name the parameters
names(allRuns) <- { if(dim(allRuns)[2] == 4) c("LL", "LP", "rateMove", "rateJump") else c("LL", "LP", "rateMove", "rateJump", "detectRate") }


plot(allRuns$rateMove,allRuns$rateJump,pch=".",xlim=c(0.005,0.075),ylim=c(0,1))
abline(h=0.4)
abline(v=0.04)

