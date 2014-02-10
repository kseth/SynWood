num.rows <- 11 
num.cols <- 11
row.dist <- 1
maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)

startingInfested <- round(dim(maps)[1]/2)+1

limitHopSkip <- 2 
lowerLimitSkip <- NULL 
limitJump <-10 
lowerLimitJump <- limitHopSkip 

stratmat <- generate_stratified_mat(maps, limitHopSkip, limitJump, lowerLimitJump=lowerLimitJump, lowerLimitSkip=lowerLimitSkip)

## very high movement should fill up the map so that all houses become filled
rateMove <- 1
rateJumpInMove <-0.5 
seed <- 100000

partitionSizes <- c(5, 8, 12)

# distance classes for the general variogram
genIntervals <- c(0, 3, 6)

# radii for the concentric circles
circleRadii <- genIntervals

# bin the map into different distances classes
dist_out <- makeDistClasses(X = as.vector(maps[, "X"]), 
				Y = as.vector(maps[, "Y"]), 
				genIntervals)

# partition the map
map.partitions <- divideMap(maps, partitionSizes, typeDivide="kmeans", language="R")

# make the concentric circles
circles <- conc.circles(maps$X, maps$Y, circleRadii, startingInfested) 


