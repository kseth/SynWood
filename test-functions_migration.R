# tests
source("functions_migration.R")
library("spam")

#========================
# Test that Hop/Skip/Jump matrix is made correctly
#========================
test_that("making hop skip jumps correctly",{
num.rows <- 3 
num.cols <- 3
row.dist <- 1
maps <- makeGrid(num.rows = num.rows, num.cols = num.cols, row.dist = row.dist)
limitHopSkip <- 2 
lowerLimitSkip <- 1 
limitJump <- 3
lowerLimitJump <- 2

stratmat <- generate_stratified_mat(maps, limitHopSkip, limitJump, lowerLimitJump=lowerLimitJump, lowerLimitSkip=lowerLimitSkip)

# all diagonals must be zero
expect_true(sum(diag(stratmat$hopMat)) == 0)
expect_true(sum(diag(stratmat$skipMat)) == 0)
expect_true(sum(diag(stratmat$jumpMat)) == 0)

# corners can hop to 2 spots
expect_true(sum(stratmat$hopMat[1, ]) == 2)
expect_true(sum(stratmat$hopMat[9, ]) == 2)

# corners can skip to 3 spots
expect_true(sum(stratmat$skipMat[1, ]) == 3)
expect_true(sum(stratmat$skipMat[9, ]) == 3)

# corners can jump to 3 spots
expect_true(sum(stratmat$jumpMat[1, ]) == 3)
expect_true(sum(stratmat$jumpMat[9, ]) == 3)

# no hop, skip, jump overlap
totalmat <- stratmat$hopMat + stratmat$skipMat + stratmat$jumpMat
expect_true(sum(diag(totalmat)) == 0)

test_mat <- spam(rep(1, 81), ncol = 9)
diag(test_mat) <- 0
test_mat <- as.spam(test_mat)

expect_equal(totalmat, test_mat)
})

test_that("MakeIndexFromNinfest ok",{
	  set.seed(777)
	  nInfs<-rpois(10,1)
ind<-MakeIndexFromNinfest(nInfs)
expect_equal(length(ind),sum(nInfs))
tabInd <- table(ind)
tabIndFrom_nInfs<-nInfs[as.numeric(names(tabInd))]
attributes(tabInd)<-NULL

expect_equal(tabIndFrom_nInfs,tabInd)

nInfsBack <- MakeNinfestFromIndex(ind,length(nInfs))

expect_equal(nInfsBack,nInfs)
})


#============================
# Test that get_stats_grid works fine with integral input
# Where the input is the number of positive per site (positive houses in block, positive bugs in house, etc.)
#============================
test_that("get_stats_grid returns reasonably for 'poisson' input", {

num.rows <- 33
num.cols <- 33
maps <- makeGrid(num.rows, num.cols, 10)
partitionSizes <- c(10, 20)
map.partitions <- divideMap(maps, partitionSizes, typeDivide = "kmeans")

numHouses <- num.rows*num.cols

infested <- rpois(numHouses, lambda = 5)
maxInfest <- rep(max(infested), numHouses)
infested <- infested - 1
infested[which(infested < 0)] <- 0

out <- get_stats_grid(infested, maxInfest, map.partitions)
expect_true(!any(is.na(out)))

## case where input is binary, not poisson (check to see that it works!)
infested2 <- as.integer(round(runif(numHouses, 0, 1)))
maxInfest2 <- rep(1, numHouses)
out2 <- get_stats_grid(infested2, maxInfest2, map.partitions)
expect_true(!any(is.na(out2)))

## test that only getting 2 partition Lmoments also works
out3 <- get_stats_grid(infested, maxInfest, map.partitions, iPartLMoments = 1:2)
expect_true(!any(is.na(out)))

## case where input is binary, not poisson (check to see that it works!)
## test that only getting 2 partition Lmoments also works
infested2 <- as.integer(round(runif(numHouses, 0, 1)))
maxInfest2 <- rep(1, numHouses)
out4 <- get_stats_grid(infested2, maxInfest2, map.partitions, iPartLMoments = 1:2)
expect_true(!any(is.na(out2)))

})

#============================
# Test that noKernelMultiGilStat works as expected
# For high + low rate move
#============================
# test_that("noKernelMultiGilStat returns names correctly in statsTable",{
source("baseModel.R")
whichPairwise = c("semivariance", "moran", "geary", "ripley")	
implStats <- c("semivariance", "grid", "circles",  "num_inf") # would need to either add "atRisk" or remove it from 
								# noKernelMultiGilStat

# for all names at once
out <- noKernelMultiGilStat(
			    stratHopSkipJump = stratmat, 
			    blockIndex = NULL, 
			    infestH = startingInfested, 
			    timeH=rep(-2), 
			    endTime = 1000, 
			    rateMove = rateMove, 
			    rateHopInMove = 1-rateJumpInMove,
			    rateSkipInMove = 0, 
			    rateJumpInMove = rateJumpInMove, 
			    Nrep = 10, 
			    coords = maps[, c("X", "Y")], 
			    simul=TRUE, 
			    getStats = TRUE, 
			    seed = seed, dist_out = dist_out, 
			    typeStat = implStats, 
			    whichPairwise = whichPairwise,
			    map.partitions = map.partitions, 
			    conc.circs = circles, 
			    rateIntro = 0)
expect_equal(dim(out$statsTable),c(18,10))
allNamesStats <- c()
for(i in 1:length(whichPairwise)){
  allNamesStats<-c(allNamesStats,paste0(whichPairwise[i],1:length(genIntervals)))
}
allNamesStats <- c(allNamesStats,
expect_equal(rownames(out$statsTable),allNamesStats)

# })
test_that("noKernelMultiGilStat num infs calculation",{
	  source("baseModel.R")
out <- noKernelMultiGilStat(
			    stratHopSkipJump = stratmat, 
			    blockIndex = NULL, 
			    infestH = 1, 
			    timeH=rep(-2), 
			    endTime = 1000, 
			    rateMove = rateMove, 
			    rateHopInMove = 1-rateJumpInMove,
			    rateSkipInMove = 0, 
			    rateJumpInMove = rateJumpInMove, 
			    Nrep = 1, 
			    coords = maps[, c("X", "Y")], 
			    simul=TRUE, 
			    getStats = TRUE, 
			    seed = seed, dist_out = NULL, typeStat = c("num_inf"), map.partitions = NULL, conc.circs = NULL, rateIntro = 0)


## all houses should be infested
expect_true(!(any(out$infestedDens == 0)))
expect_equal(out$statsTable, length(maps$X))

## very low movement should fill up the map so that only the initial house is positive
rateMove <- 0.0000000001

out <- noKernelMultiGilStat(stratHopSkipJump = stratmat, blockIndex = NULL, infestH = 1, timeH=rep(-2), endTime = 1, rateMove = rateMove, rateHopInMove=1-rateJumpInMove,rateSkipInMove = 0, rateJumpInMove = rateJumpInMove, Nrep = 1, coords = maps[, c("X", "Y")], simul=TRUE, getStats = TRUE, seed = seed, dist_out = NULL, typeStat = c("num_inf"), map.partitions = NULL, conc.circs = NULL, rateIntro = 0)

## only house 1 is infested
expect_equal(out$statsTable, 1)
expect_equal(sum(out$infestedDens), 1)
expect_equal(which(out$infestedDens==1), 1)
})


### make a basic simulation of dispersal
## parameters
cote<-50
limDist<-3
## grid definition
xs<-rep(seq(1,cote),cote)
ys<-rep(seq(1,cote),each=cote)
dists<-as.matrix(dist(cbind(xs,ys))) # matrix of distances between points
# plot(xs,ys)

basicSimulation<-function(xs,ys,zs,dists){
  ## simulation
  nEnd<-50
  for(i in 1:nEnd){
    # draw an init
    alreadyPos<-which(zs==1)
    if(length(alreadyPos)>1){
      init<-sample(alreadyPos,1)
    }else{
      init<-alreadyPos
    }
    # draw a final
    zs[sample(which(dists[init,]<limDist),1)]<-1
    # plot_reel(xs,ys,zs,base=0)
  }
  pos<-which(zs==1)
  return(list(zs=zs,pos=pos))
}
### test of distance matrix
# test_that("dist matrix in C ok",{ # for some reason expect_error doesn't like it
expect_error(makeDistMat(c(1,2,3),c(1,2)))
distsNoRowNames<-dists
attributes(distsNoRowNames)<-NULL
attributes(distsNoRowNames)$dim<-attributes(dists)$dim
expect_equal(distsNoRowNames,makeDistMat(xs,ys))
# })

### test of at_risk
# test_that("at_risk computation",{
	  trs<-cumsum(1:10)

set.seed(1234)
ncoefsAtRisk<-length(trs)-1
ncolAtRiskStats<-length(trs)+ncoefsAtRisk
atRiskStats<-mat.or.vec(2,ncolAtRiskStats) # the full matrix with multiple stats

# one center
zs<-0*xs
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
out<-basicSimulation(xs,ys,zs,dists)
expect_equal(sum(out$zs),38)
## simple indicator
at_risk<-get_at_risk_indicator(out$pos,dists,trs)
stat<-apply(at_risk,2,sum)
plot(stat)
correct<- c(38,180,382,722,1295,2067,2465,2500,2500,2500)
expect_equal(stat,correct)
## simple stat
at_risk_stat <- get_at_risk_stat(out$pos,dists,trs)
expect_equal(at_risk_stat,correct)

## overall with fit
par(mfrow=c(1,2))
plot(correct)
# C fit
at_risk_fit<-get_stats_at_risk(1,out$pos,dists,trs,atRiskStats,ncoefsAtRisk)
# compute polynom at points
get.predict.at_risk<-function(trs,at_risk_fit){
	result<- 0*trs
	for(iat in 1:length(at_risk_fit)){
		for(itr in 1:(length(trs))){
			# cat("iat:",iat,"coef:",at_risk_fit[iat],"tr:",trs[itr],"result:",result[itr],"\n")
			result[itr] <- result[itr]+at_risk_fit[iat]*trs[itr]^(iat-1)
		}
	}
	return(result)
}
coefsFit<-at_risk_fit[2,(length(trs)+1):(length(trs)+ncoefsAtRisk)]
pred.at_risk<-get.predict.at_risk(trs,coefsFit)

lines(pred.at_risk,col=4)
# R fit
Rfit<-lm(at_risk_stat ~ poly(trs, ncoefsAtRisk-1, raw=TRUE))
lines(predict(Rfit))

# check C vs R fit
plot(Rfit$coefficients,coefsFit)
abline(a=0,b=1)
expect_equal(as.vector(Rfit$coefficients),coefsFit)

# plot the spatial repartition
par(mfcol=c(2,5))
for(i in 1:length(trs)){
	plot(xs,ys,col=at_risk[,i]+2,asp=1,pch=19,cex=0.2)
	points(xs[out$pos],ys[out$pos],pch=3)
}

# three centers
set.seed(1234)
zs<-0*xs
zs[1]<-1
zs[length(zs)]<-1
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
out<-basicSimulation(xs,ys,zs,dists)
expect_equal(sum(out$zs),30)
at_risk3<-get_at_risk_indicator(out$pos,dists,trs)
# par(mfcol=c(3,4))
# for(i in 1:length(trs)){
# plot(xs,ys,col=at_risk3[,i]+2,asp=1,pch=19,cex=0.2)
# points(xs[out$pos],ys[out$pos],pch=3)
# }
lines(apply(at_risk3,2,sum))
stat<-apply(at_risk3,2,sum)
correct<-c(30,132,325,707,1371,1993,2409,2500,2500,2500)
expect_equal(stat,correct)

# })

### test of percolation computations
test_that("percolation computations",{
	  # ashape(xs[pos],ys[pos],alpha=1) # bug, doesn't want to do it

set.seed(1234)
## tests with != starting points 
# basic one cluster
zs<-0*xs
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
# plot_reel(xs,ys,zs,base=0)

out<-basicSimulation(xs,ys,zs,dists)
percGroups<-percolation_circle(dists[out$pos,out$pos],3)
# par(mfrow=c(1,2))
# plot(xs,ys,col=zs+2,type="n")
# text(xs[pos],ys[pos],0:(length(pos)-1))
# plot(xs,ys,pch=".")
# text(xs[pos],ys[pos],percGroups)

expect_equal(percGroups,rep(0,length(out$pos)))

# three clusters
# seed<-sample(10000,1)
# cat("seed:",seed,"\n")
# set.seed(seed)
set.seed(3853)
zs<-0*xs
zs[1]<-1
zs[length(zs)]<-1
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
out<-basicSimulation(xs,ys,zs,dists)
correct<-c(rep(0,6),rep(6,22),rep(28,15))
percGroups<-percolation_circle(dists[out$pos,out$pos],3)
# par(mfrow=c(1,2))
expect_equal(percGroups,correct)
}) # end test_that


#========================
# Old Tests
#========================
## was broken by missing maps, to restore
# test_that("semi-variogram computations",{
# 	  varioTest <- variog(coords = maps[,c("X","Y")], data = maps$infest1, breaks = genIntervals)
# 
# 	  DistClasses<-makeDistClasses(maps$X,maps$Y,genIntervals)
# 	  DistClasses<-makeDistClassesWithStreets(maps$X,maps$Y,genIntervals, blockIndex)
# 	expect_true(all.equal(DistClasses$classSize,varioTest$n))
# 	expect_true(all.equal(DistClasses$classSize, DistClasses$classSizeSB + DistClasses$classSizeAS))
# 
# 	vario1<-variogFromIndices(DistClasses$CClassIndex,maps$infest1,DistClasses$classSize) # variog for init
# 	# print(vario1)
# 	expect_true(all.equal(vario1$variog,varioTest$v))
# 	expect_true(all.equal(vario1$sdvariog,varioTest$sd))
# })
#
#
## was broken by maps not found
# test_that("probmatR and probmatC equivalent",{ 	
# 	threshold <- 2000
# 	### R
# 	dist_mat <- nearest.dist(x=sp, y=NULL, method="euclidian", delta=threshold, upper=NULL);          
# 	dist_mat <- as.matrix(dist_mat)
# 	probMatR <- generate_prob_mat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, threshold, sp, dist_mat, SB, FALSE)
# 	cumulProbMatR <- generate_prob_mat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, threshold, sp, dist_mat, SB, TRUE)
# 
# 	### C
# 	DistClasses<-makeDistClassesWithStreets(maps$X,maps$Y,genIntervals, blockIndex)
# 	probMatC <- generate_prob_mat_C(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, DistClasses$dists, blockIndex, L, FALSE)
# 	cumulProbMatC <- generate_prob_mat_C(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, DistClasses$dists, blockIndex, L, TRUE)
# 
# 	### tests
# 	expect_true(all.equal(cumulProbMatR, cumulProbMatC))
# 	expect_true(all.equal(probMatR, probMatC))
# })
#
# test_that("getPosteriorMaps (and multiGilStats behind the scene)",{
#	source("param.r")
#	source("maps_basic_regression.R")# for maps$X, maps$Y, maps$blockIndex
#	set.seed(1)
#	# and maps$infest3
#	DistClasses<-makeDistClassesWithStreets(maps$X,maps$Y,genIntervals, maps$blockIndex)
#	cumulProbMat <- generate_prob_mat_C(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, DistClasses$dists, blockIndex, cumul=TRUE)

#	thetas<-rnorm(10,mean=rateMove,sd=0.01)
#	thetas<-thetas[which(thetas>0)]
#	test<-getPosteriorMaps(maps,thetas,cumulProbMat=cumulProbMat,maps$infest3,nbit=2*52,repByTheta=10)
#	print(maps$test)
#	expect_true(all.equal(test,maps$test))

	# par(mfrow=c(1,3))
	# plot_reel(maps$X,maps$Y,maps$infest3,base=0,top=1)
	# plot_reel(maps$X,maps$Y,maps$test,base=0,top=1)
	# plot_reel(maps$X,maps$Y,test,base=0,top=1)
	# dump("maps",file="maps_basic_regression.R")
# })

# =============================
# Test of the  DrawFromLinked in C
# =============================

test_that("DrawFromLinked in C ok with simple microInMacro",{
	  set.seed(777)
# source("SynWood/models.R",chdir=TRUE)
X <- seq(0,100)
Y <- seq(200,300)
maps<-as.data.frame(cbind(X,Y))
maps$microInMacro <- 1

stratHopSkipJump <- generate_stratified_mat(coords=maps[, c("X", "Y")], limitHopSkip=10, limitJump=50, lowerLimitJump=10, blockIndex=NULL, lowerLimitSkip=5)

dest <- -1
macroOrigin <- 50

nRep <- 10000
dests<-rep(-1,nRep)

# hops
colInd <- stratHopSkipJump$hopMat@colindices
rowPoint <- stratHopSkipJump$hopMat@rowpointers

for(i in 1:nRep){
  out<-.C("DrawFromLinked",
	  rowPointer = as.integer(rowPoint-1),
	  colPointer = as.integer(colInd-1),
	  macroOrigin = as.integer(macroOrigin-1),
	  dest = as.integer(dest),
	  microInMacro = as.integer(maps$microInMacro)
	  )
  dests[i] <- out$dest+1
}
# everything present
okSet <- colInd[rowPoint[macroOrigin]:(rowPoint[macroOrigin+1]-1)]
expect_true(setequal(dests,okSet)) # may fail every nRep times...

# everything present in right proportion
okN<-maps$microInMacro[okSet]
okProba<-okN/sum(okN)

testN<- table(dests)
testProba <- testN/sum(testN)
expect_true(sum(abs(testProba-okProba))/length(okProba)<100/nRep)
# plot(maps$X,maps$Y)
# with(maps[macroOrigin,],points(X,Y,pch=19,col="blue"))
# with(maps[dests,],points(X,Y,col="yellow"))

# with wild microInMacro
maps$microInMacro <- 1+rpois(dim(maps)[1],lambda=3)
for(i in 1:nRep){
  out<-.C("DrawFromLinked",
	  rowPointer = as.integer(rowPoint-1),
	  colPointer = as.integer(colInd-1),
	  macroOrigin = as.integer(macroOrigin-1),
	  dest = as.integer(dest),
	  microInMacro = as.integer(maps$microInMacro)
	  )
  dests[i] <- out$dest+1
}
# everything present
okSet <- colInd[rowPoint[macroOrigin]:(rowPoint[macroOrigin+1]-1)]
expect_true(setequal(dests,okSet)) # may fail every nRep times...

# everything present in right proportion
okN<-maps$microInMacro[okSet]
okProba<-okN/sum(okN)

testN<- table(dests)
testProba <- testN/sum(testN)
expect_true(sum(abs(testProba-okProba))/length(okProba)<100/nRep)

})

