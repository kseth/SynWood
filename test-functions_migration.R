# tests
source("functions_migration.R")
source("param.r")

test_that("percolation computations",{

### make a basic simulation of dispersal
## parameters
cote<-50
limDist<-3
set.seed(1234)
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

  ## stats
  pos<-which(zs==1)
  # ashape(xs[pos],ys[pos],alpha=1) # bug, doesn't want to do it

  percGroups<-percolation_circle(dists[pos,pos],3)

  # par(mfrow=c(1,2))
  # plot(xs,ys,col=zs+2,type="n")
  # text(xs[pos],ys[pos],0:(length(pos)-1))
  # plot(xs,ys,pch=".")
  # text(xs[pos],ys[pos],percGroups)

  return(list(zs=zs,pos=pos,percGroups=percGroups))
}

## tests with != starting points 
# basic one cluster
zs<-0*xs
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
# plot_reel(xs,ys,zs,base=0)

out<-basicSimulation(xs,ys,zs,dists)
expect_equal(out$percGroups,rep(0,length(out$pos)))

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
expect_equal(out$percGroups,correct)
}) # end test_that

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
# test_that("test probMat for Grant Matrix",{
# 	  X<-seq(0,100,20)
# Y<-seq(0,100,20)
# blockIndex<-c(1,1,2,2,3,3)
# DistClasses<-makeDistClassesWithStreets(X,Y,genIntervals, blockIndex)
# probMat<-generate_prob_mat_Grant_C(42,0.1,10e-4,DistClasses$dists,blockIndex,cumul=FALSE)
# testNotCumul<-matrix(c(
# 0.000,1.000,0.001,0.001,0.001,0.001,
# 1.000,0.000,0.100,0.001,0.001,0.001,
# 0.001,0.100,0.000,1.000,0.001,0.001,
# 0.001,0.001,1.000,0.000,0.100,0.001,
# 0.001,0.001,0.001,0.100,0.000,1.000,
# 0.001,0.001,0.001,0.001,1.000,0.000),nrow=6)
# expect_true(all.equal(testNotCumul,probMat))
# probMat<-generate_prob_mat_Grant_C(42,0.1,10e-4,DistClasses$dists,blockIndex,cumul=TRUE)
# testCumul<-apply(testNotCumul,2,cumsum)
# expect_true(all.equal(testCumul,probMat))
# 
# })

# Time test
# DistClasses<-makeDistClassesWithStreets(maps$X,maps$Y,genIntervals, blockIndex)
# start<-Sys.time()
# probMat<-generate_prob_mat_Grant_C(42,0.1,10e-4,DistClasses$dists,blockIndex,cumul=TRUE)
# cat(Sys.time()-start)

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

