# tests
source("functions_migration.R")

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
ncoefsAtRisk<-5
# one center
zs<-0*xs
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
out<-basicSimulation(xs,ys,zs,dists)
expect_equal(sum(out$zs),38)
## simple stat
at_risk<-at_risk_stat(out$pos,dists,trs)
## overall with fit
atRiskStats<-mat.or.vec(2,length(trs))

at_risk_fit<-get_stats_at_risk(1,out$pos,dists,trs,atRiskStats,ncoefsAtRisk)
expect_equal(at_risk_fit[1,],rep(1,length(trs)+ncoefs))
expect_equal(at_risk_fit[2,],rep(0,length(trs)+ncoefs))


par(mfcol=c(2,5))
for(i in 1:length(trs)){
	plot(xs,ys,col=at_risk[,i]+2,asp=1,pch=19,cex=0.2)
	points(xs[out$pos],ys[out$pos],pch=3)
}

stat<-apply(at_risk,2,sum)
plot(stat)
correct<- c(38,180,382,722,1295,2067,2465,2500,2500,2500)
expect_equal(stat,correct)

# three centers
set.seed(1234)
zs<-0*xs
zs[1]<-1
zs[length(zs)]<-1
zs[which(xs==round(cote/2) & ys == round(cote/2))]<-1
out<-basicSimulation(xs,ys,zs,dists)
expect_equal(sum(out$zs),30)
at_risk3<-at_risk_stat(out$pos,dists,trs)
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

