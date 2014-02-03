source("../spatcontrol/spatcontrol.R", chdir=TRUE)
library(pracma)
library(geometry)

#========================
# Helper functions
#========================
## takes a grid.from.sample object and makes it usable for crI purposes

convertgridfromsample <- function(gridfromsample, setNA=0){

	xs <- gridfromsample$xs
	ys <- gridfromsample$ys
	z <- gridfromsample$zs
	
	x <- rep(xs, length(ys)) 
	y <- rep(ys, each=length(xs))
	z <- as.vector(z)

	bad <- which(is.na(z))
	z[bad] <- setNA

	return(list(x=x, y=y, z=z))
}

## find euclidian distance between two points (x1, y1), (x2, y2)
distance <- function(x1, y1, x2, y2){

	xd <- x2 - x1
	yd <- y2 - y1
	d <- sqrt(xd*xd + yd*yd)
	return(d)
}

## entry - a vector of 3 indexes
## x, y - vectors of x, y positions of all indexes
## finds area of the triangle given by points of entry
triangle_area<-function(entry, x, y){

	if(all(x[entry]==x[entry[1]]) || all(y[entry]==y[entry[1]])){
		warning("flat triangles in data")
		return(0) #flat triangle  case
	}

	# herons formula
	a <- distance(x[entry[1]], y[entry[1]], x[entry[2]], y[entry[2]])
	b <- distance(x[entry[2]], y[entry[2]], x[entry[3]], y[entry[3]])
	d <- distance(x[entry[3]], y[entry[3]], x[entry[1]], y[entry[1]])
	s <- (a+b+d)/2
	area <- sqrt(s*(s-a)*(s-b)*(s-d))

	return(area)
}

## check the skew-normality of a univariate sample 
skewnormal_check<-function(entry){
	cp <- sn.mle(y=entry, plot.it=FALSE)$cp
	dp <- cp.to.dp(cp)
	pval <- ks.test(entry, psn, dp=dp, engine="T.Owen")$p.value
	return(pval)
}

## check the normality of a univariate sample 
normal_check<-function(entry){
	mea<-mean(entry)
	s<-sd(entry)
	pval <- ks.test(entry, pnorm, mean=mea, sd=s)$p.value
	return(pval)
}

#===================================
# computes the synthetic likelihood at all points on the grid of parameter sets
# Nota: use at all points the structure of the SL at the true value
#===================================
# function to get the synthetic likelihood structure from stats at true value 
#! value: list by type of likelihood of lists by likelihood summary stats of matrix 
#     of likelihoods by parameter sets
GetSynLikStructure <- function(trueStats,colNums,typeSL="mvn"){
	# compute the synthetic likelihood structure at the true value
	statsFits <- list()
	if(length(colNums) > 1){
		if("smvn" %in% typeSL){
			# cat("calculating skewnormal params\n")
			statsFits[["smvn"]] <- try(msn.mle(y=trueStats[, colNums])$dp)
		}

		if("mvn" %in% typeSL){
			# cat("calculating normal params\n")
			statsFits[["mvn"]] <- try(robust.vcov(sY=t(trueStats[,colNums])))
		}
	}else{
		if("smvn" %in% typeSL){
			# cat("calculating skewnormal params\n")
			param <- try(sn.mle(y=trueStats[, colNums], plot.it=FALSE)$cp)
			param <- cp.to.dp(param)
			statsFits[["smvn"]] <- param
		}

		if("mvn" %in% typeSL){
			# cat("calculating normal params\n")
			statsFits[["mvn"]] <- c(mean=mean(trueStats[,colNums]), sd=sd(trueStats[,colNums]))
		}
	}
		return(statsFits)
}

# function summarizing a bunch of synthetic likelihoods
VectSummaryStatsOfLLs<-function(lls){
	med <- median(lls, na.rm=TRUE)
	mea <- mean(lls, na.rm=TRUE)
	std <- sd(lls, na.rm=TRUE)
	qua <- quantile(lls, probs=c(0.025, 0.975), na.rm=TRUE)
	rseDefined <- std == 0 ||  mea == 0
	if(!is.na(rseDefined)){
		if(rseDefined){ 
			rse <- 0 
		}else{ 
			rse <- std/mea 
		}
	}else{
		rse<-NA
	}
	sim_ll<-c(med,mea,std,qua,rse)
	return(sim_ll)
}
# lls a matrix, stats in different columns lines are "repeats"
SummaryStatsOfLLs<-function(lls,nCores){
	liks <- exp(lls)
	liksList <- split(liks,row(liks))
	summary_lls<-mclapply(liksList,VectSummaryStatsOfLLs,mc.cores=nCores)
	summary_lls<-t(simplify2array(summary_lls))
	return(summary_lls)
}
## calculate the summary for skew normal synthetic likelihood
LLsStatsSMVN <- function(stats,colNums,statsFit){
	if(length(colNums) > 1){
		lls <- dmsn(stats[, colNums], dp=statsFit,log=TRUE)
	}else{
		lls <- dsn(stats[, colNums], dp=statsFit,log=TRUE)
	}	
	return(lls)
}

## calculate the summary for normal synthetic likelihood
LLsStatsMVN <- function(stats,colNums,statsFit){
	if(length(colNums) > 1){
		lls <- synLik(t(stats[,colNums]),er=statsFit,
				  sY=NULL, trans=NULL)
		lls <- as.numeric(lls) ## remove attributes
	}else{
		lls <- dnorm(as.vector(stats[, colNums]), 
			     mean=statsFit["mean"], 
			     sd=statsFit["sd"],log=TRUE)
	}
	return(lls)
}
### Isolate the parameter set in paramSets closer to given couple
WhichInParamSets<-function(vect,paramSets){
	if(length(vect)==1){
		return(which.min(abs(paramSets-vect)))
	}else{
		expect_equal(length(vect),dim(paramSets)[2])	
		diff <- 0*paramSets
		for(i in 1:length(vect)){
			diff[,i] <- abs(paramSets[,i]-vect[i])
		}
	}
	overAllDist<-apply(diff,1,sum)
	return(which.min(overAllDist))
}
# Tests
a<- cbind(rep(1:5,5),rep(1:5,each=5))
b<-WhichInParamSets(c(2,5),a)
expect_equal(b,22)
b<-WhichInParamSets(c(1,1),a)
expect_equal(b,1)
b<-WhichInParamSets(c(5,5),a)
expect_equal(b,25)


## calculate the summary for normal synthetic likelihood
LLsStatsInRefStats <- function(ref_stats,stats,colNums,typeSL){
	statsFit<-GetSynLikStructure(ref_stats,colNums,typeSL=typeSL)[[typeSL]]
	stats<-rbind(stats[,colNums],ref_stats[,colNums])
	if(class(statsFit)!="try-error"){
		if(typeSL == "mvn"){
			lls<-LLsStatsMVN(stats,1:dim(stats)[2],statsFit)
		}else if(typeSL == "smvn"){
			lls<-LLsStatsSMVN(stats,1:dim(stats)[2],statsFit)
		}
	}else{
		lls<-rep(NA,dim(stats)[1])
	}
	return(lls)
}
# compute the synthetic likelihood 
# at all simulations at trueStats according to otherStats
# only use in otherStats and trueStats the columns colNums (choice of stats)
# typeSL allow to choose the type:
#    ("mvn": multivariate normal; "smvn": skew multivariate normale)
# if trueStats is not given, can specify a parameter set the 
#    that can serve as reference if giving: trueVals and paramSets
# trueVals is a vector of parameters to be used as trueValue
# paramSets a matrix with parameter values for otherStats columns are different
#    parameters and lines are different parameter sets
SynLikExistingStats <- function(trueStats=NULL, otherStats=NULL,
			    colNums=(1:dim(otherStats[[1]])[2]),
			    trueVals=NULL,paramSets=NULL,
			    typeSL="mvn",
			    trueInAll=TRUE,
			    summaryLLs=TRUE){
	nVal <- length(otherStats)

	## manage to always get something as a likelihood structure
	cat("Computing Synthetic Likelihood structure at TV\n")
	if(!is.null(trueStats)){
		statsFits<-GetSynLikStructure(trueStats,colNums,typeSL=typeSL)
	}else{
		expect_equal(nVal,dim(paramSets)[1])
		ok <- FALSE
		keepId<- 1:nVal
		while(!ok | length(keepId)<nVal/10){
			## identify the stats to use as comming from true value
			idTVInKeepId <- WhichInParamSets(trueVals,paramSets[keepId,])
			idTV <- keepId[idTVInKeepId]
			cat("Using (",paramSets[idTV,],") as a proxy for (",trueVals,")\n")
			trueStats <- otherStats[[idTV]]

			statsFits<-GetSynLikStructure(trueStats,colNums,typeSL=typeSL)
			ok<-TRUE
			for(type in typeSL){
				if(class(statsFits[[type]])=="try-error"){
					ok <- FALSE
					keepId <- keepId[-idTVInKeepId]
				}else{
					# check that the fit is not whatever
					lls <- synLik(t(trueStats[,colNums]),er=statsFits[["mvn"]],sY=NULL,trans=NULL)
					maxProba<-exp(max(lls))
					if(!is.finite(maxProba) || maxProba ==0 ){
						ok <- FALSE
						keepId <- keepId[-idTVInKeepId]
					}
				}
			}
		}
	}

	sim_ll<-list()
	for(type in typeSL){
		cat("Computing",type,"Synthetic Likelihood for each point\n")
		sim_ll[[type]]<-list()
		if(trueInAll){
			lls <- t(simplify2array(mclapply(otherStats,LLsStatsInRefStats,
								    trueStats,colNums,type,mc.cores=nCores)))
			nRepTV <- dim(trueStats)[1]
			sim_ll[[type]]$lls <- lls[,(1:nRepTV)]

			nRepOther <- dim(otherStats[[1]])[1]
			sim_ll[[type]]$loc <- lls[,nRepTV+(1:nRepOther)]
		}else{
			if(type == "mvn"){
				LLsFn <- LLsStatsMVN
			}else if(type == "smvn"){
				LLsFn <- LLsStatsSMVN
			}
			sim_ll[[type]]$lls <- t(simplify2array(mclapply(otherStats,LLsFn,
								    colNums,statsFits[[type]],mc.cores=nCores)))
		}
		# make summary statistics of the lls
		if(summaryLLs){
			namesStats <- c("mean","median","sd","loCrI","hiCrI","rse")
			cat("Computing summary stats of the likelihoods\n")
			sim_ll[[type]]$summary <- SummaryStatsOfLLs(sim_ll[[type]]$lls,nCores)
			colnames(sim_ll[[type]]$summary) <- namesStats
			if(!is.null(sim_ll[[type]]$loc)){
				sim_ll[[type]]$summaryLoc <- SummaryStatsOfLLs(sim_ll[[type]]$loc,nCores)
				colnames(sim_ll[[type]]$summaryLoc) <- namesStats
			}
		}
	}
	return(sim_ll)
}
SynLikTrueInAll <- SynLikExistingStats

# compute the synthetic likelihood 
# at all points of otherStats according to trueStats
# only use in otherStats and trueStats the columns colNums (choice of stats)
# typeSL allow to choose the type:
#    ("mvn": multivariate normal; "smvn": skew multivariate normale)
# if trueStats is not given, can specify a parameter set the 
#    that can serve as reference if giving: trueVals and paramSets
# trueVals is a vector of parameters to be used as trueValue
# paramSets a matrix with parameter values for otherStats columns are different
#    parameters and lines are different parameter sets
SynLikAllInTrue <- function(...){
	return(SynLikExistingStats(...,trueInAll=FALSE))
}

#========================
# checking the normality/skew-normality
#========================
stats.norm.check <- function(stats, whichstats, plot=T, alpha=0.05){

	whichstats <- unlist(whichstats) # the names of the stats at each position

	ks_norm <- apply(stats, 2, normal_check)
	ks_skew <- apply(stats, 2, skewnormal_check) ## warning, this step is extremely slow

	if(plot){
		plot(ks_skew, ks_norm, xlab="skew-norm p values", ylab="norm p values", xlim=c(0, 1), ylim=c(0, 1))
		abline(v=alpha, col="red")
		abline(h=alpha, col="red")
		abline(a=0, b=1)
	}

	bad_skew <- which(ks_skew < alpha)
	bad_norm <- which(ks_norm < alpha)
	bad <- union(bad_skew, bad_norm)
	names_bad <- unlist(lapply(bad, function(x, whichstats){
	       				y <- which(whichstats == x)
					y <- paste(names(y), collapse = " ")
					return(y)	}, whichstats=whichstats))

	bad_return <- data.frame(skew=ks_skew[bad], norm=ks_norm[bad], names=names_bad)

	return(list(skew.p=ks_skew, norm.p=ks_norm, sig=bad_return))
}

#========================
# checking the covariance/correlation structure
#========================
## helper annotation function
annotateStatCor<-function(statcolsname, ntotstats, values=eval(parse(text=statcolsname)),y=1.07){
	par(xpd=NA)
	shift<- 1/(2*(ntotstats-1))
	x<-min(values-1)/(ntotstats-1)-shift
	lines(c(x,x),c(-shift,1+shift))	
	lines(c(-shift,1+shift),c(x,x))	
	text(x,y,statcolsname,pos=4,offset=0)
}

## annotates all the name_stats using the top helper function
## name_stats, names of each of the stats
## ntotstats, the total # of stats (not just length(name_stats) but the sum stats in each class)
## ypos - yposition of annotation (recommended 1.07)
annote.image.stats<-function(name_stats,values=NULL, ntotstats, ypos=1.07){
	for(i in 1:length(name_stats)){
		annotateStatCor(name_stats[i], ntotstats,values=values[i], y=max(ypos[i],ypos[1],na.rm=TRUE))
	}
}

## checks the correlation of the stats against themselves
## example given below
correlation.check <- function(stats, name_stats=NULL,firstColstats=NULL, ypos=1.07, colorsCor=colorRampPalette(c("green","orange","red"))(25), ...){
	ntotstats<-dim(stats)[2]
	correlation<-cor(stats)

	image(abs(correlation), xlab= "correlation", col=colorsCor, xaxt="n", yaxt="n", useRaster=TRUE, ...)
	annote.image.stats(name_stats, ntotstats, values=firstColstats,ypos=ypos)
	return(invisible(correlation))
}
# example:
# name_stats <- c("grid_var_stats", "circ_stats", "semivar_newnew_stats", "semivar_oldnew_stats", "moran_stats", "geary_stats", "ripley_stats", "num_inf_stats")
# ypos <- c(1.07, 1.07, 1.03, 1.07, 1.03, 1.1, 1.03, 1.07)
# names(ypos) <- name_stats
# correlation.check(t(real_sims_stats), name_stats, ypos)
# correlation.check(sim_stats, name_stats, ypos)

# er is a result of call to robust.vcov
# returns precision matrix and partialCorrelation of the stats
get.partial.cor <- function(er){
	prec<- er %*% t(er) #make the precision matrix
	invI <- 0*prec
	diag(invI)<-1/sqrt(diag(prec))

	partCor<-  -invI %*% prec%*% invI #make the partial correlation matrix
	diag(partCor)<- -diag(partCor)
	return(list(partCor=partCor, prec=prec))
}

# partial correlation on real value simulations
# example given below
partial.correlation.check <- function(stats, name_stats=NULL,firstColstats=NULL, ypos=1.07, vcov.obj=NULL, colorsCor=colorRampPalette(c("green","orange","red"))(25), ...){
	ntotstats<-dim(stats)[2]

	if(is.null(vcov.obj)){
		nObs<-dim(stats)[1]
		nStats<-dim(stats)[2]
		if(nObs<nStats){
			warning("Need at least",nStats,"observations to compute the partial correlation\n")
			return(NULL)
		}else{
			vcov.obj <- robust.vcov(t(stats))
		}
	}

	er <- vcov.obj$E
	out <- get.partial.cor(er)
	
	image(abs(out$partCor), xlab= "partial correlation", col=colorsCor, xaxt="n", yaxt="n", useRaster=TRUE, ...)
	annote.image.stats(name_stats, ntotstats, values=firstColstats,ypos=ypos)

	ret <- list(partCor=out$partCor, vcov.obj=vcov.obj)
	return(invisible(ret))
}	
# example:
# name_stats <- c("grid_var_stats", "circ_stats", "semivar_newnew_stats", "semivar_oldnew_stats", "moran_stats", "geary_stats", "ripley_stats", "num_inf_stats")
# ypos <- c(1.07, 1.07, 1.03, 1.07, 1.03, 1.1, 1.03, 1.07)
# names(ypos) <- name_stats
# partial.correlation.check(t(real_sims_stats), name_stats, ypos)
# or 
# partial.correlation.check(t(real_sims_stats), name_stats, ypos, vcov.obj=robust.vcov(real_sims_stats))
# also
# partial.correlation.check(sim_stats, name_stats, ypos)
# plot all four correlations at once:
# dev.new()
# par(mfrow=c(2,2), mar=c(1,1,1,1))
# correlation.check(t(real_sims_stats), name_stats, ypos)
# correlation.check((sim_stats), name_stats, ypos)
# partial.correlation.check(t(real_sims_stats), name_stats, ypos)
# partial.correlation.check((sim_stats), name_stats, ypos)
# dev.print(dev=pdf,"correlation_partcor.pdf")

#=====================
# Functions to compute credible intervals by manipulating the likelihood profile 
#=====================
## finding a pseudo-trapezoidal volume in 3space
## x, y, z - coordinates of points in 3 space
## tri - the triangulation to be used
##	the format of the tri object should be the format of the delaunay object
##	a matrix of (number of triangle) by 3 giving the 3 points in each triangle
trapz3d <- function(x, y, z, tri=NULL){

	if(is.null(tri))
		tri <- delaunayn(data.frame(x=x, y=y))

	individual_areas <- apply(tri, 1,
			  		  triangle_area, x=x, y=y)

	heights <- apply(tri, 1,
				 function(entry, heights){
					a <- heights[entry[1]]
					b <- heights[entry[2]]
					c <- heights[entry[3]]
					mean <- (a+b+c)/3
					return(mean)
				 			}, heights=z)

	tri_volumes <- heights*individual_areas
	volume <- sum(tri_volumes)
	attr(volume,"tri_volumes") <- tri_volumes
	return(volume)						  
}

## where x1, x2 are the input variables or planar coordinates (params)
## and y is the output variable (ll) or height coordinate
## xy should contain on each line parameters sets corresponding to the entries of y
## x1, x2 should be strictly increasing (if xy not defined)
## y should be a matrix of lls for dim(x1, x2) - or all combinations of x1 and x2
## prI gives the intervals to be found
## ... are parameters passed to the image display 
## if x1, x2, y are scatterplot data
## 	use grid.from.sample(x1, x2, y, steps=100, tr=1, kern=gaussianKernel, xlim=c(min(x1), max(x1)), ylim=c(min(x2), max(x2)))
##	and pass the out$xs, out$ys, out$zs
twoDim_precI <- function(xy=NULL,x1=NULL, x2=NULL, y, prI=c(0.95), plotLog=F, col=NULL,...){
	if(is.null(xy)){
		grid_smooth <- list(xs=x1, ys=x2, zs=y)
		pred_smooth <- convertgridfromsample(grid_smooth)
		px1 <- pred_smooth$x
		px2 <- pred_smooth$y
		py <- pred_smooth$z
	}else{
		px1 <- xy[,1]
		px2 <- xy[,2]
		x1 <- sort(unique(px1))
		x2 <- sort(unique(px2))
		py <- as.vector(y)
		py[which(is.na(py))] <- 0
	}
	if(is.null(col)){
		col <- colorRampPalette(c("red", "orange", 
					  "yellow", 
					  "green"))(length(x1))
	}


	# TODO: can handle any xy and not just on a grid

	tri <- delaunayn(data.frame(px1=px1, px2=px2)) #triangulate the grid
	volume_out <- trapz3d(px1, px2, py, tri) #find the volume	
	indiv_tri_vols <- attr(volume_out, "tri_volumes") #find the volumes of the individual triangles
	attributes(volume_out) <- NULL #get rid of the attributes

	indiv_tri_vols <- indiv_tri_vols/volume_out #normalize tri volumes
	py <- py/volume_out #normalize lls

	# print(volume_out)

	zmat <- t(matrix(py, nrow=length(x2), byrow=TRUE))
	if(!plotLog)
		image(x=x1, y=x2, z=zmat, useRaster=TRUE, ...)
	else
		image(x=x1, y=x2, z=log(zmat), useRaster=TRUE, ...)

	#order the tri_volumes in reverse order
	rev_tri_volumes <- rev(sort(indiv_tri_vols))
	match_tri_index <- match(rev_tri_volumes, indiv_tri_vols)
	cum_rev_tri_volumes <- cumsum(rev_tri_volumes)

	whichInterval <- findInterval(cum_rev_tri_volumes, prI, rightmost.closed=TRUE)

	cutoff_ll <- mat.or.vec(length(prI), 1)
	
	for(i in 0:(length(prI)-1)){
		mch <- max(which(whichInterval == i))
		# print(mch)
		cutoff_index <- match_tri_index[mch]
		# print(median(py[tri[cutoff_index, ]]))
		cutoff_ll[i+1] <- median(py[tri[cutoff_index, ]])
	}

	contour(x=x1, y=x2, z=zmat, levels=cutoff_ll[!is.na(cutoff_ll)], labels=prI, add=TRUE, labcex=0.8)

	names(cutoff_ll) <- prI

	out <- list(cutoff_ll=cutoff_ll, gx=x1, gy=x2, gz=zmat)
	return(out)
}

#========================
# Old functions left in for continuity
#========================
## x is the x coordinates
## y is the density at the x coordinates
## all y >= 0
## 1-alpha gives the precision inteval
## sigfig is number of significant figures to evaluate to
## returns only the cutoff ll
density_crI <- function(x, y, alpha=0.05, sigfig=5){

	#sort (x, y) such that <x> is in order
	orderx <- order(x)
	x <- x[orderx]
	y <- y[orderx]

	alpha <- signif(alpha, sigfig)

	hi <- max(y)
	lo <- min(y)
	mid <- (hi+lo)/2

	ycheck <- y
	ycheck[y < mid] <- 0

	area <- signif(trapz(x, ycheck), sigfig)
	oldarea <- area
	newarea <- 1-alpha

	if(1-area == alpha)
		return(list(ll=mid))
	if(1-area < alpha){
		mid <- c(mid, (mid+hi)/2)
		lo <- mid[length(mid)-1]
	}else{
		mid <- c(mid, (mid+lo)/2)
		hi <- mid[length(mid)-1]
	}
	while(lo < hi &&  oldarea!=newarea){
	
		oldarea <- area
		ycheck <- y
		ycheck[y < mid[length(mid)]] <- 0
		area <- signif(trapz(x, ycheck), sigfig)
		newarea <- area

		if(1-area == alpha)
			return(list(ll=mid[length(mid)]))
		if(1-area < alpha){
			mid <- c(mid, (mid[length(mid)]+hi)/2)
			lo <- mid[length(mid)-1]
		}else{
			mid <- c(mid, (mid[length(mid)]+lo)/2)
			hi <- mid[length(mid)-1]
		}
	}	

	return(list(ll=mid[length(mid)-1]))
}

## where x1, x2 are the input variables or planar coordinates (params)
## and y is the output variable (ll) or height coordinate
## y >= 0
## 1-alpha gives the precision inteval
## sigfig is number of significant figures to evaluate to
## triangulation of x1, x2 as specified in trapz3d - defaults to NULL and delaunay triangulation
## returns only the ll cut off
volume_crI <- function(x1, x2, y, alpha=0.05, sigfig=5, tri=NULL){
	
	## triangulate to pass to volume function
	if(is.null(tri))
		tri <- delaunayn(data.frame(x1=x1, x2=x2))

	alpha <- signif(alpha, sigfig)

	hi <- max(y)
	lo <- 0
	mid <- (hi+lo)/2

	ycheck <- y
	ycheck[y < mid] <- 0

	vol <- signif(trapz3d(x1, x2, ycheck, tri=tri), sigfig)
	oldvol <- vol
	newvol <- 1-alpha
	print(vol)
	if(isTRUE(all.equal(1-vol, alpha)))
		return(list(ll=mid))
	if(1-vol < alpha){
		mid <- c(mid, (mid+hi)/2)
		lo <- mid[length(mid)-1]
	}else{
		mid <- c(mid, (mid+lo)/2)
		hi <- mid[length(mid)-1]
	}
	while(lo < hi && !isTRUE(all.equal(oldvol, newvol))){
	
		oldvol <- vol
		ycheck <- y
		ycheck[y < mid[length(mid)]] <- 0
		vol <- signif(trapz3d(x1, x2, ycheck, tri=tri), sigfig)
		newvol <- vol
		print(vol)
		if(isTRUE(all.equal(1-vol, alpha)))
			return(list(ll=mid[length(mid)]))
		if(1-vol < alpha){
			mid <- c(mid, (mid[length(mid)]+hi)/2)
			lo <- mid[length(mid)-1]
		}else{
			mid <- c(mid, (mid[length(mid)]+lo)/2)
			hi <- mid[length(mid)-1]
		}
	}	

	return(list(ll=mid[length(mid)-1]))
}

#=========================
# Fit and get coverage
#=========================
# given probabilities for a credible interval return 
get_expected_coverage<-function(probs=0.95,likObs,likExp,plot=FALSE){
  # probs: 1-alpha of the coverage beeing tested
  # likObs: likelihood of the observed data
  # likExp: likelihood of the data generated by the distribution fitted to the observation 
  # plot: should we plot results?

  ## get the density corresponding to the probs quantile
  lhTr<-quantile(likExp,1-probs)

  ## is probs of our likelihoods of stats from the TV above the previous density
  get_cov<-function(lhTr,likObs){1-count(likObs<lhTr)/length(likObs)}
  cov<-sapply(lhTr,get_cov,likObs)

  # cat("cov is",cov,"in spite of",probs,"\n")
  # cat("or out",(1-cov)/(1-probs),"too much\n")

  if(plot){
    ## plotting
    par(mfrow=c(1,2))
    hist(likObs,breaks=100)
    abline(v=lhTr)
    hist(likExp,breaks=100,xlim=range(likObs))
    abline(v=lhTr)
  }

  return(cov)
}

#' Estimation of coverage to be expected for a given set of parameters 
#'
#' @description From a matrix of statistics drawn from at a 
#' given set of parameters returns an expected. 
#' This is similar to cross-validation in the ABC world.
#' @param TVstats a matrix with the value of the statistics 
#'        obtained with a number of draws from the expected value
#'        columns are expected to correspond to different statistiques
#'        lines are correspond to the different repetitions
#' @param NdrawsFromFit the number of draws that should be taken from the
#'        fitted multivariate distribution
#' @param probs a vector of credible interval probability sizes
#' @param typefit a vector of the multivariate distributions to be tried
#'        the ones implemented are:
#'	  \itemize{
#'          \item{"normal"} for multivariate skew
#'          \item{"skew"} for multivariate skew normal
#'	  }
#' @return a data frame with in lines the different probabilities and 
#'       in columns the expected coverages according to each fitting 
#'       multivariate distribution
#'
#' @details 
#'        The indication of coverage should be understood 
#'        as a multivariate 
#'        coverage: the probability that the set of parameter used to 
#'        generate the statistics would be in the joint credible region
#'        when fitting this model with the synthetic likelihood.
#'        In practice, this function fits a multivariate distribution to the 
#'        statistics, calculate the likelihoods of these statistics 
#'        according to the fitted distribution, draw a number 
#'        of values from the fitted distribution, calculate the
#'        likelihoods of these drawn values and finally compare the 
#'        distributions of the likelihoods of the stats and the
#'        likelihoods of the drawn values.
#'        In particular, the multivariate distribution of the statistics 
#'        should have 5% of its likelihood under the 5% lower quantile of 
#'        the likelihoods drawn from the distribution. 
#' @author Corentin M. Barbu
#' @seealso zm(), session.zoom().
#' @examples 
#' require(mvtnorm)
#' 
#' # generate multivariate normal statistics
#' sigma <- matrix(c(4,2,2,3), ncol=2) # covariance matrix
#' stats <- rmvnorm(n=1000, mean=c(1,2), sigma=sigma) # statistics generated by the model
#' local.coverage(stats) # without surprise the normal and skew normal perform equally well
#'
#' @export local.coverage

local.coverage <- function(TVstats, NdrawsFromFit=dim(TVstats)[1], probs=c(0.90, 0.95, 0.99), typeFit=c("normal", "skew")){

	if(is.vector(TVstats)){
		warning("only 1 statistic passed")
		return()
	}

	cat("Fit and get coverage...\n")
	### with multivariate normality

	probsObj <- probs

	if("normal" %in% typeFit){
		cat("Multivariate normal...\n")
		# get structure and the likelihoods for the TV
		sl<-synLik(t(TVstats),t(TVstats),  trans = NULL)
		er<-attr(sl,"er")
		Q<-t(er$E)%*% er$E # precision matrix
		attributes(sl)<-NULL
		sls<-sl # likelihoods without all the mess

		# draw stats from the fitted field
		sampMVF<-rmvnorm.prec(NdrawsFromFit,mu=er$mY,Q=Q)

		# get corresponding likelihoods
		sl<-synLik(t(TVstats),t(sampMVF),  trans = NULL,er=er)
		attributes(sl)<-NULL
		slMVF<-sl

		###### get joint coverage 
		covNorm<-get_expected_coverage(probs,sls,slMVF)
		probsObj <- cbind(probsObj,covNorm)
	}

	##### with skew normality
	if("skew" %in% typeFit){
		cat("Multivariate skew normal...\n")
		# fit
		require("sn")
		fit.dp <- msn.mle(y=TVstats)$dp
		# get likelihood of observed in fit
		ssnls <- dmsn(x=TVstats, dp=fit.dp, log=T)

		# Draw 10000 from the fitted multivariate space
		sampMVF <- rmsn(NdrawsFromFit, dp=fit.dp)

		# likelihood of drawn
		slSMVF <- dmsn(x=sampMVF, dp=fit.dp, log=T)

		###### get joint coverage 
		covSNorm<-get_expected_coverage(probs,ssnls,slSMVF)
		probsObj <- cbind(probsObj,covSNorm)
	}


	return(probsObj)
}

#===========================
# For each probability alpha passed into alpha (alpha = 0.05 corresponds to 95% credible intervals)
# 	Return the coverage of the MCMCs given in allRuns
#	allRuns is a vector of one parameter from all MCMCs concatenated together
#	allLengths is the length of each MCMC
#	realMean is the actual value that we want to be covered
#===========================

cred_cov_anal <- function(allRuns, allLengths, realMean, alpha=c(0.05, 0.10, 0.20)){

	probs_left <- c(alpha/2) # left quantiles
	probs_right <- c(1-alpha/2) # right quantiles
	
	quantile_left <- quantile(allRuns[1:allLengths[1]], probs = probs_left)
	quantile_right <- quantile(allRuns[1:allLengths[1]], probs = probs_right)
	
	#using the counts, determine how many times we fall outside of the confidence interval
	counts <- mat.or.vec(length(alpha), 1)
	names(counts) <- alpha
	
	for(index in 1:length(alpha)){
			if(realMean >= quantile_left[index] && realMean <= quantile_right[index])
				counts[index] <- counts[index]+1
		}

	#calculate the measures over the rest of the runs
	for(run in 2:length(allLengths)){

		quantile_left <- quantile(allRuns[(allLengths[run-1]+1):allLengths[run]], probs = probs_left)
		quantile_right <- quantile(allRuns[(allLengths[run-1]+1):allLengths[run]], probs = probs_right)

		for(index in 1:length(alpha)){
			if(realMean >= quantile_left[index] && realMean <= quantile_right[index])
				counts[index] <- counts[index]+1
		}
	}

	coverage <- counts/length(allLengths)
	return(data.frame(alpha=alpha,credibility=1-alpha,coverage=coverage))	
}

library("synlik") # as should be installed from https://bitbucket.org/cbarbu/synlik
