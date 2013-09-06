source("../spatcontrol/spatcontrol.R", chdir=TRUE)
library(pracma)
library(geometry)

#========================
# Helper functions
#========================
## helper annotation function
annotateStatCor<-function(statcolsname, ntotstats, values=eval(parse(text=statcolsname)),y=1.07){
	par(xpd=NA)
	shift<- 1/(2*(ntotstats-1))
	x<-min(values-1)/(ntotstats-1)-shift
	lines(c(x,x),c(-shift,1+shift))	
	lines(c(-shift,1+shift),c(x,x))	
	# text(x,y,statcolsname,pos=4)
}

## annotates all the name_stats using the top helper function
## name_stats, names of each of the stats
## ntotstats, the total # of stats (not just length(name_stats) but the sum stats in each class)
## ypos - yposition of annotation (recommended 1.07)
annote.image.stats<-function(name_stats, ntotstats, ypos){
	for(name in name_stats)
		annotateStatCor(name, ntotstats, y=ypos[name])
}

## takes a grid.from.sample object and makes it usable for crI purposes

convertgridfromsample <- function(gridfromsample, setNA=0){

	z <- gridfromsample$zs
	xs <- gridfromsample$xs
	x <- rep(xs, ncol(z)) 
	ys <- gridfromsample$ys
	y <- rep(ys, each=nrow(z))
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
## checks the correlation of the stats against themselves
## example given below
correlation.check <- function(stats, name_stats, ypos, colorsCor=colorRampPalette(c("green","orange","red"))(25), ...){
	ntotstats<-dim(stats)[2]
	correlation<-cor(stats)

	image(abs(correlation), xlab= "correlation", col=colorsCor, xaxt="n", yaxt="n", useRaster=TRUE, ...)
	annote.image.stats(name_stats, ntotstats, ypos)
	return(correlation)
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
partial.correlation.check <- function(stats, name_stats, ypos, vcov.obj=NULL, colorsCor=colorRampPalette(c("green","orange","red"))(25), ...){
	ntotstats<-dim(stats)[2]

	if(is.null(vcov.obj))
		vcov.obj <- robust.vcov(t(stats))

	er <- vcov.obj$E
	out <- get.partial.cor(er)
	
	image(abs(out$partCor), xlab= "partial correlation", col=colorsCor, xaxt="n", yaxt="n", useRaster=TRUE, ...)
	annote.image.stats(name_stats, ntotstats, ypos)

	ret <- list(partCor=out$partCor, vcov.obj=vcov.obj)
	return(ret)
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


## x is the parameter values
## y is the density at the parameter values
## steps is the number of steps at which to bin the area
## prI is the quantiles of precision wanted
## ... passed to the smooth.spline function
oneDim_precI <- function(x, y, steps=length(x)/2, prI=c(0.5, 0.75, 0.95), ...){

	pred_smooth <- smooth.spline(x, y,...)
	px <- seq(min(x), max(x), length.out=steps+1)
        py <- predict(pred_smooth, px)$y

	py[which(py < 0)] <- 0 #set things with negative density to zero density	

	xpairs <- matrix(c(px[1:(length(px)-1)], px[2:length(px)]), ncol = 2)
	ypairs <- matrix(c(py[1:(length(py)-1)], py[2:length(py)]), ncol = 2)

	trap_widths <- apply(xpairs, 1, function(xpair){return(abs(xpair[2]-xpair[1]))})
	trap_heights <- apply(ypairs, 1, mean)
	trap_areas <- trap_widths * trap_heights

	print(sum(trap_areas))
	py <- py / sum(trap_areas) #normalize py
	trap_heights <- trap_heights / sum(trap_areas) #normalize trap_heights
	ypairs <- ypairs/sum(trap_areas) # normalize the py pair values
	trap_areas <- trap_areas / sum(trap_areas) # normalize trap_areas

	plot(px, py, col = "grey", type = "l", xlab = "parameter", ylab = "density")
	points(x, y, pch = ".")

	#order the trap areas in reverse order
	rev_trap_areas <- rev(sort(trap_areas))
	rev_trap_index <- rev(order(trap_areas))
	cum_rev_trap_areas <- cumsum(rev_trap_areas)

	whichInterval <- findInterval(cum_rev_trap_areas, prI, rightmost.closed=TRUE)

	cutoff_ll <- mat.or.vec(length(prI), 1)
	
	for(i in 0:(length(prI)-1)){
		cutoff_ll[i+1] <- which.min((trap_heights[rev_trap_index])[whichInterval == i])
		cutoff_ll[i+1] <- min(ypairs[rev_trap_index, ][whichInterval==i, ][cutoff_ll[i+1], ])
		abline(h=cutoff_ll[i+1], lty = 2)
		text(min(px), cutoff_ll[i+1], label = prI[i+1])
	}

	names(cutoff_ll) <- prI

	out <- list(cutoff_ll=cutoff_ll, px=px, py=py, x=x, y=y)
	return(out)
}

## where x1, x2 are the input variables or planar coordinates (params)
## and y is the output variable (ll) or height coordinate
## x1, x2 should be strictly increasing
## y should be a matrix of lls for dim(x1, x2) - or all combinations of x1 and x2
## prI gives the intervals to be found
## ... are parameters passed to the image display 
## if x1, x2, y are scatterplot data
## 	use grid.from.sample(x1, x2, y, steps=100, tr=1, kern=gaussianKernel, xlim=c(min(x1), max(x1)), ylim=c(min(x2), max(x2)))
##	and pass the out$xs, out$ys, out$zs
twoDim_precI <- function(x1, x2, y, prI=c(0.5, 0.75, 0.95), plotLog=F, ...){

	grid_smooth <- list(xs=x1, ys=x2, zs=y)
	pred_smooth <- convertgridfromsample(grid_smooth)
	px1 <- pred_smooth$x
	px2 <- pred_smooth$y
	py <- pred_smooth$z

	tri <- delaunayn(data.frame(px1=px1, px2=px2)) #triangulate the grid
	volume_out <- trapz3d(px1, px2, py, tri) #find the volume	
	indiv_tri_vols <- attr(volume_out, "tri_volumes") #find the volumes of the individual triangles
	attributes(volume_out) <- NULL #get rid of the attributes

	indiv_tri_vols <- indiv_tri_vols/volume_out #normalize tri volumes
	py <- py/volume_out #normalize lls

	print(volume_out)

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
		print(mch)
		cutoff_index <- match_tri_index[mch]
		print(median(py[tri[cutoff_index, ]]))
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

# given probabilities for a credible interval return 
get_expected_coverage<-function(probs=0.95,likObs,likExp,plot=FALSE){
  ## get the density corresponding to the 5% quantile
  lhTr<-quantile(likExp,1-probs)

  ## is 95% of our likelihoods of stats from the TV above the previous density
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

#=========================
# Fit and get coverage
#=========================
post_cov_analysis <- function(TVstats, NdrawsFromFit=dim(TVstats)[1], probs=c(0.90, 0.95, 0.99), typeFit=c("normal", "skew")){

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
