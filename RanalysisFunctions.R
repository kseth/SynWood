###############
# cleaning help functions
###############
### general
# replace NULL by 0 in all factor column of a data frame and change the column to numeric if possible.
# "NULL" can be any word (including NA) and 0 can also be anything else.

set_to<-function(x,init=c("NULL"),final=0){
    # set all in init to final
    # if possible to change the column to numeric do it
    # init can be multiple, final has to be unique

    if(class(x) != 'data.frame')
        stop('x must be a data.frame')
    if(length(final)>1){
        warning("final is of length >1, will only use first item\n")
    }

    isfacts <- sapply(x, is.factor)
    for(colname in names(x)[isfacts]){
        x[,colname] <- factor(x[,colname], levels=unique(c(levels(x[,colname]), final)))
        sel<-which(x[,colname]%in%init)
	if("NA" %in% init){
		cat("NA in init")
		sel<-c(sel,which(is.na(x[,colname])))
	}

        if(length(sel>0)){
            x[sel,colname]<-as.character(final)
        }
        nbNA<-suppressWarnings(length(which(is.na(as.numeric(as.character(x[,colname]))))))
        if(nbNA==0){
            x[,colname]<-as.numeric(as.character(x[,colname]))
        }
    }

    return(x)
}
# return all the lines presenting a somewhere duplicated item
# (duplicated don't return the first of duplicated items)
dup.all<-function(vect){
	linesDuplicates<-which(duplicated(vect)) # every non first version of the duplicate

	if(is.null(dim(vect))){
		dupItems<-unique(vect[linesDuplicates])# duplicated items
	}else{
		dupItems<-unique(vect[linesDuplicates,])# duplicated items
	}
	linesAllDuplicates<-which(vect %in% dupItems)

	return(linesAllDuplicates)
}

# change name of a column using its name
changeNameCol <- function(Table,oldname,newname){
	names(Table)[which(names(Table)==oldname)]<-newname
	return(Table)
}
# remove a column
# Table a data.frame
# colNames a vector of column names
removeCol <- function(Table,colNames){
	return(Table[,names(Table)[-which(names(Table) %in% colNames)]])
}

# unify "colName.x" and "colName.y" in a "colName"
# colName: name of the column without .x/.y extension
# db: database (dataframe)
# unifun: function to apply to the two initial columns to get the new one, default mean
unify_column<-function(colName,db,unifun=mean,
		colNameA=paste(colName,".x",sep=""),
		colNameB=paste(colName,".y",sep="")){
	colNameB=paste(colName,".y",sep="")
	db[[colName]]<-apply(cbind(db[[colNameA]],db[[colNameB]]),1,unifun,na.rm=TRUE)
	db<-db[,-which(names(db) %in% c(colNameA,colNameB))]
	return(db)
}

### specific to project

# Clean P,D,L,V from spaces
cleanPDLV<-function(Table){
	Table$P<-gsub(" ","",Table$P)
	Table$D<-gsub(" ","",Table$D)
	Table$L<-gsub(" ","",Table$L)
	Table$V<-gsub(" ","",Table$V)
	Table$unicode<-paste(Table$P,Table$D,Table$L,Table$V,sep=".")
	return(Table)
}
# clean unicodes from spaces and lower case
cleanUnicodes<-function(unicodes){
	unicodes<-gsub(" ","",unicodes)
	unicodes<-toupper(bugs$UNICODE)
	return(unicodes)
}

# remove the trailing a-Z in unicodes to be able to merge with arequipa_gps_google 
make_unicode_gps<-function(unicodes){
	return(gsub("[A-Z]$","",unicodes,ignore.case=TRUE))
}



################
# convert to utms
################
library("PBSmapping")
LL.to.our.utms<-function(coord){ # columns must be Longitude and Latitude in this order
	coord<-as.data.frame(coord)
	names(coord)<-c("X","Y")
	
	# only try to convert not NAs
	sel<-which(!is.na(coord[,1]) & !is.na(coord[,1]))
	coordDefined<-coord[sel,] 

	attributes(coordDefined)$projection<-"LL"

	utm.coord<-as.data.frame(cbind(rep(NA,dim(coord)[1]),rep(NA,dim(coord)[1])))
	names(utm.coord)<-c("X","Y")
	utm.coord[sel,]<-as.data.frame(convUL(coordDefined))

	our.utms<-utm.coord*1000
	if(max(utm.coord$Y,na.rm=TRUE)<0){
	our.utms$Y<-our.utms$Y+10000000 # CAREFUL may not be needed with new versions of R
	}

	attributes(our.utms)<-attributes(utm.coord)

	return(our.utms)
}
# H may be "S" or "N" if UTMs in South or North hemisphere
our.utms.to.LL<-function(our.utms,zone=FALSE,H=FALSE){ # columns must be X and Y in this order
	names(our.utms)<-c("X","Y")

	if(H=="S" || H=="s"){
		our.utms$Y<-our.utms$Y-10000000# CAREFUL may not be needed with new versions of R
	}else if(H!="N" && H!="n"){
		stop("Missing H (\"S\" or \"N\") for our.utms.to.LL\n");
	}
	our.utms$X<-our.utms$X/1000
	our.utms$Y<-our.utms$Y/1000

	# avoid NAs
	sel<-which(!is.na(our.utms[,1]) & !is.na(our.utms[,1]))
	our.utmsDefined<-our.utms[sel,] 

	attributes(our.utmsDefined)$projection<-"UTM"
	attributes(our.utmsDefined)$zone<-zone

	NAvect<-rep(NA,dim(our.utms)[1])
	back.coord<-as.data.frame(cbind(NAvect,NAvect))
	names(back.coord)<-c("X","Y")
	back.coord[sel,]<-as.data.frame(convUL(our.utmsDefined))
	attributes(back.coord)$projection<-"LL"

	return(back.coord)
}
library(testthat)
## testing / Example
# arequipa with NA
latlongArequipa<-mat.or.vec(3,2)
latlongArequipa[1,]<-c(-71.50535,-16.39649)
latlongArequipa[2,]<-c(-71.50530,-16.39620)
latlongArequipa[3,]<-c(NA,NA)
utmArequipa<-mat.or.vec(3,2)
utmArequipa[1,]<-c(232411.8,8185554)
utmArequipa[2,]<-c(232416.8,8185587)
utmArequipa[3,]<-c(NA,NA)
calcUtmArequipa<-LL.to.our.utms(latlongArequipa)
calcLLArequipa<-our.utms.to.LL(calcUtmArequipa,zone=19,H="S")
expect_true(all(calcUtmArequipa[-3,]-utmArequipa[-3,]<1))
expect_true(all(calcLLArequipa[-3,]-latlongArequipa[-3,]<10e-6))

# philly
latlongPhilly<-mat.or.vec(2,2)
latlongPhilly[1,]<-c(-75.06446,40.03222)
latlongPhilly[2,]<-c(-75.17687,39.98272)
utmPhilly<-mat.or.vec(2,2)
utmPhilly[1,]<-c(494500.2,4431335)
utmPhilly[2,]<-c(484898.5,4425854)
calcUtmPhilly<- LL.to.our.utms(latlongPhilly)
calcLLPhilly<-our.utms.to.LL(calcUtmPhilly,zone=18,H="N")
expect_true(all(calcLLPhilly-latlongPhilly<10e-6))
expect_true(all(calcUtmPhilly-utmPhilly<1))

###############
# analysis helper function
###############
# count the TRUE in a vector
count<-function(boolVect){  
	return(length(which(boolVect)))
}
# return the first element of a vector, useful for aggregate()
first<-function(vect,...){
	return(vect[1])
}
# return the last element of a vector, useful for aggregate()
last<-function(vect,...){
	return(vect[length(vect)])
}
## improved list of objects with size in memory
.ls.objects <- function (pos = 1, pattern, order.by,
			 decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
				       fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
		      as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

# then use lsos() to list the biggest objects in R


##############
# plotting functions
##############
# allow to plot boxplots for custom boundaries of the box and whiskers
# x should be a table with in columns the classes you want to boxplot
boxplot.free<-function(x,breaks=c(0.025,0.25,0.5,0.75,0.975),rm.out=TRUE,...){
	out<-boxplot(x,plot=FALSE)
	out$stats[1,]<-apply(x,2,quantile, probs = breaks[1], na.rm = TRUE, names= FALSE)
	out$stats[2,]<-apply(x,2,quantile, probs = breaks[2], na.rm = TRUE, names= FALSE)
	out$stats[3,]<-apply(x,2,quantile, probs = breaks[3], na.rm = TRUE, names= FALSE)
	out$stats[4,]<-apply(x,2,quantile, probs = breaks[4], na.rm = TRUE, names= FALSE)
	out$stats[5,]<-apply(x,2,quantile, probs = breaks[5], na.rm = TRUE, names= FALSE)
	# cat("outstats:",out$stats,"\n")
	if(rm.out){
		out$out<-c()
		out$group<-c()
	}
	bxp(out,...)
	return(out)
}
# example:
# fakeData<-cbind(rnorm(100,mean=1,sd=1),rnorm(100,mean=3,sd=2))
# boxplot.free(fakeData)

# allow to plot random confidence intervals in a barplot
# y: vector with y of the bars
# yminus: vector with y of the lower end of confidence intervals
# ymax: vector with y of the higher end of confidence intervals
# see barplot() for further commands
library(Hmisc)
barplot.ci<-function(y,yminus,ymax,ylim=c(min(y,yminus,ymax),max(y,yminus,ymax)),...){

	xbar<-barplot(y,ylim=ylim,...)

	errbar(xbar,y,yminus,ymax,add=TRUE,type="n")
	
	return(xbar)
}
# Example
# barplot.ci(1:3,(1:3)-0.5,(1:3)+0.5)

# ### How to plot barplot by groups
# db<-read.csv("DataGroupedBoxplot.csv")
# 
# # get a 2*2 table with groups in columns and subgroups in lines
# data <- tapply(df$total_dist, list(df$groupname,df$bin), sum)
# 
# # plot the 2*2 table
# barplot(data,beside=T,col=c("red","blue")
# ,main="European Parliament Elections",xlab="Group",ylab="Seats")
# 
# legend("topleft",rownames(data),fill=c("#ee7700","#3333ff"))

### equivalent of heat.colors but black to yellow
### like used in plot_reel
yellow.colors<-function(n){
  return(rgb((1:n)/n,(1:n)/n,0))
}

### plot of x,y points with a z scale of color black to yellow
# base: under that coded black
# top: above that coded yellow
plot_reel<-function(x,y,z,base=-1,top=1,asmax=top,asmin=base,zp=FALSE,...){
	# if z has NA remove it
	sel<-which(!is.na(z))
	z<-z[sel]
	y<-y[sel]
	x<-x[sel]
	# z_plot<-((z+1)/2)^(1/2)
	# all same color range
	z[z < base]<- base
	z[z > top]<- top
	# general
	z_plot<-((z-min(z,asmin))/(max(z,asmax)-min(z,asmin)+0.0001))# ^(1/2) # +0.01 avoid NAN
	# cat("min:",min(z),"max z:",max(z));
	# z_plot<-((z-base)/(top-base))#
	# z_plot<-((z+1)/2)^(1/2)

	# ordering to have the biggest on the top of the picture
	x<-x[order(z_plot)]
	y<-y[order(z_plot)]
	z_plot<-z_plot[order(z_plot)]

	red<-z_plot
	green<-z_plot
	blue<-z_plot*0
	zcol<-rgb(red,green,blue)
	if(exists("zoom_loaded") && zp==TRUE){
		zplot(x,y,pch=15,asp=1,type="n",zs=0,...)
		zlines(x,y,pch=15,type="p",col=zcol,cex=0.6,zs=0)
	}else{
		plot(x,y,pch=15,asp=1,type="n",...)
		lines(x,y,pch=15,type="p",col=zcol,cex=0.6)
	}
	#zcol
}

## pie Charts
# install.packages("plotrix")
library(plotrix)
# X,Y X and Y vectors of coordinates of the center of pie charts
# shares: a matrix with in columns the shares for each color
# radius: a vector with the radius of each pie
# col: a vector of the colors to be used for each column of "shares"
# ... additional elements are past to plot() if add=FALSE
floating.pies<-function(X,Y,shares,radius,col,smallOnTop=TRUE,add=TRUE,...){
	# test length for all data
	if(length(unique(length(X),length(Y),dim(shares)[1],length(radius)))!=1){ 
		stop("Not all data vectors have the same length")
	}
	# order by size of radius
	if(smallOnTop){
		newOrder<-order(radius,decreasing=TRUE)
		X<-X[newOrder]
		Y<-Y[newOrder]
		shares<-shares[newOrder,]
		radius<-radius[newOrder]
	}
	# plot background if necessary
	if(!add){
		minx<-min(X)-max(radius)
		miny<-min(Y)-max(radius)
		maxx<-max(X)+max(radius)
		maxy<-max(Y)+max(radius)
		plot(c(minx,maxx),c(miny,maxy),type="n",...)
	}

	# plot the pies
	nbpie<-length(X)
	for(p in 1:nbpie){
		# floating.pie doesn't handle 0 values for the shares
		localShares<-shares[p,]
		notNull<-which(localShares>0)
		localShares<-localShares[notNull]
		localColors<-col[notNull]
		if(length(notNull)<2){
			draw.circle(X[p],Y[p],radius=radius[p],col=localColors)
		}else{
			floating.pie(X[p],Y[p],localShares,radius=radius[p],col=localColors)
		}
	}
}
# # Example
# Xs<-rnorm(20,sd=10)
# Ys<-rnorm(20,sd=10)
# share1<-rpois(20,3)
# share2<-rpois(20,3)
# radius<-sqrt(share1+share2)
# colShares<-c("white","black")
# floating.pies(Xs,Ys,cbind(share1,share2),radius,colShares,add=FALSE,xlab="X",ylab="Y")



# plot shapefiles of streets
library("maptools")
## regenerate the R object file:
# mapLim<-readShapeSpatial("LocalitiesBoundaries/Arequipa_limits")
# mapLim@data
# mapLim@polygons

# # remove the points at the end of the Code if needed
# names(mapLim@data)<-c("Name","Code")
# mapLim@data$Code<-gsub("\\.$","",mapLim@data$Code)
# 
# # split provincia, district and localily level codes
# mapLim@data$P<-gsub("^[0-9]*\\.","",mapLim@data$Code)

### function to plot a given set of localities
# locNames should be a vector of character of the form 1.9.4
# * and other characters for regex can be used
# additional standar arguments to plot can be passed 
# in particular add=TRUE to add a locality on the top of a plot
plot.loc.arequipa<-function(locNames,...){
	# get the number of the polygons
	listNum<-c()
	for(name in locNames){
		listNum<-c(listNum,grep(name,mapLim@data$Code))
	}
	plot(mapLim[listNum,],...)
	return(invisible(listNum))
}
## save this as a separate file
# save(mapLim,plot.loc.arequipa,file="ArequipaLim.R")

## can then be load using simply
# load("ArequipaLim.R")

## then can be used directly
# # example
# par(mfrow=c(1,2))
# plot(mapLim) # to plot all the localities
# locNames<-c("1.8.12","1.8.30")
# plot.loc.arequipa(locNames,add=TRUE,col="blue")
# 
# plot.loc.arequipa("1.8.*")
# plot.loc.arequipa(locNames,add=TRUE,col="blue")
# # and to plot something (like text/data) at the centroid of the polygons you can use
# text(mapLim@polygons[[1]]@labpt[1],mapLim@polygons[[1]]@labpt[2],mapLim@data$Code[[1]])

#### Extrapolate points on a grid
# than can be ploted using persp() or for really nice things persp3d in package rgl (see Transect_functions.R for examples) 
expKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T*exp(-matdist@entries/f);
	}else{
		K<- T*exp(-matdist/f);
	}
	# K<- T*exp(-log(2)*matdist/f); ## normalized to have 0.5 when dist=f
	return(K);
}
gaussianKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T*exp(-(matdist@entries*matdist@entries)/f^2);
	}else{
		K<- T*exp(-(matdist*matdist)/f^2);
	}
	# K<- T*exp(-log(2)*(matdist*matdist)/f^2);## normalized to have 0.5 when dist=f
	return(K);
}
# system.time(for(i in 1:10){A<-gaussianKernel(1,Dmat,20)})
cauchyKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T/(1+(matdist@entries*matdist@entries)/f^2);
	}else{
		K<- T/(1+(matdist*matdist)/f^2);
	}
	return(K);
}
geometricKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T/(1+matdist@entries/f);
	}else{
		K<- T/(1+matdist/f);
	}
	return(K);
}
adjust.lim<-function(limsmall,limbig,steps=NULL,stepsize=NULL){
	# limsmall: c(min,max) for coord with smallest range
	# limbig: c(min,max) for coord with biggest range
  limsmall<-c(min(limsmall),max(limsmall))
  if(is.null(stepsize)){
	steps.small<-steps
	stepsize<-(limsmall[2]-limsmall[1])/steps.small
  }

	limbig<-c(min(limbig),max(limbig))
	size.big.init<-(limbig[2]-limbig[1])
	steps.big<-ceiling(size.big.init/stepsize)
	shift.big<-(steps.big*stepsize-size.big.init)/2
	limbig<-c(limbig[1]-shift.big,limbig[2]+shift.big)

	return(list(limsmall=limsmall,limbig=limbig,stepsize=stepsize))
}
library("spam")
grid.from.sample<-function(known.x,known.y,known.z,
	Kernel=expKernel,f=NULL,T=1,
	xlim=NULL,
	ylim=NULL,
	pixel.size=NULL,
	steps=NULL,
	tr=NULL
	){
	# # tr should always be set
	# if(is.null(tr)){
	# 	tr<-min(max(known.x)-min(known.x),max(known.y)-min(known.y))/(steps/2);
	# }
	if(is.null(xlim)){
		xlim=c(min(known.x)-tr,max(known.x)+tr);
	}
	if(is.null(ylim)){
		ylim=c(min(known.y)-tr,max(known.y)+tr);
	}
	if(is.null(f)){
		f<-tr/4
	}


	# get ToGuess locations
	if(abs(xlim[2]-xlim[1])>abs(ylim[2]-ylim[1])){
	  if(is.null(steps))
	    out<-adjust.lim(ylim,xlim,stepsize=pixel.size)
	  else
	    out<-adjust.lim(ylim,xlim,steps=steps)
		ylim<-out$limsmall
		xlim<-out$limbig
		stepsize<-out$stepsize
	}else{
	  if(is.null(steps))
	    out<-adjust.lim(xlim,ylim,stepsize=pixel.size)
	  else
	    out<-adjust.lim(xlim,ylim,steps=steps)
		xlim<-out$limsmall
		ylim<-out$limbig
		stepsize<-out$stepsize
	}
	cat("xlim:",xlim,"ylim",ylim,"stepsize:",stepsize,"\n")
	# vectors with the xs and ys of the grid
	xs<-seq(xlim[1],xlim[2],stepsize)
	ys<-seq(ylim[1],ylim[2],stepsize)

	# coordinates for all each point
	ToGuess.x<-rep(xs,length(ys))
	ToGuess.y<-as.vector(sapply(ys,rep,length(xs)))

	# get distance matrix
	matdist<-nearest.dist(x=cbind(ToGuess.x,ToGuess.y),y=cbind(known.x,known.y),method="euclidian",delta=tr,upper=NULL)
	
	weightsKnownInToGuessRaw<-Kernel(T,matdist,f); # raw weights

	# get normalized by ToGuess weights
	sumR<-drop(weightsKnownInToGuessRaw%*%rep(1,dim(weightsKnownInToGuessRaw)[2]))
	isolated<-which(sumR<0.01)
	sumRsimple<-sumR
	sumRsimple[isolated]<-1
	NormMat<-diag.spam(1/sumRsimple)
	weightsKnownInToGuess<-NormMat%*%weightsKnownInToGuessRaw

	ToGuess.z<-weightsKnownInToGuess%*%known.z
	ToGuess.z[isolated]<- NA

	return(list(x=ToGuess.x,y=ToGuess.y,z=ToGuess.z,xs=xs,ys=ys,zs=t(matrix(ToGuess.z,nrow=length(ys),byrow=TRUE))
			# ,raw.weights=weightsKnownInToGuessRaw,dists=matdist
			));
}




