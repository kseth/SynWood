library(spam);
library(sp)
library(msm)
library("truncnorm")
library("PBSmapping")
library(testthat)
library(Hmisc)
library(plotrix)
library("maptools")
library(LaplacesDemon)
library(locfit)
library("binom")
library(fields)
#===============================
# General purpose functions
#===============================
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


## return the matrix A extended to the given dimensions
## keeping the values in it
resized<-function(A,nr=nrow(A),nc=ncol(A)){
	B<-as.matrix(mat.or.vec(nr,nc));
	B[1:(dim(A)[1]),1:dim(A)[2]]<-A
	return(B);
}

## improved list of objects with size in memory
## use lsos() to list the biggest objects in R
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
#-----------------
# convert to utms
#-----------------
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


#--------------------------------------
# Functions extending spam capabilities
#--------------------------------------
# sum on rows of a matrix
rowsum.spam <- function(A){
	return(A%*%rep(1,A@dimension[2]));
}
# generic function that send back to the rowsum system
rsum<-function(x,...)
UseMethod("rowsum")

fn_at_not_null.spam<-function(A,num_row,funct){
	return(funct(A@entries[A@rowpointers[num_row]:A@rowpointers[num_row+1]]));
}
apply_by_row_not_null.spam<-function(A,funct,void.as=NA,...){
	# extremely much faster than apply.spam
	# careful, when row has no value this logically returns NA
	# unless void.as.NA is set to something different than NA
	results<-mat.or.vec(A@dimension[1],1);
	
	for(num_row in 1:A@dimension[1]){
		if(A@rowpointers[num_row]!=(A@rowpointers[num_row+1])){
			results[num_row]<-funct(A@entries[A@rowpointers[num_row]:(A@rowpointers[num_row+1]-1)],...);
		}else{
			results[num_row]<-NA
		}
	}
	if(!is.na(void.as)){
		results[is.na(results)]<-void.as
	}
	
	return(results);
}

importOk<-try(dyn.load("spatcontrol.so"),silent=TRUE)
if(class(importOk)!="try-error"){
	random_spam_entries_by_row<-function(A){
		# randomly exchange the *defined entries* of a spam matrix by rows
		# allows to permut the neighbors
		out<-.C("resample_spam_entries_by_row",
			entries=as.numeric(A@entries),
			pnb_ent=as.integer(length(A@entries)),
			rpoint=as.integer(A@rowpointers),
			pnb_rows=as.integer(A@dimension[1]))
		B<-A
		B@entries<-out$entries

		return(B);
	}
}else{
	random_spam_entries_by_row<-function(A){
		# randomly exchange the *defined entries* of a spam matrix by rows
		for(r in 1:A@dimension[1]){ # for all rows
			select<-(A@rowpointers[r]:(A@rowpointers[r+1]-1))
			nb<-length(select)
			A@entries[select]<-sample(A@entries[select]);
		}
		return(A);
	}
}

# multiply by f the entries of mat that correspond to an entry in indexMat
if(class(importOk)!="try-error"){
SpecificMultiply.spam<-function(f,mat,indexMat){
	out <- .C("specific_multiply",
		  nbent_dm = as.integer(mat@dimension[1]),
		  entries_dm = as.numeric(mat@entries),
		  nbent_dm = as.integer(length(mat@entries)),
		  rpoint_dm = as.integer(mat@rowpointers),
		  colind_dm = as.integer(mat@colindices),
		  entries_indexMat = as.numeric(indexMat@entries),
		  nbent_indexMat = as.integer(length(indexMat@entries)),
		  rpoint_indexMat = as.integer(indexMat@rowpointers),
		  colind_indexMat = as.integer(indexMat@colindices),
		  f = as.numeric(f))
	Dmat<-mat;
	Dmat@entries<-out$entries_dm;

	return(Dmat);
}
}else{
	SpecificMultiply.spam<-function(f,mat,indexMat){
		return(mat+indexMat*mat*(f-1))
	}
}

# change all entries in A bigger than maxtobesetnull to 0
if(class(importOk)!="try-error"){
filter.spam<-function(A,maxtobesetnull){
	out <- .C("filter_spam",
		  entries = as.numeric(A@entries),
		  nbent = as.integer(length(A@entries)),
		  limit = as.numeric(maxtobesetnull))
	A@entries<-out$entries
	return(as.spam(A));
}
}else{
	cat("\nWARNING filter_spam.so not available, will use pure R code, may be slower.\n")
	filter.spam<-function(A,maxtobesetnull){
		A@entries[A@entries>maxtobesetnull]<- 0
		return(A);
	}
}
symmetric.random.entries.by.line.spam.c<-function(A){
	out <- .C("symmetric_resample_spam_entries_by_row",
		  entries = as.numeric(A@entries),
		  coli = as.integer(A@colindices),
		  nbent = as.integer(length(A@entries)),
		  rpoint = as.integer(A@rowpointers),
		  nbrows = as.integer(dim(A)[1]))
	A@entries<-out$entries
	return(A);
}
# pseudo log likelihood for a multivariate normal vector r given it's center b and it's precision matrix Q
dmvnorm.canonical<-function(r,b,Q,logout=TRUE,cholQ=NULL,...){
	if(is.null(cholQ)){
		cholQ <- chol.spam(Q,...);
	}
	# log.det<- -2*determinant(cholQ)$modulus;
	log.det<- determinant(cholQ)$modulus*2
	mu_r <- drop(solve.spam(cholQ, b)) 
	# print("mu_r:")
	# print(mu_r)
	X<-(r-mu_r);
	# print("X:")
	# print(X)
	# log.det<-sum(log(diag(cholQ)))
	in_exp<-drop(t(X)%*%Q%*%X);

	d <- ncol(Q)
	logpi<-d*log(2*pi);

	cat("in_exp:",in_exp,"log.det",log.det,"logpi",logpi,"\n");
	lnP <- -1/2 * (in_exp - log.det+ logpi);

	if (logout==TRUE){
		return(lnP);
	}else{
		return(exp(lnP));
	}
}

dmvnorm.prec<-function(r,mu,Q,logout=TRUE,...){
	# recicling chol and tweaking it (see help(chol)) could still save a lot of time
	cholQ <- chol.spam(Q,...);
	# log.det<- -2*determinant(cholQ)$modulus;
	# log.det<- determinant.spam(Q,Rstruct=cholQ)$modulus
	log.det<- determinant(cholQ)$modulus*2
	X<-(r-mu);
	# log.det<-sum(log(diag(cholQ)))
	in_exp<-drop(t(X)%*%Q%*%X);

	d <- ncol(Q)
	logpi<- d*log(2*pi);

	# cat("in_exp:",in_exp,"log.det",log.det,"logpi",logpi,"\n");
	lnP <- -1/2 * (in_exp - log.det+ logpi);

	if (logout==TRUE){
		return(lnP);
	}else{
		return(exp(lnP));
	}
}
pseudo_inv<-function(A){
	A.svd<-svd(A);
	Ddiag<-as.spam(diag(A.svd$d));
	A.plus<-A.svd$v%*%(1/Ddiag)%*%t(A.svd$u)
	return(A.plus);
}


#--------------------------------------
# Plotting
#--------------------------------------
strongColors<-c("black","red","green3","blue","skyblue","magenta","yellow","purple","yellow","grey","orange","slategrey","navyblue","darkgreen")

plot.palette<-function(colVect=palette()){
  N<-length(colVect)
  plot(1:N,1:N,type="n")
  abline(v=1:N,col=colVect,lty=1)
}
# Ex: plot.palette(strongColors)
# 
#     or equivalently
# 
#     palette(strongColors)
#     plot.palette()

class.colors<-function(C,colVect=colors()){
  # return a vector of colors corresponding to the values of C understood
  # as classes 
  # allow call in plot(...,col=class.colors(C))
  Cc<-as.factor(C)
  couleurs<-sample(colVect,length(levels(Cc)),replace=TRUE)
  return(couleurs[Cc])
}
plot.classes<-function(X,Y=NULL,C,asp=1,pch=15,...){
  # plot (X,Y) point with colors according to their C classes
  plot(X,Y,col=class.colors(C),pch=pch,asp=asp,...)
}
# plot.classes(db$X,db$Y,db$GroupNum)


# barplot with confidence intervals
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

# allow to plot boxplots for custom quantiles for the box and whiskers
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
	z_plot<-((z-min(z,asmin))/(max(z,asmax)-min(z,asmin)+0.01))# ^(1/2) # +0.01 avoid NAN
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
		zplot(x,y,pch=15,asp=1,zs=0,...)
		zlines(x,y,pch=15,type="p",col=zcol,cex=0.6,zs=0)
	}else{
		plot(x,y,pch=15,asp=1,...)
		lines(x,y,pch=15,type="p",col=zcol,cex=0.6)
	}
	#zcol
}
make.col.persp<-function(z,nbcol=100,color.function=jet.colors){
	# make col vector for persp from base package to get colors according to z
	# nrz: number of columns in z
	nrz <- nrow(z)
	ncz <- ncol(z)

	nbcol <- 100
	color <- color.function(nbcol)
	# Compute the z-value at the facet centers
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)
	zcol<-color[facetcol]
	return(zcol)
}



printdev<-function(...){
	# if(!exists("X11_plot")){
	if(interactive()){
		out<-try(X11())
		if(class(out)=="try-error"){
			X11_plot<<-FALSE
		}else{
			dev.off();
			X11_plot<<-TRUE
		}
		# }
	}else{
		X11_plot<<-FALSE
	}

	if(X11_plot){
		dev.print(...);
	}
	return(X11_plot);
}

#### Extrapolate points on a grid
# than can be ploted using persp() or for really nice things persp3d in package rgl (see Transect_functions.R for examples) 

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
grid.from.sample<-function(known.x,known.y,known.z,
	kern=expKernel,f=NULL,T=1,
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
	
	weightsKnownInToGuessRaw<-kern(T,matdist,f); # raw weights

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



#===============================
# Functions specific to the autocorrelation analysis
#===============================
# moran's I on spam matrices
moran.spam<-function(dmt,values){
	# dmt a spam matrix with the weights between points
	# values the values for each location

	# correction of dmt for correct calculus of I
	dim_mp<-length(values)
	one_vect<-rep(1,dim_mp)
	#need 0 on the diagonal of weight
	diag(dmt)<-0
	#need weight normalized by row
	rsum<-dmt%*%one_vect
	rsum[rsum==0]<-1;
	rsum_mat<-dmt*0;
	diag(rsum_mat)<-1/rsum
	dmt<-rsum_mat%*%dmt

	# matrix of residuals
	res_z<- values-mean(values)
	res_z_mat<-diag.spam(1,length(values));
	diag.spam(res_z_mat)<-res_z

	### Moran's I
	# covariance term
	int<-dmt%*%res_z_mat;
	int<-res_z_mat%*%int

	MI<-length(res_z_mat)*one_vect%*%int%*%one_vect/(sum(dmt@entries)*res_z%*%res_z)
	return(drop(MI))
}
gen.mats.neigh<-function(distances,x,y,group=NULL){
	dist_matA <-nearest.dist(x=cbind(x,y), y=NULL, method="euclidian", delta=max(distances), upper=NULL);          

	dmtA<-dist_matA
	dimension<-dim(dist_matA)[1]
	dmtA@entries<-rep(1,length(dmtA@entries))# [dmtA@entries!=0]<-1 # 1 only when dist_matA not 0
	if(!is.null(group)){
		SBA <- nearest.dist(x=cbind(group,rep(0,length(group))), method="euclidian", upper=NULL,delta=0.1)
		SBA@entries<-rep(1,length(SBA@entries))
		SBA@entries<-rep(1,length(SBA@entries))
		SBA<-SBA*dmtA;

		ASA<-dmtA-SBA; # get 1 whereever the distances matrix is defined(under threshold) and not same block
		ASA<-as.spam(ASA)
	}
	mats_neigh<-list()
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\ndist:",limiteinf,"->",limitesup);

		dmtr<-dmtA # will be the matrix of neighbours for a ring(1/0)
		unselect<-which(dist_matA@entries>=limitesup | dist_matA@entries<limiteinf)
		dmtr@entries[unselect]<-rep(0,length(unselect))# 1 when dist_matA defined (distance<max(distance))
		dmtr<-as.spam(dmtr)
		cat(" mean nb Neigh",round(sum(dmtr/2)/dimension,digits=2))

		if(!is.null(group)){
			# ASr<-logic_and_spam(ASA,dmtr);
			# SBr<-logic_and_spam(SBA,dmtr);
			ASr<-ASA*dmtr;
			SBr<-SBA*dmtr;
			cat(" (SG:",round(sum(SBr/2)/dimension,digits=2),";DG:",round(sum(ASr/2)/dimension,digits=2),")")
			mats_neigh[[i]]<-list(dmtr=dmtr,SBr=SBr,ASr=ASr);
		}else{
			mats_neigh[[i]]<-list(dmtr=dmtr);
		}
	}
	cat("\n")
	attributes(mats_neigh)$breaks<-distances
	return(mats_neigh)
}

# compute structured morans I 
structured.moransI<-function(mats_neigh=NULL,raw.values,nb_rep_sign=0,rm.NA=TRUE,rand.sym=FALSE){
	# distances: breaks for the classes of distance to be tested for Moran's I
	# mats_neigh: for each class of distances the matrices of neighborhood (general, same block, different blocks)
	# values: analysed values
	# include_streets_anal: if compute the structured by streets moran's I
	# nb_rep_sign: number of repetitions if computing the significance
	# rm.NA: should NA be removed/ignored?

	## test validity of values
	zNA<-which(is.na(raw.values))
	zNoNA<-which(!is.na(raw.values))
	if(length(zNA)>0){
		if(!rm.NA){
			stop("values contains NA, if should be ignored, please use rm.NA=TRUE")
		}else{
			values<-raw.values[zNoNA]
		}
	}else{
		values<-raw.values
	}
	breaks<-attributes(mats_neigh)$breaks
	if(length(mats_neigh[[2]])==3){
		include_streets_anal<-TRUE
	}else{
		include_streets_anal<-FALSE
	}
	
	## calculate moran's I
	mI1 <- list() ; # list of moran test results general
	mI2 <- list() ; # list of moran test results in blocks
	mI3 <- list() ; # list of moran test results across streets
	difft<-NULL;
	nb_neigh<-mat.or.vec(length(breaks),3)

	if(nb_rep_sign>0){
		difft<-mat.or.vec(length(breaks),nb_rep_sign)
	}
	for (i in 2:(length(breaks))){
		limiteinf=breaks[i-1];
		limitesup=breaks[i];
		cat("\nlimiteinf=",limiteinf,"limitesup=",limitesup);
		# NB: any function can be put in the glist argument of nb2listw

		label<-paste(limiteinf,limitesup,sep="-")
		if(length(zNA)>0){
			dmtr<-mats_neigh[[i]][[1]][zNoNA,zNoNA]
		}else{
			dmtr<-mats_neigh[[i]][[1]]
		}
		mI1[[label]]<-moran.spam(dmtr,values);
		cat(" Moran's I:", mI1[[label]]);
		nb_neigh[i,1]<-sum(dmtr/2) # div by 2 as the matrix is sym
		if(include_streets_anal){
			if(length(zNA)>0){
				SBrtrue<-mats_neigh[[i]][[2]][zNoNA,zNoNA]
			}else{
				SBrtrue<-mats_neigh[[i]][[2]]
			}
			if(length(zNA)>0){
				ASr<-mats_neigh[[i]][[3]][zNoNA,zNoNA]
			}else{
				ASr<-mats_neigh[[i]][[3]]
			}
			mI2[[label]]<-moran.spam(SBrtrue,values);
			mI3[[label]]<-moran.spam(ASr,values);
			nb_neigh[i,2]<-sum(SBrtrue/2)
			nb_neigh[i,3]<-sum(ASr/2)
			if(nb_rep_sign>0){
				# get a 0 in entries of SBr any time dmtr is 1
				phantom<-dmtr
				phantom@entries<-rep(0,length(phantom@entries))
				SBrtrue<-phantom+as.spam(SBrtrue)
				for(repetition in 1:nb_rep_sign){
					#randomize SBrtrue
					if(rand.sym){
						SBr<-as.spam(symmetric.random.entries.by.line.spam.c(SBrtrue)) # randomize sym
					}else{
						SBr<-as.spam(random_spam_entries_by_row(SBrtrue)) # randomize not sym
					}
					#put complement in ASr
					ASr<-as.spam(dmtr-SBr);
					mI2t<-moran.spam(SBr,values);
					mI3t<-moran.spam(ASr,values);
					difft[i,repetition]<-mI2t-mI3t;
				}
			}
		}
	}
	cat("\n")
	sMI<-list(mI1=mI1,mI2=mI2,mI3=mI3,difft=difft,nb.neigh=nb_neigh)
	attributes(sMI)$breaks<-breaks
	return(sMI)
}
legend_position<-NULL
plot.structured.moransI<-function(mI,add=FALSE,neigh.struct=FALSE,plot=TRUE){
  	breaks<-attributes(mI)$breaks
	mI1<-mI[[1]]
	if(length(mI[[2]])>0){
		mI2<-mI[[2]]
		mI3<-mI[[3]]
		include_streets_anal<-TRUE
		if(!is.null(mI[[4]])){
		  plot.signif<-TRUE
		}else{
		  plot.signif<-TRUE
		}
	}else{
		include_streets_anal<-FALSE
	}

	## plot moran's I
	morans_I1 = {}
	morans_I2 = {}
	morans_I3 = {}
	med_position = {}
	legend_position = {}
	for (i in 1:length(mI1)){
		med_position=c(med_position,breaks[i]-eval(parse(text=(names(mI1[i]))))/2);
		legend_position=c(legend_position,names(mI1[i]));
		morans_I1 = c(morans_I1,mI1[[i]])
		if(include_streets_anal){
			morans_I2 = c(morans_I2,mI2[[i]])
			morans_I3 = c(morans_I3,mI3[[i]])
		}
	}

	plot.mI<-list(mI1=morans_I1,mI2=morans_I2,mI3=morans_I3,med_position=med_position)

	if(plot){
	  if(add==FALSE){
	    plot(c(breaks[1],breaks[(length(breaks))]),c(0.8*min(morans_I1,morans_I2,morans_I3),1.1*max(morans_I1,morans_I2,morans_I3)),type='n',xaxt='n',xlab="Distance class (m)",ylab="Morans's I")
	    axis(1,at=med_position,labels=legend_position)
	  }
	  lines(med_position,morans_I1,col=1) # black general
	  if(include_streets_anal){
	    lines(med_position,morans_I2,col=4) # blue within blocks
	    lines(med_position,morans_I3,col=2) # red inter blocks
	    legend("topright",c("Intra block","General","Inter blocks"),col=c("blue","black","red"),lty=1)
	    if(plot.signif){
	      dev.new()
	      plot.signif.structured.moransI(med_position,mI[[4]],morans_I1,morans_I2,morans_I3,legend_position)
	    }
	  }
	}
		
	legend_position<<-legend_position;
	attributes(plot.mI)$breaks<-breaks

	return(invisible(plot.mI))
}
plot.signif.structured.moransI<-function(med_position,difft,morans_I1,morans_I2,morans_I3,legend_position,main=NULL,pch=3,col=4,xlab="Distance class (m)",ylab=expression(paste(I[sg],"-",I[dg])),...){
	## plot diff AS/SB
	diffr<-morans_I2-morans_I3
	nb_rep_sign<-dim(difft)[2]
	if(is.null(main)){
		main=paste(nb_rep_sign,"permutations")
	}
	plot(c(min(med_position),max(med_position)),c(min(difft,na.rm=TRUE),max(diffr,na.rm=TRUE)),type="n",xlab=xlab,xaxt="n",ylab=ylab,main=main,pch=pch,...)
	axis(1,at=med_position, lab=legend_position)
	lines(med_position,diffr,col=col,type="b",pch=pch,...)
	# ## get the 95% confidence lines 
	# ordered_difft<-mat.or.vec(length(med_position),nb_rep_sign)
	# for(i in 2:length(med_position)){
	# 	ordered_difft[i,]<-difft[i,order(difft[i,])]
	# }
	boxplot.free(t(difft[-1,]),add=TRUE,xaxt="n",at=med_position,pars=list(boxwex=(med_position[2]-med_position[1])/3),pch=pch)
	
	return()
}
get.med.position.axis<-function(distances){
	left.lim<-distances[1:(length(distances)-1)]
	right.lim<-distances[2:length(distances)]
	med.position<-(right.lim+left.lim)/2
	labels.axis<-paste(as.character(left.lim),"-",as.character(right.lim),sep="")

	axis(1,at=seq(1:length(left.lim)), lab=labels.axis)
	return(med.position)
}
zgen<-function(detection,zNA,Q=NULL,Kv=NULL,Ku=NULL,mu=NULL,c.comp=NULL,est.v=NULL,est.u=NULL,est.w=NULL,force.mu=FALSE){
  # generate the positiveness according to estimates
  # Necessary:
  # detection: quality of detection for each point
  # zNA: indices of the the non-observed  

  # Recommended:
  # Q: spatial precision matrix between points (necessary if including a spatial component)
  # Kv: estimated non-spatial precision
  # Ku: estimated spatial precision factor
  # mu: estimated mean of the spatial field
  # c.comp: estimated cumulated effect of the cofactors (per point)

  # Not recommended but for testing purposes:
  # est.v: fixed estimated non-spatial component of the field (per point)
  # est.u: fixed estimated spatial component of the filed (per point)
  # est.w: fixed estimated field (per point, if given everything above is ignored)
  # force.mu: force the mean of the non-spatial component to be at the given mu value

	dimension<-length(detection)
	u.p<-rep(0,dimension)
	c.p<-rep(0,dimension)
	v.p<-rep(0,dimension)
	if(is.null(est.w)){
		if(is.null(est.u)){
			if(!is.null(Ku) && !is.null(Q) && !is.null(mu)){
				u.p <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mu,dimension), Q=Ku*Q));
				if(force.mu){
					u.p <-(u.p-mean(u.p)+mu)
				}
				cat("mean u.p",mean(u.p),"sd u.p",sd(u.p),"\n")
			}else{
				cat("Use no u component\n");
			}
		}else{
			u.p<-est.u
		}
		if(!is.null(c.comp)){
			c.p<-c.comp
			cat("mean c.p",mean(c.comp),"sd c.p",sd(drop(c.comp)),"\n")
		}else{
			cat("Use no c component\n");
		}
		if(is.null(est.v)){
			if(!is.null(Kv)){
				Qvp<-diag.spam(Kv,dimension)
				if(is.null(est.v)){
					est.v<-rep(0,dimension)
				}
				v.p <- drop(rmvnorm.prec.pseudo(n=1, mu=est.v, Q=Qvp));
				cat("mean v.p",mean(v.p),"sd v.p",sd(v.p),"\n")
			}else{
				cat("Use no v component\n");
			}
		}else{
			v.p+est.v
		}
		w.p<-u.p+c.p+v.p;
	}else{
		w.p<-est.w
	}

	y.p <- rnorm(n=dimension, mean=w.p, sd=1);
	cat("mean y.p",mean(y.p),"sd y.p",sd(y.p),"\n")
	z.p <-generate_z(y.p,detection,zNA)
	cat("mean z.p",mean(z.p[!is.na(z.p)]),"\n")

	return(list(z=z.p,y=y.p,w=w.p,v=v.p,c=c.p,u=u.p))
}
gen.map<-function(db,mu=-1,Ku=1,Kv=10,f=10,T=1,kern=expKernel,obs.qual=1,c.val=NULL,randomizeNA=FALSE,threshold=50,epsilon=0.01){
  # simpler wrapper calling zgen
  # for map generations

  # db: data.frame with in columns at least the X,Y of the points
  # Ku: spatial component precision factor
  # Kv: non-spatial component precision factor
  # f: characteristic distance of the kernel
  # T: factor for streets (<1 => barrier effect of streets) imply 
  #    the presence of a GroupNum column in db
  # kern: the kernel to use (expKernel,gaussianKernel, cauchyKernel or geometricKernel). If f and T have been estimated the kernel here should correspond to the one used to estimate them.
  # obs.qual: vector of the quality of the observers, default to 1 (all perfect). It can be set to a real in [0,1], applied to all points; a vector the same length than db, each then applied to its point; or be a vector with named values corresponding to the levels in db$IdObserver, then applied to these same Observers. In the last case, if observers do not have an estimate in obs.qual there value is fixed to the mean of obs.qual

  nPoints<-dim(db)[1]
  # get the Q precision matrix
  cat("Generating the distance matrix...\n")
  dist_mat <-nearest.dist(x=db[,c("X","Y")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
  diag(dist_mat)<- rep(0,dim(dist_mat)[1])
  dbout<-db[,c("X","Y")]

  # Missing information
  if(randomizeNA){
    cat("Generating missingness...\n")
    nNA<-sum(db$observed==0)
    zNA<-sample(nNA,1:nPoints)
  }else{
    cat("Keeping missingness\n")
    zNA<-which(db$observed==0)
  }
  dbout$observed<-rep(0,dim(dbout)[1])
  dbout$observed[-zNA]<-1

  # spatial groupings
  if(!is.null(db$GroupNum)){
    dbout$GroupNum<-db$GroupNum # cannot generate that yet
    cat("Keeping the barriers structure\n")
    SB <- nearest.dist(x=cbind(db$GroupNum,rep(0,length(db$GroupNum))), method="euclidian", upper=NULL,delta=0.1)
    SB@entries<-rep(1,length(SB@entries))

    dmt<-dist_mat
    dmt@entries<-rep(1,length(dmt@entries))# [dmt@entries!=0]<-1 # 1 only when dist_mat not 0
    SB<-as.spam(SB*dmt);
  }else{
    if(T!=1){
      stop("T !=1 and no GroupNum in db. Aborting gen.mat().\n")
    }
    SB<-1
  }
  Q<-QfromfT(dist_mat,SB,f,T,kern=kern,addeps=epsilon)

  # cofactors
  if(is.null(c.val)){
    c.comp<-NULL
  }else{
    cat("Keeping the cofactors structure\n")
    cofs<-names(c.val)
    c.map<-as.matrix(db[,cofs])
    dbout<-cbind(dbout,db[,cofs])
    c.comp<-drop(c.map%*%c.val)
  }

  # get the observation quality vector
  dbout$IdObserver<-db$IdObserver
  if(length(obs.qual)==nPoints){
    cat("Keeping observation quality\n")
    obs.vect<-obs.qual
  }else if(length(obs.qual)==1){
    cat("Uniform observation quality:",obs.qual,"\n")
    obs.vect<- rep(obs.qual,nPoints)
  }else if(!is.null(names(obs.qual))){
    cat("Generating observation quality on IdObserver structure\n")
    if(is.null(db$IdObserver)){
      cat("Missing IdObserver in db. Abort.\n")
      return(NULL)
    }
    # make index of observers each with its sensitivity (ObsValues)
    ObsNames<-levels(as.factor(db$IdObserver))
    ObsValues<-rep(0,length(ObsNames))
    names(ObsValues)<-ObsNames
    ObsValues[names(obs.qual)]<-obs.qual
    ObsValues[setdiff(ObsNames,names(obs.qual))]<-mean(obs.qual)

    # apply ObsValues to the visited houses
    obs.vect<-rep(0,nPoints)
    for(numo in 1:length(ObsNames)){
      nameo<-names(ObsValues)[numo]
      obs.vect[which(db$IdObserver==nameo)]<- ObsValues[nameo]
    }
  }else{
    stop("obs.val not formated as expected. Aborting.\n")
  }
  dbout$bgen<-obs.vect

  # generate spatial autocorrelation and final values
  cat("Generating autocorrelation...\n")
  simul<-zgen(obs.vect,zNA,Q=Q,Kv=Kv,Ku=Ku,mu=mu,c.comp=c.comp)
  dbout$positive<-simul$z
  dbout$ygen<-simul$y
  dbout$wgen<-simul$w
  dbout$vgen<-simul$v
  dbout$cgen<-simul$c
  dbout$ugen<-simul$u

  attributes(dbout)$Q<-Q
  attributes(dbout)$Kv<-Kv
  attributes(dbout)$Ku<-Ku
  attributes(dbout)$mu<-mu

  return(dbout)
}
# simul<-gen.map(db)
# plot_reel(db$X,db$Y,simul$z)


generated.morans.struct<-function(distances,mats_neigh,nbRep,est.detection=NULL,est.Q=NULL,est.Ku=NULL,est.Kv=NULL,est.mu=NULL,est.c.comp=NULL,est.v=NULL,true.val=NULL,est.u=NULL,est.w=NULL,force.mu=FALSE,trueStatus=FALSE){
	if(is.null(est.detection)){
		dimension<-max(c(length(est.u),length(est.v),length(est.w),dim(est.Q)[1]))

		est.detection<-rep(1,dimension)
	}
	if(is.null(true.val)){
		sel<-1:length(est.detection)
		zNA<-c()
	}else{
		sel<-(true.val!=9)
		zNA<-which(true.val==9)
	}
	IdMoinsIs<-mat.or.vec(nbRep,length(distances)-1)
	MI1<-mat.or.vec(nbRep,length(distances)-1)
	MI2<-mat.or.vec(nbRep,length(distances)-1)
	MI3<-mat.or.vec(nbRep,length(distances)-1)
	for(i in 1:nbRep){
		cat(i,"th generation\n")
		generated<-zgen(est.detection,zNA,est.w=est.w,est.u=est.u,est.v=est.v,est.Q=est.Q,est.Ku=est.Ku,est.Kv=est.Kv,est.mu=est.mu,est.c.comp=est.c.comp,force.mu=force.mu);
		if(trueStatus){
			y.p<-generated$y
			cat("l y.p:",length(y.p));
			mIrefPseudoObs<-structured.moransI(mats_neigh,y.p,nb_rep_sign=0);
		}else{
			z.p<-generated$z
			mIrefPseudoObs<-structured.moransI(mats_neigh,z.p[sel],nb_rep_sign=0);
		}
		plotIval<-plot.structured.moransI(distances,mIrefPseudoObs,plot=FALSE);
		IdMoinsIs[i,]<-plotIval$mI2-plotIval$mI3
		MI2[i,]<-plotIval$mI2
		MI3[i,]<-plotIval$mI3
		MI1[i,]<-plotIval$mI1
	}
	# meanIdMoinsIs<-apply(IdMoinsIs,2,mean)
	# limIdMoinsIs<-apply(IdMoinsIs,2,quantile,prob=c(0.025,0.975))

	boxplot.free(MI1,ylim=c(0,max(MI1)),xaxt="n")
	plotTrueIval<-NULL
	if(!is.null(true.val)){
		mIrefData<-structured.moransI(mats_neigh,true.val[sel]);
		plotTrueIval<-plot.structured.moransI(distances,mIrefData,plot=FALSE);
		trueIdMoinsIs<-plotTrueIval$mI2-plotTrueIval$mI3
		lines(plotTrueIval$mI1,col=4)
	}
	get.med.position.axis(distances)
	boxplot.free(IdMoinsIs,ylim=c(0,max(IdMoinsIs)),xaxt="n")
	if(!is.null(true.val)){
		lines(trueIdMoinsIs,col=4)
	}
	get.med.position.axis(distances)
	return(list(MI1=MI1,MI2=MI2,MI3=MI3,ref=plotTrueIval))
}

replot.gen.struct<-function(distances,MIs,ylim1=NULL,ylim2=NULL){
	MI1<-MIs$MI1
	MI2<-MIs$MI2
	MI3<-MIs$MI3
	if(is.null(ylim1)){
		ylim1<-c(0,max(MI1))
	}
	boxplot.free(MI1,ylim=ylim1,xaxt="n")
	if(!is.null(MIs$ref)){
		lines(MIs$ref$mI1,col=1)
	}
	dumb<-get.med.position.axis(distances)
	IdMoinsIs<-MI2-MI3
	if(is.null(ylim2)){
		ylim2<-c(0,max(IdMoinsIs))
	}
	boxplot.free(IdMoinsIs,ylim=ylim2,xaxt="n")
	if(!is.null(MIs$ref)){
		lines(MIs$ref$mI2-MIs$ref$mI3,col=1)
	}
	dumb<-get.med.position.axis(distances)
	return()
}

MeanCovPairAtDist<-function(mask,CovMat){
	linkedCovMat<-mask*CovMat
	linkedCovMat<-as.spam(linkedCovMat)
	meanCovPair<-mean(linkedCovMat@entries)

	return(meanCovPair)
}

StructCorrel<-function(distances,Q=NULL,CovMat=NULL,mats_neigh=NULL){
	if(is.null(CovMat)){
		if(!is.null(Q)){
			CovMat<-solve(Q)
		}else{
			stop("Need Q or CovMat\n")
		}
	}
	if(length(mats_neigh[[2]])==3){
		include_streets_anal<-TRUE
	}else{
		include_streets_anal<-FALSE
	}

	SC <- data.frame() ; # list of moran test results general
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\nlimiteinf=",limiteinf,"limitesup=",limitesup);
		# NB: any function can be put in the glist argument of nb2listw

		label<-paste(limiteinf,limitesup,sep="-")
		dmtr<-mats_neigh[[i]][[1]]
		SC[i,1]<-MeanCovPairAtDist(dmtr,CovMat);
		if(include_streets_anal){
			SBr<-mats_neigh[[i]][[2]]
			ASr<-mats_neigh[[i]][[3]]
			SC[i,2]<-MeanCovPairAtDist(SBr,CovMat);
			SC[i,3]<-MeanCovPairAtDist(ASr,CovMat);
		}
		# print(SC)
	}

	return(SC)
}
StructNeigh<-function(distances,mats_neigh){
	
	nb_mat<-length(mats_neigh[[2]])
	nb_neigh<-as.data.frame(mat.or.vec(length(distances),nb_mat+1))
	if(nb_mat==3){
		include_streets_anal<-TRUE
		names(nb_neigh)<-c("nbNeighTotal","med_position","nbNeighSB","nbNeighAS")
	}else{
		include_streets_anal<-FALSE
		names(nb_neigh)<-c("nbNeighTotal","med_position")
	}

	SC <- data.frame() ; # list of moran test results general
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\nlimiteinf=",limiteinf,"limitesup=",limitesup);
		# NB: any function can be put in the glist argument of nb2listw

		label<-paste(limiteinf,limitesup,sep="-")
		row.names(nb_neigh)[i]<-label
		dmtr<-mats_neigh[[i]][[1]]
		nb_neigh[i,1]<-sum(dmtr/2) # div by 2 as the matrix is sym
		nb_neigh[i,2]<- (limiteinf+limitesup)/2
		if(include_streets_anal){
			SBr<-mats_neigh[[i]][[2]]
			ASr<-mats_neigh[[i]][[3]]
			nb_neigh[i,3]<-sum(SBr/2)
			nb_neigh[i,4]<-sum(ASr/2)
		}
		# print(SC)
	}

	return(nb_neigh[-1,])
}




#===============================
# Map generation
#===============================
make_block_points<-function(nb_per_sidex,nb_per_sidey,xmin,xmax,ymin,ymax,ratio_street_dist=2){
	# return the coordinates of the houses in a block
	if(nb_per_sidex<2||nb_per_sidey<2){
		cat("Cannot handle a nb_per_side=(",nb_per_sidex,nb_per_sidey,")\n");
		return();
	}
	nb_per_block<-2*(nb_per_sidex+nb_per_sidey)-4

	block_points<-mat.or.vec(nb_per_block,2);
	
	intervalx<-(xmax-xmin)/(nb_per_sidex-1+ratio_street_dist);
	intervaly<-(ymax-ymin)/(nb_per_sidey-1+ratio_street_dist);
	street_intx<-intervalx*ratio_street_dist/2
	street_inty<-intervaly*ratio_street_dist/2
	# cat("intervalx",intervalx,"intervaly",intervaly,"street_int",street_int,"\n");
	# first line
	block_points[1:nb_per_sidex,1]<-seq(1,nb_per_sidex)*intervalx+xmin-intervalx+street_intx;
	block_points[1:nb_per_sidex,2]<-rep(ymin+street_inty);

	# middle points
	for(num_line in 1:(nb_per_sidey-2)){
		yline<-ymin+street_inty+(num_line)*intervaly;
		block_points[nb_per_sidex+2*num_line-1,1]<-xmin+street_intx;
		block_points[nb_per_sidex+2*num_line-1,2]<-yline;
		block_points[nb_per_sidex+2*num_line,1]<-xmax-street_intx;
		block_points[nb_per_sidex+2*num_line,2]<-yline;
	}

	# last line
	block_points[((nb_per_block-nb_per_sidex+1):nb_per_block),1]<-block_points[1:nb_per_sidex,1];
	block_points[((nb_per_block-nb_per_sidex+1):nb_per_block),2]<-rep(ymax-street_inty);

	return(block_points);
}

make_data_points<-function(nb_per_sidex,nb_blocks_x,nb_blocks_y,points_dist,nb_per_sidey=nb_per_sidex,ratio_street_dist=2){
	# make a map with points in blocks, starting at 0,0
	data_points={};
	num_block=1;
	block_width=nb_per_sidex*(points_dist+1)
	block_height=nb_per_sidey*(points_dist+1)
	for(num_line_block in 1:nb_blocks_x){
		for(num_col_block in 1:nb_blocks_y){
			block<-make_block_points(nb_per_sidex,nb_per_sidey,
				block_width*(num_line_block-1),
				block_width*num_line_block,
				block_height*(num_col_block-1),
				block_height*num_col_block,
				ratio_street_dist=ratio_street_dist);
			data_points<-rbind(data_points,cbind(block,rep(num_block,nrow(block))));
			num_block=num_block+1;
		}
	}
	data_points<-as.data.frame(data_points)
	names(data_points)<-c("X","Y","GroupNum")
	
	return(data_points)
}

## generation of cofactors
gen_c.map<-function(nbfact,nbpoints,rates=rep(0.5,nbfact)){
  c.map={}
  for(fact in 1:nbfact){
    prob<-c(1-rates[fact],rates[fact])
    c.col<-sample(0:1,nbpoints,replace=TRUE,prob=prob)
    c.map<-cbind(c.map,c.col);
  }
  return(c.map);
}
# gen_c.map(3,10,rates=c(0.1,0.9,0.5))

#===============================
# Functions specific to the GMRF
#===============================
#-------------------------------
# general sampling functions
#-------------------------------

# files for the saving of parameters
betafile<-"betasamples.txt"
coffile<-"cofactors.txt"
monitorfile<-"sampled.txt"

# adapt the standard deviation of the proposal
adaptSDProp <- function(sdprop, accepts, lowAcceptRate=0.15, highAcceptRate=0.4,tailLength=20){
  # according to Gelman1996
  # adapt the sampling deviation so that acceptance rate fall within:
  # 0.15 and 0.4 (careful to begin the sampling after that)
  acceptRate <- mean(tail(accepts,tailLength))
  # cat("accept rate:",acceptRate);

  attributes(sdprop)$noupdate<-FALSE
  if(acceptRate < lowAcceptRate){
    newsdprop<-sdprop*0.9
    # cat("update sdprop",sdprop,"to",newsdprop);
    return(newsdprop)
  }else if(acceptRate > highAcceptRate){
    newsdprop<-sdprop*1.1
    # cat("update sdprop",sdprop,"to",newsdprop);
    return(newsdprop)
  }else{
    attributes(sdprop)$noupdate<-TRUE
    return(sdprop)
  }
}
# Example:
# 	logsdTprop<-adaptSDProp(logsdTprop,acceptT)
#	adaptOK<-adaptOK && attributes(logsdTprop)$noupdate

setClass("cb.diag.out",representation(colAnalysed = "vector", geweke= "vector",limitGeweke="numeric",raftery="matrix"),contains="data.frame")

gibbsit <- function(data, varnames = NULL, q = 1/40, r =
	  1/80, s = 19/20, epsilon = 1/1000, spacing = 1, resfile = "",
	  warning=TRUE,NminOnly=FALSE){
# Version 1.1: August 15, 1995
# as retrieved on http://madison.byu.edu/bayes/gibbsit.txt August 3 2011
# slightly adapted for R by Corentin Barbu corentin.barbu AT gmail.com
# changes are between ## CB
#
# An S code translation of the Gibbsit FORTRAN program by Adrian Raftery
# and Steven Lewis.  This S function was developed by Steven Lewis from
# an earlier function originally written by Karen Vines.
#
# Developer:  Steven M. Lewis (slewis@stat.washington.edu).
# This software is not formally maintained, but I will be happy to
# hear from people who have problems with it, although I cannot
# guarantee that all problems will be corrected.
#
# Permission is hereby granted to StatLib to redistribute this
# software.  The software can be freely used for non-commercial
# purposes, and can be freely distributed for non-commercial
# purposes only.  The copyright is retained by the developer.
#
# Copyright 1995  Steven M. Lewis, Adrian E. Raftery and Karen Vines
#
#
# References:
#
# Raftery, A.E. and Lewis, S.M. (1992).  How many iterations in the
# Gibbs sampler?  In Bayesian Statistics, Vol. 4 (J.M. Bernardo, J.O.
# Berger, A.P. Dawid and A.F.M. Smith, eds.). Oxford, U.K.: Oxford
# University Press, 763-773.
# This paper is available via the World Wide Web by linking to URL
#   http://www.stat.washington.edu/tech.reports and then selecting
# the "How Many Iterations in the Gibbs Sampler" link.
# This paper is also available via regular ftp using the following
# commands:
#   ftp ftp.stat.washington.edu (or 128.95.17.34)
#   login as anonymous
#   enter your email address as your password
#   ftp> cd /pub/tech.reports
#   ftp> get raftery-lewis.ps
#   ftp> quit
#
# Raftery, A.E. and Lewis, S.M. (1992).  One long run with diagnos-
# tics: Implementation strategies for Markov chain Monte Carlo.
# Statistical Science, Vol. 7, 493-497.
#
# Raftery, A.E. and Lewis, S.M. (1995).  The number of iterations,
# convergence diagnostics and generic Metropolis algorithms.  In
# Practical Markov Chain Monte Carlo (W.R. Gilks, D.J. Spiegelhalter
# and S. Richardson, eds.). London, U.K.: Chapman and Hall.
# This paper is available via the World Wide Web by linking to URL
#   http://www.stat.washington.edu/tech.reports and then selecting
# the "The Number of Iterations, Convergence Diagnostics and Generic
# Metropolis Algorithms" link.
# This paper is also available via regular ftp using the following
# commands:
#   ftp ftp.stat.washington.edu (or 128.95.17.34)
#   login as anonymous
#   enter your email address as your password
#   ftp> cd /pub/tech.reports
#   ftp> get raftery-lewis2.ps
#   ftp> quit
#
#
# Input arguments:
# data     = matrix of output from MCMC run; rows are the results for each
#            iteration (assumed to be in a consecutive order), columns
#            being the different variables.
# varnames = vector of names to be tested.  These names must correspond to
#            column names in the data matrix.
## CB
# q 	   = the quantile to be estimated
# r 	   = the desired margin of error of the estimate
# s	   = the probability of obtaining an estimate in the interval (q-r,q+r)
# epsilon  = Precision required for estimate of time to convergence
#  Example values of q, r, s:                                          
#     0.025, 0.005,  0.95 (for a long-tailed distribution)             
#     0.025, 0.0125, 0.95 (for a short-tailed distribution);           
#     0.5, 0.05, 0.95;  0.975, 0.005, 0.95;  etc.                      
#                                                                      
#  The result is quite sensitive to r, being proportional to the       
#  inverse of r^2.                                                     
#                                                                      
#  For epsilon, we have always used 0.001.  It seems that the result   
#  is fairly insensitive to a change of even an order of magnitude in  
#  epsilon.                                                            
## CB
#                                                                      
#  One way to use the program is to run it for several specifications  
#  of r, s and epsilon and vectors q on the same data set.  When one   
#  is sure that the distribution is fairly short-tailed, such as when  
#  q=0.025, then r=0.0125 seems sufficient.  However, if one is not    
#  prepared to assume this, safety seems to require a smaller value of 
#  r, such as 0.005.                                                   
#
# spacing  = Frequency at which the MCMC output was recorded; this will
#            usually be 1.  If only every other MCMC iteration was retained,
#            then spacing should be set to 2, etc.
# resfile  = Name of an output file to which the results are to be written.
#
#
# Output:
# A matrix whose rows are the variables tested and whose columns contain the
# following results:
# 1st col. = Kthin, the thinning parameter required to make the chain first-
#            order Markov.
# 2nd col. = Nburn, the number of iterations needed for the burn-in.
# 3rd col. = Nprec, the number of iterations required to achieve the specified
#            precision.
# 4th col. = Nmin, the number of iterations required if they were independent.
# 5th col. = I_RL, the ratio of Nburn plus Nprec to the number of iterations
#            required using independent sampling
#            (see Raftery and Lewis, StatSci 1992).
# 6th col. = Kind, the thinning parameter required to make the chain into an
#            independence chain.

	phi <- qnorm((s + 1)/2)
	nmin <- as.integer(ceiling((q * (1 - q) * phi^2)/r^2))
## CB
	if(NminOnly==TRUE){
		return(nmin);
	}
	if(length(data[,1])>45000){
		if(warning==TRUE){
			cat("\nR limitations on integer manipulations doesn't allow to handle all this matrix\n analysis limited to the last 45000 iterations\n")
		}
		data<-data[-(1:(length(data[,1])-45000)),]
	}	
	if(is.null(varnames)){
		data<-as.data.frame(data)
		varnames=dimnames(data)[[2]]
		if(is.null(varnames)){
			varnames<- as.character(seq(1:dim(data)[2]))
		}
	}
## CB
	resmatrix <- matrix(NA, nr = length(varnames), nc = 6)
	iteracnt <- nrow(data)
	if(resfile != "")
		cat("Results when q = ", q, ", r = ", r, ", s = ", s,
			", epsilon = ", epsilon, ":\n\n", file = resfile, sep
			 = "",append=TRUE)
	for(r1 in 1:length(varnames)) {
		quant <- quantile(data[, (curname <- varnames[r1])], probs = q)
		dichot <- as.integer(data[, curname] <= quant)	#
#	First find the actual thinning parameter, kthin.
		testtran <- array(as.integer(0), dim = c(2, 2, 2))
		kwork <- 0
		bic <- 1
		while(bic >= 0) {
			kwork <- as.integer(kwork + 1)
			testres <- dichot[seq(1, iteracnt, by = kwork)]
			thindim <- length(testres)
			tttemp <- as.integer(testres[1:(thindim - 2)] + testres[
				2:(thindim - 1)] * 2 + testres[3:thindim] * 4)
			testtran[1, 1, 1] <- sum(as.integer(tttemp == 0))
			testtran[2, 1, 1] <- sum(as.integer(tttemp == 1))
			testtran[1, 2, 1] <- sum(as.integer(tttemp == 2))
			testtran[2, 2, 1] <- sum(as.integer(tttemp == 3))
			testtran[1, 1, 2] <- sum(as.integer(tttemp == 4))
			testtran[2, 1, 2] <- sum(as.integer(tttemp == 5))
			testtran[1, 2, 2] <- sum(as.integer(tttemp == 6))
			testtran[2, 2, 2] <- sum(as.integer(tttemp == 7))
			g2 <- 0
			for(i1 in 1:2) {
				for(i2 in 1:2) {
				  for(i3 in 1:2) {
				    if(testtran[i1, i2, i3] != 0) {
				      fitted <- (sum(testtran[i1, i2, 1:2]) *
				        sum(testtran[1:2, i2, i3]))/(sum(
				        testtran[1:2, i2, 1:2]))
				      g2 <- g2 + testtran[i1, i2, i3] * log(
				        testtran[i1, i2, i3]/fitted) * 2
				    }
				  }
				}
			}
			bic <- g2 - log(thindim - 2) * 2
		}
		kthin <- as.integer(kwork * spacing)	#
#	Now determine what the thinning parameter needs to be to achieve
#	independence, kmind.
		indptran <- matrix(as.integer(0), nr = 2, nc = 2)
		firsttime <- TRUE
		bic <- 1
		while(bic >= 0) {
			if(!firsttime) {
				kwork <- as.integer(kwork + 1)
				testres <- dichot[seq(1, iteracnt, by = kwork)]
				thindim <- length(testres)
			}
			indptemp <- as.integer(testres[1:(thindim - 1)] +
				testres[2:thindim] * 2)
			indptran[1, 1] <- sum(as.integer(indptemp == 0))
			indptran[2, 1] <- sum(as.integer(indptemp == 1))
			indptran[1, 2] <- sum(as.integer(indptemp == 2))
			indptran[2, 2] <- sum(as.integer(indptemp == 3))
			if(firsttime) {
#	Save this particular transfer matrix for later use in calculating the
#	length of the burn-in and for the required precision.
				finaltran <- indptran
				firsttime <- FALSE
			}
			den <- rep(apply(indptran, 1, sum), 2) * rep(apply(
				indptran, 2, sum), c(2, 2))
			g2 <- sum(log((indptran * (dcm1 <- thindim - 1))/den) *
				indptran, na.rm = TRUE) * 2
			bic <- g2 - log(dcm1)
		}
		kmind <- as.integer(kwork * spacing)	#
#	Next find the length of burn-in and the required precision.
		alpha <- finaltran[1, 2]/(finaltran[1, 1] + finaltran[1, 2])
		beta <- finaltran[2, 1]/(finaltran[2, 1] + finaltran[2, 2])
		tempburn <- log((epsilon * (alpha + beta))/max(alpha, beta))/(
			log(abs(1 - alpha - beta)))
		nburn <- as.integer(ceiling(tempburn) * kthin)
		tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2)/(((
			alpha + beta)^3) * r^2)
		nprec <- as.integer(ceiling(tempprec) * kthin)
		iratio <- (nburn + nprec)/nmin
		kind <- max(floor(iratio + 1), kmind)
		resmatrix[r1,  ] <- c(kthin, nburn, nprec, nmin, round(iratio,
			digits = 2), kind)
		if(resfile != "")
			cat("  ", curname, "\tKthin = ", kthin, "\tNburn = ",
				nburn, "\tNprec = ", nprec, "\tNmin = ", nmin,
				"\tI_RL = ", resmatrix[r1, 5], "\tKind = ",
				kind, "\n", file = resfile, sep = "", append =
				TRUE)
	}
	dimnames(resmatrix) <- list(varnames, c("Kthin", "Nburn", "Nprec",
		"Nmin", "I_RL", "Kind"))
	return(resmatrix)
}
Geweke.Diagnostic <- function (db,frac1=0.1,frac2=0.5) {
  x <- as.matrix(db)
  startx <- 1
  endx <- nrow(x)
  xstart <- c(startx, endx - frac2 * {endx - startx})
  xend <- c(startx + frac1 * {endx - startx}, endx)

  y.variance <- y.mean <- vector("list", 2)
  for (i in 1:2) {
    y <- x[xstart[i]:xend[i], ]
    y.mean[[i]] <- colMeans(as.matrix(y))
    yy <- as.matrix(y)
    y <- as.matrix(y)
    max.freq <- 0.5
    order <- 1
    max.length <- 200
    if (nrow(yy) > max.length) {
      batch.size <- ceiling(nrow(yy)/max.length)
      yy <- aggregate(ts(yy, frequency = batch.size), nfreq = 1, 
		      FUN = mean)
    }
    else {
      batch.size <- 1
    }
    yy <- as.matrix(yy)
    fmla <- switch(order + 1, spec ~ one, spec ~ f1, spec ~ 
		   f1 + f2)
    if (is.null(fmla)) 
      stop("Invalid order.")
    N <- nrow(yy)
    Nfreq <- floor(N/2)
    freq <- seq(from = 1/N, by = 1/N, length = Nfreq)
    f1 <- sqrt(3) * {
      4 * freq - 1
    }
    f2 <- sqrt(5) * {
      24 * freq * freq - 12 * freq + 1
    }
    v0 <- numeric(ncol(yy))
    for (j in 1:ncol(yy)) {
      cat("Gew:",names(db)[j],"\n")
      zz <- yy[, j]
      zz<-zz[which(is.finite(zz))]
      if (var(zz) == 0) {
	v0[j] <- 0
      }
      else {
	yfft <- fft(zz)
	spec <- Re(yfft * Conj(yfft))/N
	spec.data <- data.frame(one = rep(1, Nfreq), 
				f1 = f1, f2 = f2, spec = spec[1 + {
				  1:Nfreq
				}], inset = I(freq <= max.freq))
		glm.out <- try(glm(fmla, family = Gamma(link = "log"), 
				   data = spec.data), silent = TRUE)
	if (!inherits(glm.out, "try-error")) 
	  v0[j] <- predict(glm.out, type = "response", 
			   newdata = data.frame(spec = 0, one = 1, f1 = -sqrt(3), 
						f2 = sqrt(5)))
      }
    }
    spec <- list(spec = v0)
    spec$spec <- spec$spec * batch.size
    y.variance[[i]] <- spec$spec/nrow(y)
  }
  z <- {
    y.mean[[1]] - y.mean[[2]]
  }/sqrt(y.variance[[1]] + y.variance[[2]])
  print(cbind(y.mean[[1]],y.mean[[2]],y.variance[[1]],y.variance[[2]]))

  return(z)
}
basicGeweke<-function(db,frac1=0.1,frac2=0.5){
  x<-as.matrix(db)
  startx <- 1
  endx <- nrow(x)
  xstart <- c(startx, endx - frac2 * {endx - startx})
  xend <- c(startx + frac1 * {endx - startx}, endx)
  ymean<-yvar<-mat.or.vec(dim(db)[2],2)
  for(samp in 1:2){
    for(fact in 1:dim(db)[2]){
      y<-db[xstart[samp]:xend[samp],fact]
      ymean[fact,samp]<-mean(y)
      yvar[fact,samp]<-var(y)
    }
  }
  zscore<-(ymean[,1]-ymean[,2])/sqrt(yvar[,1]+yvar[,2])
  pvalue<-pnorm(-abs(zscore))
  return(pvalue)
}

cb.diag<-function(sampBrut,baseLimitGeweke=0.05,KthinInit=1,logfile=""){
	enoughIt<-FALSE;
	nbItEff<-0
	resG<-0
	discard<-0
	limitGeweke<-0
	nbItBrut<-length(sampBrut[,1])
	sampThined<-sampBrut[seq(1,nbItBrut,KthinInit),];
	# eliminate constant parameters

	colAnalysed<-which(apply(sampBrut,2,var)!=0)
	if(length(colAnalysed)<1){
		stop("cb.diag:all columns of sampBrut are constant\n")
	}else{
		sampThined<-sampThined[,colAnalysed]
	}
	
	nbIt<-length(sampThined[,1])
	# cat("sampThined:\n");
	# print(str(sampThined))

	## find the minimum number of iterations according to the raftery
	resRafLeft<-gibbsit(sampThined,warning=FALSE,resfile=logfile)
	# print(resRafLeft)
	# resRafMean<-gibbsit(sampThined,q=20/40,warning=FALSE) ## much too difficult to reach
	# print(resRafMean)
	resRafRight<-gibbsit(sampThined,q=39/40,warning=FALSE,resfile=logfile)
	# print(resRafRight)
	resRaf<-rbind(resRafLeft,resRafRight)
	resRafMax<-apply(resRaf,2,max,na.rm=TRUE)
	nbItMin<-as.integer(resRafMax[2]+resRafMax[3])*KthinInit;
	Kthin<-resRafMax[1]*KthinInit;
	Kind<-resRafMax[6]*KthinInit;
	cat("new Kthin:",Kthin,"Kind",Kind,"\n",file=logfile,append=TRUE);
	burnIn<-resRafMax[2]*KthinInit;
	Nmin<-resRafMax[4];
	if(nbItMin<nbItBrut){
		cat("Raftery positive (",nbItBrut,">",nbItMin,")\n",file=logfile,append=TRUE);
		## if enough iterations according to the raftery, test if Geweke ok on non burnin 
		discard<-burnIn;
		result<-TRUE;
		nbRep<-1
		nMaxRep<-5
		# account for the multiple variables 
		# that need to pass the test
		limitGeweke<-1-(1-baseLimitGeweke)^(nMaxRep/dim(sampThined)[2])
		while( nbRep<=nMaxRep && min(resG)<limitGeweke && class(result)!="try-error"){
			resGbrut<-1+limitGeweke # something bigger than limitGeweke
			discard<-burnIn+round((nbIt-burnIn)*(nbRep-1)*0.1);
			gewekeSample<-sampThined[-(1:discard),]
			cat("test Geweke with",discard,"discarded and limit=",limitGeweke,"\n",file=logfile,append=TRUE);
			result<-try(resGbrut<-Geweke.Diagnostic(gewekeSample,frac1=0.5),silent=FALSE) # return p-values, 
			resG<-pnorm(-abs(resGbrut))
			print(resG)

			nbRep<-nbRep+1
		}
		if(min(resG)<limitGeweke){
			## if the convergence is not ok according to Geweke
		  	## increase the number of iterations substentially
			cat("Geweke not ok:",resG,"(at least one <",limitGeweke,")\n",file=logfile,append=TRUE);
			newNbIt<-ceiling(nbItBrut*1.5);
		}else{
			## if the convergence is ok test if we have enough samples
			cat("Geweke ok :",resG,"(all >",limitGeweke,"), burnIn=",discard*KthinInit,"\n",file=logfile,append=TRUE);
			resESS<-ESS(sampThined[-(1:discard),])
			nbItEff<-ceiling(max(resESS))
			if(nbItEff<Nmin){
				cat("Need some more iterations, only",nbItEff,"efficient iterations, we need:",Nmin,")\n",file=logfile,append=TRUE);
				## not yet enough iterations to get the precision we need
				newNbIt<-as.integer(ceiling(Nmin*nbItBrut/nbItEff));
			}else{
				cat("All test passed, enough iterations (",nbItEff," it. eff. >",Nmin,")\n",file=logfile,append=TRUE);
				enoughIt<-TRUE
				newNbIt<-nbItBrut
			}
		}
	}else{
		cat("Raftery negative (",nbItBrut,"<",nbItMin,")\n",file=logfile,append=TRUE);
		newNbIt<-nbItMin;
	}
	output<-as.data.frame(t(c(enoughIt,newNbIt,nbItEff,KthinInit*max(discard,burnIn),Kthin,Kind,min(resG))))
	names(output)<-c("ok","newNbIt","nbItEff","burnIn","Kthin","Kind","min(resG)")

	return(invisible(new("cb.diag.out",output,colAnalysed=colAnalysed,geweke=resG,limitGeweke=limitGeweke,raftery=resRaf)));
}
# Ex: cb.diag(a<-matrix(rnorm(10000),ncol=10))


# Kernel, NB: specific treatment of spam is 2 to 10 times more efficient
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
# # presentation of the kernels:
# xabs<-seq(0,threshold)
# plot(xabs,expKernel(1,xabs,fprior),type="l",lty=1,xlab="distance (m)",ylab="strength of the association")
# lines(xabs,gaussianKernel(1,xabs,fprior),type="l",lty=2)
# lines(xabs,cauchyKernel(1,xabs,fprior),type="l",lty=3)
# lines(xabs,geometricKernel(1,xabs,fprior),type="l",lty=4)
# limf<-qlnorm(c(0.025,0.975),1+log(fprior),sdlf)
# lines(xabs,expKernel(1,xabs,limf[1]),type="l",lty=1,col="grey")
# lines(xabs,gaussianKernel(1,xabs,limf[1]),type="l",lty=2,col="grey")
# lines(xabs,cauchyKernel(1,xabs,limf[1]),type="l",lty=3,col="grey")
# lines(xabs,geometricKernel(1,xabs,limf[1]),type="l",lty=4,col="grey")
# lines(xabs,expKernel(1,xabs,limf[2]),type="l",lty=1,col="grey")
# lines(xabs,gaussianKernel(1,xabs,limf[2]),type="l",lty=2,col="grey")
# lines(xabs,cauchyKernel(1,xabs,limf[2]),type="l",lty=3,col="grey")
# lines(xabs,geometricKernel(1,xabs,limf[2]),type="l",lty=4,col="grey")
# 
# dev.print(device=pdf,paste("kernels_comparison_f",fprior,"_prior.pdf",sep=""))
# 
# legend("topright",c("Exponential","Gaussian","Cauchy","Geometric"),lty=c(1,2,3,4))

# from a natural mean and a standard deviation 
# get the mean of a corresponding log normal distribution
meansd2meansdlognorm <-function(mean=mean,sd=NA,sdlog=NA){
	if(is.na(sdlog)){
		sdlog<-sqrt(log((sd/mean)^2+1))
	}
	meanlog<-log(mean)-sdlog^2/2
	return(list(meanlog=meanlog,sdlog=sdlog))
}
meansdlognorm2meansd <-function(meanlog=meanlog,sdlog=NA){
	mean<-exp(meanlog+sdlog^2/2)
	sd<-sqrt(exp(2*meanlog+sdlog^2)*(exp(sdlog^2)-1))
	return(list(mean=mean,sd=sd))
}
# # testing
# fprior<-100
# sdfprior<-50
# flnparam<-meansd2meansdlognorm(mean=fprior,sd=sdfprior)
# 
# # mf<-log(fprior)-sdlf^2/2;
# mf<-flnparam$meanlog
# sdlf<-flnparam$sdlog
# fparam<-meansdlognorm2meansd(meanlog=mf,sdlog=sdlf)
# 
# cat("fprior",fprior,"back",fparam$mean,"sdfprior:",sdfprior,"back:",fparam$sd)

makeRuv <- function(dim,Q,K) {
  Ku <- K[1];
  Kv <- K[2];
  R <- Ku*Q;
  diag.spam(R) <- diag.spam(R) + Kv;
  R <- cbind(R, diag.spam(-1*Kv,dim,dim));
  R <- rbind(R, cbind(diag.spam(-1*Kv,dim,dim), diag.spam((Kv+1), dim, dim)));
  return(R);
}

QfromfT<-function(Dmat,SB,f=f.r,T=T.r,addeps=epsilon,origin=NULL,kern=NULL,dimension=dim(Dmat)[1]){
	if(is.spam(SB)){
		Q<- SpecificMultiply.spam(1/T,-kern(T,Dmat,f),SB);#the precision matrix of u 
	}else{
		if(SB==1){
			Q<- -kern(1,Dmat,f)
		}else{
			stop("Bad SB in QfromfT\n")
		}
	}
	diag(Q)<-0 # diagonal starts at 0

	rsumQ <- -Q %*% rep(1,dimension); # the precision for each yi
	rsumQ<-rsumQ+addeps # NB we want addeps to be in the normalisation so that we don't have too heavy outliers, addeps is the variance of isolated houses and logically is included in the normalization of the isolation component
	diag(Q) <- rsumQ;

	return(Q);
}
lik.f<-function(f,mf,sdlf,log=TRUE){
	LLH<-dlnorm(f,mf,sdlf,log=log)
	# LLH<-0*f # flat prior
	return(LLH);
}

acceptf<-{};
sample_f <- function(u,Ku,T,logsdfprop,f,mf,sdlf,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL,kern=kern){
        # sampling not symmetrically (but simple and correct) arround T
    f_prop<-rlnorm(1,log(f),logsdfprop);
    hasting_term=dlnorm(f,log(f_prop),logsdfprop,log=TRUE)-dlnorm(f_prop,log(f),logsdfprop,log=TRUE);
    
	Qprop<-QfromfT(Dmat,SB,f=f_prop,T=T,kern=kern);

	cholQprop<-get.cholMat(Qprop,cholQ)
	LLHuprop<-fast.llh.ugivQ(u,Qprop,Ku,cholQ=cholQprop);

	llhfprop<-lik.f(f_prop,mf,sdlf);
	llhf<-lik.f(f,mf,sdlf);

	LLHproposal <- LLHuprop+llhfprop; 
	LLH <-LLHu+llhf;
	lnr <- LLHproposal-LLH+hasting_term;

	cat("f:",f," LLH:",LLH,"(",LLHu,"+",llhf,") f prop:",f_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhfprop,") h_t:",hasting_term," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		f <- f_prop;
		Q<-Qprop;
		cholQ<-cholQprop;
		LLHu<-LLHuprop;
		acceptf<<-c(acceptf,1);
		cat(" accept 1\n");
	}else{
		acceptf<<-c(acceptf,0);
		cat("accept 0\n");
	}
	return(list(f=f,Q=Q,LLHu=LLHu,cholQ=cholQ));
}
lik.T<-function(T,mT,sdlT,log=TRUE){
	LLH<-dlnorm(T,mT,sdlT,log=log)
	# LLH<-0*f # flat prior
	return(LLH);
}

acceptT<-{};
sample_T <- function(u,Ku,f,T,logsdTprop,mT,sdT,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL,kern=kern){
  # sampling not symmetrically (but simple and correct) arround T	
  T_prop<-rlnorm(1,log(T),logsdTprop);
  hasting_term=dlnorm(T,log(T_prop),logsdTprop,log=TRUE)-dlnorm(T_prop,log(T),logsdTprop,log=TRUE);

  # calculate the LLH (common for both samplings)
  Qprop<-QfromfT(Dmat,SB,f=f,T=T_prop,kern=kern);

  cholQprop<-get.cholMat(Qprop,cholQ)
  LLHuprop<-fast.llh.ugivQ(u,Qprop,Ku,cholQ=cholQprop);
  llhTprop<-lik.T(T_prop,mT,sdlT);
  llhT<-lik.T(T,mT,sdlT);
  LLHproposal <- LLHuprop+llhTprop; 
  LLH <- 	       LLHu+llhT;
  lnr <- LLHproposal-LLH+hasting_term;

  cat("T:",T," LLH:",LLH,"(",LLHu,"+",llhT,") T prop:",T_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhTprop,") h_t:",hasting_term," lnr",lnr,sep="");

  if(lnr>=log(runif(1))) {
    T <- T_prop;
    Q<-Qprop;
    cholQ<-cholQprop;
    LLHu<-LLHuprop;
    acceptT<<-c(acceptT,1);
    cat(" accept 1\n");
  }else{
    acceptT<<-c(acceptT,0);
    cat("accept 0\n");
  }
  return(list(T=T,Q=Q,LLHu=LLHu,cholQ=cholQ));
}

get.cholMat<-function(Mat,cholMat=NULL){ # compute cholMat, update if possible 
	if(is.null(cholMat)){
		cholMat<-chol.spam(Mat);
	}else{
		out<-try(cholMat<-update.spam.chol.NgPeyton(cholMat,Mat));
		if(class(out)=="try-error"){ 
			cholMatbad<<-cholMat;
			Matbad<<-Mat;
			dump(list("cholMatbad","Matbad"),file="badMat.r")
			cholMat<-chol.spam(Mat);
		}
	}
	return(cholMat)
}
# general function for likelihood of u given Q
llh.ugivQ<-function(dimension,u,Q,Ku,cholQ=NULL){
	# cholQ<-chol(Q);
	cholQ<-get.cholMat(Q,cholMat=cholQ)

	exp_part<- Ku*(t(u-mean(u))%*%Q%*%(u-mean(u))); # version with mean of u has no importance
	# exp_part<- Ku*(t(u)%*%Q%*%(u)); # prior mean is at 0
	det_part<- dimension*log(Ku)+2*determinant(cholQ)$modulus
	logpi<-dimension*log(2*pi);
	LLH= -1/2*(exp_part-det_part+logpi); 
	# cat("LLH",LLH,"exppart:",exp_part,"det_part",det_part,"logpi",logpi,"\n");
	return(LLH);
}
fast.llh.ugivQ<-function(u,Q,Ku,cholQ,uQu){ # when cholQ is for sure known
  	dimension<-length(u)
	exp_part<- Ku*(t(u-mean(u))%*%Q%*%(u-mean(u))); # version with mean of u has no importance
	# exp_part<- Ku*(t(u)%*%Q%*%(u)); # prior mean is at 0
	det_part<- dimension*log(Ku)+2*determinant(cholQ)$modulus
	logpi<-dimension*log(2*pi);
	LLH= -1/2*(exp_part-det_part+logpi); 
	return(LLH);
}
llh.vgivKv<-function(dimension,v,Kv){
	LLH<-sum(dnorm(v,mean=0,sd=sqrt(1/Kv),log=TRUE))

	return(LLH);
}
llh.wgivc<-function(w,c,Kv){
  LLH<-sum(dnorm(c.compt,mean=wt,sd=sqrt(1/TKv),log=TRUE))
}

llh.ygivw<-function(y,w){
	## return the loglikelihood of y given w
	logpi<-length(y)*log(2*pi);
	LLH<- -1/2 * t(y-w)%*%(y-w) - logpi;
	return(LLH);
}
llh.zgivy<-function(y,zpos,zneg,bivect){
	probypos<-(y>0)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}
llh.zgivw<-function(w,zpos,zneg,bivect){
	probypos<-1-pnorm(0,w,1)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}
rmvnorm.prec.pseudo <- function (n, mu = rep(0, nrow(Q)), Q, Rstruct = NULL,given_rnorm=NULL, ...) {
    if (is(Rstruct, "spam.chol.NgPeyton")) 
        R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
    else R <- chol(Q, ...)
    N <- dim(Q)[1]
    if(is.null(given_rnorm)){
	    given_rnorm<-rnorm(n * N)
    }
    nu <- backsolve(R, array(given_rnorm, c(N, n)))
    return(t(nu + mu))
}
sample_u <- function(dimension,Q,K,y,cholQ=NULL){
	R <- K[1]*Q+diag.spam(1,dimension);
	center <- y;
	u <- rmvnorm.canonical(n=1, b=center, Q=R,Rstruct=cholQ);
	
	return(drop(u));
}
# enable to sample a mixture of two truncated normal with upper limit 
# of one being the lower limit of the other and this limit being 0
# this allow specifically to sample y, the "continuous reality" in the probit model
# directly according to the data

# to pick something according to the previous curve
sample_composite_ptnorm_vect <- function(xvect,bivect){
	# sample y given that it's density is a normalized sum of 
	# dnorm(xvect,1,0,+Inf)*(1-beta)+dnorm(xvect,1,-Inf,0)
	# xvect: probability of being 1 vs. 0
	# bivect: probability of being observed as 1 vs. 0
	l<-length(xvect);
	# cat("l",l);
	tunifvect<-runif(l);
	A<-pnorm(0,mean=xvect,sd=1); 	# P(y-,z-|w)=P(y-|w) as P(y-,z+)=0
	B<-(1-A)*(1-bivect);		# P(y+,z-|w)=P(y+,obs-|w)=P(obs-|y+,w)*P(y+|w)=(1-beta)*(1-P(y-|w))
	area<-A+B			# P(z-)=P(y-|w)+P(y+,obs-|w)
	samp<-B/area;			# P(y+|z-,w)=P(y+,z-|w)/P(z-)
	samp[area==0]<-1 ; 		# avoid errors when bivect=1 and xvect>0

	tfinal<-mat.or.vec(l,1);
	ypos<-which(tunifvect<=samp);
	yneg<-which(tunifvect>samp);
	# tfinal[ypos]<- rtnorm(length(ypos),mean=xvect[ypos],sd=1,lower=0,upper=Inf); # yprime=1, determine a corresponding y
	# tfinal[yneg]<- rtnorm(length(yneg),mean=xvect[yneg],sd=1,lower=-Inf,upper=0); # yprime=0, determine a corresponpding y
	# cat("\n",length(tfinal[tunifvect<=samp]),length(which(tunifvect<=samp)),length(xvect[which(tunifvect<=samp)]));
	
	if(length(ypos)>0){
		tfinal[ypos]<- rtruncnorm(1,mean=xvect[ypos],sd=1,a=0,b=Inf); # yprime=1, determine a corresponding y
	}

	if(length(yneg)>0){
		tfinal[yneg]<- rtruncnorm(1,mean=xvect[yneg],sd=1,a=-Inf,b=0); # yprime=0, determine a corresponpding y
	}
	
	return(tfinal)
}

# from there we can sample y which is sampled differently according to the observation or not of bugs in the data
sample_y_direct <-function(w,zpos,zneg,zNA,bivect){
	# return the continuous variable result of the probit model
	# it is the augmented variable of yprime
	# y prime can simply be obtained by using
	# yprime<- (y>0)
	# w is the general risk in the probit model
	# zpos the vector of references of positive z
	# zneg --------------------------- negative z
	# zNA ---------------------------- unknown z
	# bivect the vector of per house probability to find infestation if there is

	y<-0*w
	if(length(zpos)>0){
		y[zpos]<- rtruncnorm(1,mean=w[zpos],sd=1,a=0,b=Inf);
	}

	if(length(zneg)>0){
		y[zneg]<- sample_composite_ptnorm_vect(w[zneg],bivect[zneg]);
	}

	if(length(zNA)>0){
		y[zNA] <- rnorm(length(zNA),mean=w[zNA],sd=1) # w[zNA];
	}

	return(y);
}

dmtnorm<-function(x,mu=0,std=1,shift=0,w1=1,w2=1,logout=F){
	# x the value(s) to examin
	# shift the values at wich we shift from the first to the second trunc normal
	# w1/w2 the weight of the left/right part 
	# output the probability as a log
	A<-pnorm(shift,mean=mu,sd=std)*w1;
	B<-(1-pnorm(shift,mean=mu,sd=std))*w2;
	x1<-x;
	x1[x>=shift]<-Inf;
	x2<-x;
	x2[x<shift]<--Inf;
	px=dnorm(x1,mean=mu,sd=std)*w1/(A+B)+dnorm(x2,mean=mu,sd=std)*w2/(A+B);
	# cat("A",A,"B",B,"px",px,"\n");
	if(logout==T){
		px<-log(px);
	}
	return(px);
}
area.perso<-function(fun,lower,upper,eps=1e-05,visu=FALSE,...){
	# numerical integral of fun between lower and upper
	x<-seq(lower,upper,eps);
	dens<-fun(x,...);
	a<-(sum(dens)+(dens[2]-dens[length(dens)])/2)*eps;
	attributes(a)$mean<-sum(dens*x)*eps/(a);
	if(visu){
	plot(x,dens)
	}

	return(a);
}

# # testing:
# bi= 0.7 # find rate of inspectors
# x=-0.3
# A<-pnorm(0,mean=x,sd=1);
# B<-(1-pnorm(0,mean=x,sd=1))*(1-bi);
# xabs<-seq(-5,5,0.01)
# xabs1<-seq(-5,0,0.01)
# xabs2<-seq(0,5,0.01)
# 
# par(mfrow=c(1,6))
# plot(c(xabs1,xabs2),dnorm(c(xabs1,xabs2),mean=x,sd=1),type="n")
# lines(xabs1,dnorm(xabs1,mean=x,sd=1))
# lines(xabs2,dnorm(xabs2,mean=x,sd=1)*(1-bi))
# 
# plot(xabs,dmtnorm(xabs,mu=x,std=1,shift=0,w2=(1-bi)))
# 
# plot(c(xabs1,xabs2),pnorm(c(xabs1,xabs2),mean=x,sd=1),type="n")
# lines(xabs1,pnorm(xabs1,mean=x,sd=1))
# lines(xabs2,pnorm(xabs2,mean=x,sd=1)*(1-bi))
# 
# plot(c(xabs1,xabs2),pnorm(c(xabs1,xabs2),mean=x,sd=1),type="n")
# lines(xabs1,pnorm(xabs1,mean=x,sd=1)/(A+B))
# lines(xabs2,(pnorm(xabs2,mean=x,sd=1)*(1-bi)+pnorm(0,mean=x,sd=1)*bi)/(A+B))
# 
# 
# xvect<-rep(x,10000)
# bivect<-rep(bi,10000)
# system.time(y<-sample_composite_ptnorm_vect(xvect,bivect))
# 
# hist(y,freq=FALSE)
# lines(xabs1,dnorm(xabs1,mean=x,sd=1)/(A+B))
# lines(xabs2,dnorm(xabs2,mean=x,sd=1)*(1-bi)/(A+B))
# 
# # very fast (less than 2times rtnorm for 10000) and exact

# test sample_y_direct

# nbItem<-10000
# w<-rnorm(nbItem,0,3)
# bivect<-rep(0.99,nbItem)
# zpos<-sample(1:nbItem,round(20*nbItem/100))
# zNoPos<-(1:nbItem)[-zpos]
# zNApos<-sample(1:length(zNoPos),round(20*nbItem/100))
# zNA<-zNoPos[zNApos]
# zneg<-zNoPos[-zNApos]
# 
# y<-sample_y_direct(w,zpos,zneg,zNA,bivect)
# par(mfrow=c(1,3))
# plot(y[zpos],w[zpos])
# abline(a=0,b=1)
# plot(y[zNA],w[zNA])
# abline(a=0,b=1)
# plot(y[zneg],w[zneg])
# abline(a=0,b=1)

sampleK <- function(dim,Q,x,K.hyper) {
  Ku.a <- K.hyper[1];
  Ku.b <- K.hyper[2];
  Kv.a <- K.hyper[3];
  Kv.b <- K.hyper[4];
  u <- x[1:dim];
  v <- x[(dim + (1:dim))] - u;
  u.pshape <- (0.5*(dim-1) + Ku.a);
  # u.pscale <- (0.5*as.numeric(u %*% (Q %*% u)) + Ku.b^(-1))^(-1); # the bad thing you can do against convergence
  u.pscale <- (0.5*as.numeric(u-mean(u)) %*% (Q %*% (u-mean(u))) + Ku.b^(-1))^(-1); # centered u
  v.pshape <- (0.5*dim + Kv.a);
  v.pscale <- (0.5*as.numeric(v %*% v) + Kv.b^(-1))^(-1);
  Ku <- rgamma(n=1, shape=u.pshape, scale=u.pscale);
  Kv <- rgamma(n=1, shape=v.pshape, scale=v.pscale);
  K <- c(Ku,Kv);
  return(K);
}
sampleKv <- function(v,Kv.a,Kv.b){
  dim<-length(v)
  v.pshape <- (0.5*dim + Kv.a);
  v.pscale <- (0.5*as.numeric(v %*% v) + Kv.b^(-1))^(-1);
  Kv <- rgamma(n=1, shape=v.pshape, scale=v.pscale);
  return(Kv);
}
# # test sampleKv
# ntest<-10000
# Kv.r<-10
# v.r<-rnorm(500,mean=0,sd=sqrt(1/Kv.r))
# sampledKv<-rep(0,ntest)
# Kvshape <- 0.001; Kvscale <- 1000; # same for Kv
# for(i in 1: ntest){
#   Kv<-sampleKv(v.r,0*v.r,Kvshape,Kvscale)
# sampledKv[i]<-Kv
# }
# par(mfrow=c(1,2))
# plot(sampledKv)
# abline(h=Kv.r,col="green")
# abline(h=mean(sampledKv[-(1:length(sampledKv)/2)]),col="blue")
# abline(h=1/sd(v.r)^2,col=1)
# get.estimate(sampledKv[-(1:length(sampledKv)/2)],true.val=Kv.r)



acceptKv={}
sampleKvMHunifSigma<-function(Kv,v,logsdDraw){
	sig<-sqrt(1/Kv)

	sigProp<-rlnorm(1,mean=log(sig),sd=logsdDraw)
	
	LLHproposal<-sum(dnorm(v,mean=0,sd=sigProp,log=TRUE))
	LLH<-sum(dnorm(v,mean=0,sd=sig,log=TRUE))

	hasting_term=dlnorm(sig,log(sigProp),logsdDraw,log=TRUE)-dlnorm(sigProp,log(sig),logsdDraw,log=TRUE);

	lnr <- LLHproposal-LLH+hasting_term;

	cat("sigv:",sig," LLH:",LLH," sigv prop:",sigProp," LLHproposal:",LLHproposal," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		Kv <- 1/(sigProp^2);
		acceptKv<<-c(acceptKv,1);
		cat(" accept 1\n");
	}else{
		acceptKv<<-c(acceptKv,0);
		cat("accept 0\n");
	}
	return(list(Kv=Kv));
}
# # test sampleKvMHunifSigma
# ntest<-10000
# Kv.r<-10
# v.r<-rnorm(500,mean=0,sd=sqrt(1/Kv.r))
# sampledKv<-rep(0,ntest)
# Kvshape <- 0.001; Kvscale <- 1000; # same for Kv
# for(i in 1: ntest){
#    Kv<-sampleKvMHunifSigma(Kv,v.r,0.05)$Kv
# sampledKv[i]<-Kv
# }
# par(mfrow=c(1,2))
# plot(sampledKv)
# abline(h=Kv.r,col="green")
# abline(h=mean(sampledKv[-(1:length(sampledKv)/2)]),col="blue")
# abline(h=1/sd(v.r)^2,col=1)
# get.estimate(sampledKv[-(1:length(sampledKv)/2)],true.val=Kv.r)



acceptKu={}
## ok for both Kv and Kc (and could be used for Ku with dim<-dimension-1)

sampleKu <- function(dim,Q,u,Kushape,Kuscale) {
	pos.shape <- (0.5*(dim-1) + Kushape);
	pos.scale <- (0.5*as.numeric(u %*% Q %*% u) + Kuscale^(-1))^(-1);
	Ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
	return(Ku);
}

llh.c<-	function(c.val,c.comp,Kc,wnoc,y,zNA){
  w<-wnoc+c.comp;
  if(length(zNA)>0){
    w<-w[-zNA]
    y<-y[-zNA]
  }
  LLH<-llh.ygivw(y,w)+sum(dnorm(c.val,mean=0,sd=sqrt(1/Kc),log=TRUE));
  return(LLH);
}
llh.cbis<-function(c.val,c.comp,Kc,wnospat,y,zNA,Kv){
  if(length(zNA)>0){
    wnospat<-wnospat[-zNA]
    c.comp<-c.comp[-zNA]
    y<-y[-zNA]
  }
  LLH<-sum(dnorm(y,mean=c.comp,sd=1+sqrt(1/Kv)),log=TRUE)+sum(dnorm(c.val,mean=0,sd=sqrt(1/Kc),log=TRUE));
  return(LLH);
}
acceptc.val<-{}
mhsamplec<-function(c.val,c.comp,c.map,sdc.val,Kc,wnoc,y,zNA){
	nbc<-length(c.val);
	c.valprop<-rnorm(nbc,c.val,sdc.val);
	# hasting_term=dnorm(c.val,log(c.valprop),logsdc.val,log=TRUE)-dnorm(c.valprop,log(c.val),logsdc.val,log=TRUE);
	c.compprop<-drop(c.map%*%c.valprop);
	LLHproposal <- llh.c(c.valprop,c.compprop,Kc,wnoc,y,zNA); 
	LLH <-llh.c(c.val,c.comp,Kc,wnoc,y,zNA);
	lnr <- LLHproposal-LLH# +hasting_term;
	# cat("c.val:",c.val,"LLH:",LLH,"c.valprop:",c.valprop,"LLH proposal:",LLHproposal,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		c.val <- c.valprop;
		c.comp<-c.compprop;
		acceptc.val<<-c(acceptc.val,1);
		# cat(" accept 1\n");
	}else{
		acceptc.val<<-c(acceptc.val,0);
		# cat("accept 0\n");
	}

	return(list(c.val,c.comp));
}
mhsamplecbis<-function(c.val,c.comp,c.map,sdc.val,Kc,wnospat,y,zNA,Kv){
	nbc<-length(c.val);
	c.valprop<-rnorm(nbc,c.val,sdc.val);
	# hasting_term=dnorm(c.val,log(c.valprop),logsdc.val,log=TRUE)-dnorm(c.valprop,log(c.val),logsdc.val,log=TRUE);
	c.compprop<-drop(c.map%*%c.valprop);
	LLHproposal <- llh.cbis(c.valprop,c.compprop,Kc,wnospat,y,zNA,Kv); 
	LLH <-llh.cbis(c.val,c.comp,Kc,wnospat,y,zNA,Kv);
	lnr <- LLHproposal-LLH# +hasting_term;
	# cat("c.val:",c.val,"LLH:",LLH,"c.valprop:",c.valprop,"LLH proposal:",LLHproposal,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		c.val <- c.valprop;
		c.comp<-c.compprop;
		acceptc.val<<-c(acceptc.val,1);
		# cat(" accept 1\n");
	}else{
		acceptc.val<<-c(acceptc.val,0);
		# cat("accept 0\n");
	}

	return(list(c.val,c.comp));
}
samplexuv <- function(dim,Q,K,y,cholR=NULL) {
  x <- rnorm(n=(2*dim), mean=0, sd=1);
  center <- c(rep(0,dim), y);
  R <- makeRuv(dim,Q,K);
  if(is.null(cholR)){
	  cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
  }else{
	  cholR <- update.spam.chol.NgPeyton(cholR,R);
  }

  center <- backsolve(cholR, forwardsolve(cholR, center));
  x <- backsolve(cholR,x);
  x <- x + center;
  cholR<<-cholR;
  return(x);
}

fastsamplexuv <- function(dim,cholR,y) {
  x <- rnorm(n=(2*dim), mean=0, sd=1);
  center <- c(rep(0,dim), y);
  center <- backsolve(cholR, forwardsolve(cholR, center));
  x <- backsolve(cholR,x);
  x <- x + center;
  return(x);
}

sample_w_nospat<-function(y,c,Kv){
## classical formula for posterior given
## y ~ N(w,1)
## w ~ N(c,1/Kv)
  var<-1/Kv
  meanPost<-y*var/(1+var)+c/(1+var)
  sdPost<-sqrt(1/(Kv+1))
  w<-rnorm(length(w),mean=meanPost,sd=sdPost)
  return(w)
}


sample_v<-function(ynotv,Kv){
  cov<-1/(1+Kv)
  meanPost<- cov * ynotv
  sdPost<- sqrt(cov)
  v<-rnorm(length(ynotv),mean=meanPost,sd=sdPost)
  return(v)
}


samplebeta <- function(zpos,zneg,matrix,yprime,a,b) {
	yp.positive <- yprime;
	yp.negative <- yprime; 
	yp.positive[-zpos] <- 0; # z- turn y+ into 0 a.pos being then the sum of the y+ well observed
	yp.negative[-zneg] <- 0; # z+ turn y+ into 0 a.neg being then the sum of the y+ badly observed 
	# the NA are not apparent here
	a.pos <- (as.vector(t(matrix) %*% yp.positive) + a); 
	b.pos <- (as.vector(t(matrix) %*% yp.negative) + b); 
	beta <- rbeta(n=ncol(matrix), shape1=a.pos, shape2=b.pos); 
	return(beta);
}
generate_z <-function(y,b,zNA){
	# y the probit continuous value
	# b the vector with per house the discovery rate of inspectors when y>0
	# if NA, z is set to NA
	b[zNA]=0;
	pre_z<-(y>0);
	zypos<-rbinom(sum(pre_z),1,b[pre_z]);
	pre_z[pre_z]<-zypos;
	pre_z[zNA]<-NA
	return(pre_z);
}
subsetAround <- function (priordatafullmap, reportsUnicodes, threshold, ...) {
  reportslines <- match(reportsUnicodes, priordatafullmap$unicode)
  if(sum(is.na(reportslines))>0){
    stop("new points not in initial map") # in fine better to send a warning
  }
  dist_mat <- nearest.dist(x = priordatafullmap[, c("X", "Y")], 
			   y = priordatafullmap[reportslines, c("X", "Y")], method = "euclidian", 
			   delta = threshold, upper = NULL)
  selectedlines <- which(!is.na(apply_by_row_not_null.spam(dist_mat, 
							   min)))
  subprior <- priordatafullmap[selectedlines, ]
  return(subprior)
}

fit.spatautocorel<-function(db=NULL,
			    cofactors=c(),
			    pfile="parameters_extrapol.r",
			    fit.spatstruct=TRUE,
			    use.generated=FALSE,
			    make.map.cofactors=FALSE,
			    kern="exp",
			    nbiterations=100,
			    nocheck=FALSE,
			    threshold=50,
			    use.v=TRUE
			    ){
  # # db should contain in columns at least:
  # X 
  # Y
  # positive
  # # optional
  # block_num the polygon number (only constraint is being the same for all items of a group on only them)
  # fit.spatstruct : do we fit the spatial precision (autocorrelation) structure (f and T)
  # use.generated : use or not generated data on the original map
  # cofactors: names of the columns corresponding to cofactors to be used in db
  # make.map.cofactors : should the vector of cofactor be used to generate pseudo values for the cofactors in spite of using the ones in db?
  # nbiterations: nb of iterations of the MCMC chain if <0 the mcmc is automatically stopped when raftery and geweke diagnostic are satisfied (can take days)
  # threshold: distance above which the spatial linked is not assessed (considered null), Nota: this can be much lower than the distance at which there is covariance as it is linked to the partial-covariance, not the general covariance
  # the function can be "softly stopped", saving everything by 
  # uncommenting break() in manual.stop.r

  # source(pfile)
  source("parameters_extrapol.r")
  
  # check fondamental variables present
  if(is.null(db$positive)){
    cat("\nMissing \"positive\" in db. Aborting.\n")
    return(NULL)
  }


  cat("\n")
  cat("Assessing computation to perform\n")
  cat("--------------------------------\n")

  # total size
  cat("Mapping:",dim(db)[1],"points\n")
  
  # cofactors
  nbfact.gen<-length(cofactors)
  if(nbfact.gen==0){
    use.cofactors<-FALSE;
    mes<-c()
  }else{
    use.cofactors<-TRUE
    cofs<-cofactors
    mes<-paste(cofs)
  }
  cat("Account for cofactors:",use.cofactors,mes,"\n");

  # Observers
  if(is.null(db$IdObserver)){
    use.insp<-FALSE
    mes<-""
    beta<-1
  }else{
    use.insp<-TRUE
    mes<-paste("(",length(levels(as.factor(db$IdObserver))),")",sep="")
  }
  cat("Account for observers:",use.insp,mes,"\n") 

  # Non observed
  if(is.null(db$observed)){
    use.NA<-FALSE
    mes<-""
  }else{
    use.NA<-TRUE
    mes<-paste("(",length(which(db$observed==0)),")",sep="")
  }
  cat("Account for non-observed points:",use.NA,mes,"\n") 

  if(!is.null(db$p.i)){
    intercept <- qnorm(mean(db$p.i))
  }else{
    intercept <- 0
  }
  dimension <- nrow(db);
  db$status<-rep(0,dim(db)[1])
  db$status[db$observed!=1]<-9
  db$status[db$positive==1]<-1
  INTERMEDIARY<-FALSE


  cat("Account for spatial autocorrelation:")
  if(!is.null(db$X)){
    use.spat<-TRUE
    cat(use.spat)
  # affectation of the kernel
  cat(" Using kernel:",kern)
  if(kern == "exp"){
    kern<-expKernel;
  }else if(kern == "gaussian"){
    kern<-gaussianKernel;
  }else if(kern == "cauchy"){
    kern<-cauchyKernel;
  }else if(kern == "geometric"){
    kern<-geometricKernel;
  }else{
    stop("Error, unknown kernel:",kern,"Please change kern")
  }
  # print(kern)


  ## functions
  # cat("generate distance matrix\n")
  spam.options(nearestdistnnz=c(9058076,400))
  dist_mat <-nearest.dist(x=db[,c("X","Y")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
  diag(dist_mat)<- rep(0,dim(dist_mat)[1])
  cat("\nAccount for known barriers: ")
  }else{
  # affectation of the kernel
    use.spat<-FALSE
    cat(use.spat)
    kern<-expKernel
    fit.spatstruct<-FALSE
    dm<-mat.or.vec(dim(db)[1],dim(db)[1])
    dist_mat<-as.spam(dm)
  }

  if(!is.null(db$GroupNum)){
    use.streets<-TRUE
    SB <- nearest.dist(x=cbind(db$GroupNum,rep(0,length(db$GroupNum))), method="euclidian", upper=NULL,delta=0.1)
    SB@entries<-rep(1,length(SB@entries))

    dmt<-dist_mat
    dmt@entries<-rep(1,length(dmt@entries))# [dmt@entries!=0]<-1 # 1 only when dist_mat not 0

    SB<-as.spam(SB*dmt);

    AS<-dmt-SB; # get 1 whereever the distances matrix is defined(under threshold) and not same block
    AS<-as.spam(AS)
    mes<-paste("(",length(levels(as.factor(db$GroupNum)))," groups)",sep="");

    # for fastest assessment of sameblock/accross street share
    ASdist<-as.spam(dist_mat*AS)
    SBdist<-as.spam(dist_mat*SB)
  }else{
    use.streets<-FALSE
    T<-1
    SB<-1
    AS<-0
    mes<-""
    meanSBshare<-NA
    meanSBshareNoT<-NA
    SBdist<-1
  }
  if(use.spat){
  cat(use.streets,mes,"\n")
  }

  cat("Account for local noise:",use.v,"\n")

  # should have here kernel used etc...

  # computation size
  if(nbiterations<0){
    final.run<-TRUE
    mes<-"Performing as many iterations as necessary for convergence."
    adaptOK<-FALSE
  }else{
    final.run<-FALSE
    mes<-paste("Stopping after",nbiterations,"iterations of the MCMC")
    adaptOK<-TRUE # No adaptation of sampling sd when fix number of iterations
  }
  cat(mes,"\n")

  if(interactive() && !nocheck){
    continue<-readline("Is this what you want?(Y/n)\n")
    if(continue=="n" || continue== "N"){
      cat("Ok aborting.\n")
      return(NULL)
    }
  }

  if(use.insp){
    data.insp <- db$IdObserver;

    inspectors <- unique(data.insp[which(db$observed==1)]);
    inspector <- matrix(0,dimension,length(inspectors))
    for (i in 1:length(inspectors)) {
      inspector[,i] <- (data.insp==inspectors[i]);
    }
    inspector <- as.spam(inspector);

    beta<-rep(1,length(inspectors))
    bivect <- as.vector(inspector %*% beta);
    # par(mfrow=c(1,1))
  }else{
    bivect<-rep(priorinspquality,length(db$X))
  }

  ## set the name
if(use.generated){racgen<-"gen"}else{racgen<-"true"}
if(use.v){racv<-"v"}else{racv<-"uv"}
if(use.insp){racinsp<-"Insp"}else{racinsp<-"NoInsp"}

## prep db 
if(use.generated){
	zpos <- which(z.r==1);
	zneg <- which(z.r==0);
	zNA <- which(z.r==9);
}else{
	z.r<-db$status
	## initialisation specific to true db
	zpos <- which(db$status==1);
	zneg <- which(db$status==0);
	zNA <- which(db$status==9);
	select<-c(zneg,zpos)
	# plot_reel(db$X[select],db$Y[select],2*db$status[select]-1,main="data")
}

# priors and sampling settings
nparam<-3
if(use.v){
	nparam=nparam+1;
}
if(use.cofactors){
	nparam=nparam+1;
}
if(use.insp){
	nparam=nparam+1;
}

## definition logmean for log normal priors
flnparam<-meansd2meansdlognorm(mean=fprior,sdlog=sdlf)
mf<- flnparam$meanlog 

Tlnparam<-meansd2meansdlognorm(mean=Tprior,sdlog=sdlT)
mT<- Tlnparam$meanlog # 1 

K.hyper <- c(Kushape,Kuscale,Kvshape,Kvscale); 

## initial sd for cofactors decrease with the number of cofactors
if(use.cofactors){
	sdc.val<-0.1/nbfact.gen;
}

## hyperparameter finding rate
b <- c(abeta,bbeta) # hyperparameters of the detection quality


if(visu.progression){
  ## priors plotting
  par(mfrow=c(1,nparam))

  xabs<-seq(0.001,5*fprior,0.1)
  plot(xabs,lik.f(xabs,mf,sdlf,log=FALSE),main="f prior",xlab="f",ylab="density",type="l")

  xabs<-seq(0,50,0.1)
  plot(xabs,lik.T(xabs,mT,sdlT,log=FALSE),main="T prior",xlab="T",ylab="density",type="l")

  xabs<-seq(0.0001,10,0.1)
  plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density",type="l")

  if(use.v){
    xabs<-seq(0.0001,10,0.1)
    plot(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),main="Kv prior",xlab="Kv",ylab="Density",type="l")
  }

  # b finding rate when positive
  ## uniform prior
  if(use.insp){
    xabs<-seq(0,1,0.01)
    plot(xabs,dbeta(xabs,abeta,bbeta),main="Betas prior",xlab="Beta",ylab="Density",type="l")
  }
}

## set up the inspector matrix for fast access to observer information
if(use.insp){
  data.insp <- db$IdObserver;
  inspectors <- unique(data.insp[which(db$status!=9)]);
	inspector <- matrix(0,dimension,length(inspectors))
	for (i in 1:length(inspectors)) {
		inspector[,i] <- (data.insp==inspectors[i]);
	}
	inspector <- as.spam(inspector);
	beta<-rep(1,inspector@dimension[2])
	bivect <- as.vector(inspector %*% beta);
}else{
	bivect<-rep(priorinspquality,dimension)
}


## starting values
K<-c(Ku,Kv,Kc);
if(use.streets){
	Q<-QfromfT(dist_mat,SB,f,T,kern=kern);
}else{
	Q<-QfromfT(dist_mat,SB,f,T=1,kern=kern);
}

# cat("at Q build\n T:",T,"f:",f,"Ku:",Ku,"Kv:",Kv,"\n")
cholQ<-chol(Q);

u<-rep(0,dimension);
y<-rep(0,dimension);
yprime <- (y>0);
if(use.v){
w <- u # rnorm(dimension,u,sqrt((K[2])^(-1)));
v<-w-u;
x <- c(u, w);
}else{
  w<-u
  v<-0*u
  x<-c(u,w)
}


R <- makeRuv(dimension,Q,K);
cholR <- chol.spam(R,memory=list(nnzcolindices=300000));

if(use.cofactors){
  c.map<-as.matrix(db[,cofs])
  c.val<-rep(0,length(cofs));
  c.comp<-c.map%*%c.val
  wnotr<-c.comp
}else{
  LLHc<-0
  wnotr<-0*u
}

ItTestNum<-gibbsit(NULL,NminOnly=TRUE);
cat("init ItTestNum:",ItTestNum,"\n",file="convergence_tests.txt")
if(final.run){
	# lastAdaptProp: end of adjustment of the sampling variances
	nbiterations<-ItTestNum;
	lastAdaptProp<-nbiterations
	use.autostop<-TRUE
}else{
	use.autostop<-FALSE;
	lastAdaptProp <- 0
}
LLHu<-llh.ugivQ(dimension,u,Q,K[1])
LLH<-llh.zgivw(w,zpos,zneg,bivect)
if(use.insp){
	LLHb<-sum(dbeta(beta,abeta,bbeta,log=TRUE));
}else{
	LLHb<-0;
}
if(use.cofactors){
	LLHc<-sum(dnorm(c.val,0,sqrt(1/Kc),log=TRUE));
	c.valsamp<-as.matrix(mat.or.vec(nbiterations+1,nbfact.gen));
}else{
	LLHc<-sum(dnorm(rep(0,nbfact.gen),mean=0,sd=sqrt(1/Kc),log=TRUE));
}
LLHv<-sum(dnorm(v,0,Kv,log=TRUE));

## save starting values
nbtraced=19;
if(fit.spatstruct){
	sampled<-as.matrix(mat.or.vec(nbiterations+1,nbtraced));
	sampled[1,1]<-T;
	sampled[1,2]<-LLHu
	sampled[1,3]<-f;
	sampled[1,4]<-LLHu
	sampled[1,5]<-Ku;
	sampled[1,6]<-llh.zgivy(y,zpos,zneg,bivect);
	sampled[1,7]<-LLH;
	sampled[1,8]<-llh.ygivw(y,w);
	sampled[1,9]<-0
	sampled[1,10]<-Kv;
	sampled[1,11]<-mean(u);
	sampled[1,12]<-Kc;
	sampled[1,13]<-LLHv;
	sampled[1,14]<-LLHc;
	sampled[1,15]<-LLHb;
	sampled[1,16]<-0;
	sampled[1,17]<-0
	sampled[1,18]<-0
	sampled[1,19]<-mean(beta);
	namesSampled<-c("T","LLHTu","f","LLHfu","Ku","LLHy","LLH","LLHyw","i","Kv","mu","Kc","LLHv","LLHc","LLHb","LLHTotal","meanSBshare","meanSBshareNoT","meanBeta")
	
}else{
	grid.stab<-seq(1,length(w),ceiling(length(w)/5))# values of the field tested for stability, keep 5 values
	nbtraced<-2*(2+length(grid.stab))+4
	spacer<-(2+length(grid.stab))

	sampled<-as.matrix(mat.or.vec(nbiterations+1,nbtraced));
	sampled[1,1]<-mean(u)
	sampled[1,2]<-sd(u)
	sampled[1,3:spacer]<-u[grid.stab]
	sampled[1,spacer+1]<-mean(w)
	sampled[1,spacer+2]<-sd(w)
	sampled[1,(spacer+3):(2*spacer)]<-w[grid.stab]
	LLHu<-llh.ugivQ(dimension,u,Q,K[1])
	sampled[1,(2*spacer)+1]<-llh.ugivQ(dimension,u,Q,K[1])
	sampled[1,(2*spacer)+2]<-llh.ygivw(y,w);
	sampled[1,(2*spacer)+3]<-llh.zgivy(y,zpos,zneg,bivect);
	sampled[1,(2*spacer)+4]<-llh.zgivw(w,zpos,zneg,bivect);
	namesSampled<-c("meanu","sdu",rep("localu",length(grid.stab)),"meanw","sdw",rep("localw",length(grid.stab)),"llhugivQ","llhygivw","llhzgivy","llhzgivw")
}
if(length(namesSampled)!=nbtraced){
  stop("Number of names in sampled != number of columns")
}

write.table(t(namesSampled),monitorfile,sep="\t",col.names=FALSE,row.names=FALSE)
write.table(t(sampled[1,]), monitorfile, sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
lastsaved<-1

write.table(t(u), "usamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
write.table(t(w), "wsamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
write.table(t(y), "ypsamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
if(use.cofactors){
	write.table(t(cofs), coffile, sep="\t",col.names=FALSE,row.names=FALSE,append=FALSE)
	write.table(t(c.val), coffile, sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
}else{
  if(file.exists(coffile)){
    file.remove(coffile)
  }
}
if(use.insp){
	write.table(t(beta),betafile, sep="\t",col.names=FALSE,row.names=FALSE)
}else{
  if(file.exists(betafile)){
    file.remove(betafile)
  }
}


## Display Data
z<-db$positive
z[db$observed==0]<-0.5

if(use.spat){
dev.new()
par(mfrow=c(1,2))
plot_reel(db$X,db$Y,z,base=0,top=1,main="data with pale Non Observed")
legend("bottomright",c("Negative","Non Observed","Positive"),col=c("black","yellow4","yellow"),pch=15)
plot.classes(db$X,db$Y,db$GroupNum,main="Groups")
}

## initialization of small recording variables
nbploted<-2
sum.u<-rep(0,dimension);
if(use.cofactors){
	nbploted<-nbploted+1
	sum.c.val<-0*c.val;
}
if(use.v){
	nbploted<-nbploted+1
	sum.v<-rep(0,dimension);
}
sum.w<-rep(0,dimension);
sum.y<-rep(0,dimension);
sum.yp<-rep(0,dimension)
if(use.insp){
	nbploted<-nbploted+1
	sum.beta<-beta*0;
}
if(visu.progression){
	dev.new()
	par(mfcol=c(nbploted,5))
}
sumSBshares<-rep(0,length(w))

### Main Loop
mainLoopStart<-proc.time()
num.simul<-1;
Kthin<-1
lastBurnIn<-0
# set.seed(777)
# cat("rand",rnorm(1));
cat("Begin sampling...\n")
while (num.simul <= nbiterations || (!adaptOK && final.run)) {
  # cat("\n\n",num.simul,", out of:",nbiterations," ");
  # cat("c.val",c.val,"mw",mean(w),"Kv",Kv,"my",mean(y),"\n")

    # update y (latent infestation variable)
    y<-sample_y_direct(w,zpos,zneg,zNA,bivect);
    yprime <- (y>0);
    # cat("meany:",mean(y),"sums zpos",sum(zpos),"neg",sum(zneg),"NA",sum(zNA),"b",mean(bivect))

    # update "r" the spatial component and/or non-spatial noise 
    if(use.spat){ # spatial component
      if(use.v){ # local error
	if(fit.spatstruct){
	  x <- samplexuv(dimension,Q,K,y-wnotr,cholR);

	  K <- sampleK(dimension,Q,x,K.hyper);
	  Ku<-K[[1]]
	  Kv<-K[[2]];
	  cat("Ku",Ku,"Kv",Kv,"\n");
	}else{
	  x <- fastsamplexuv(dimension,cholR,y-wnotr);
	}
	u<-x[1:dimension];
	v<-x[dimension+(1:dimension)]-u;
	wnoc<-x[dimension+(1:dimension)]
      }else{ # no local error
	u<-sample_u(dimension,Q,K,y-wnotr,cholQ);
	wnoc<-u
	if(fit.spatstruct){
	  Ku<-sampleKu(dimension,Q,u,Kushape,Kuscale); 
	  K[1]<-Ku;
	  cat("Ku",Ku,"\n");
	}
      }
    }else{ # no spatial component
      if(use.v){ # local error
	u<-0*w
	v<-sample_v(y-wnotr,Kv)
	Kv<-sampleKv(v,Kvshape,Kvscale)
	wnoc<-v
      }else{ # no local error
	wnoc<-0*w
      }
    }

    if(use.cofactors){
      c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,wnoc,y,zNA);
      c.val<-c.all[[1]]
      c.comp<-c.all[[2]]
      wnotr<-c.comp
	# cat("mc",mean(c.all),"sdc.val",sdc.val,"mv",mean(y),"mwnoc",mean(wnoc),"Kc",Kc,"\n");
    }else{
      wnotr<-0*w
    }
    w<-wnoc+wnotr

  if(fit.spatstruct){
    LLHu<-llh.ugivQ(dimension,u,Q,Ku,cholQ=cholQ) 
  # LLHu<-fast.llh.ugivQ(u,Q,Ku,cholQ=cholQ) # for some mysterious reason, cannot be used here, some update of cholQ is not done correctly
  }

  ## update the autocorrelation kernel
  if(fit.spatstruct){
      out<-sample_f(u,Ku,T,logsdfprop,f,mf,sdlf,Q,LLHu,AS,SB,cholQ=cholQ,Dmat=dist_mat,kern=kern); 
      f<-out$f;
      Q<-out$Q;
      LLHu<-out$LLHu;
      cholQ<-out$cholQ;

    if(use.streets==TRUE){
      out<-sample_T(u,Ku,f,T,logsdTprop,mT,sdT,Q,LLHu,AS,SB,cholQ=cholQ,Dmat=dist_mat,kern=kern);
      T<-out$T;
      Q<-out$Q;
      LLHu<-out$LLHu;
      cholQ<-out$cholQ;
    }
  }

  if(use.insp){
    beta <- samplebeta(zpos,zneg,inspector,yprime,abeta,bbeta);
    bivect <- as.vector(inspector %*% beta);
    cat("beta (",mean(beta),",",sd(beta),") suminsp",sum(inspector),"sumyp",sum(yprime));
  }

  if(lastAdaptProp<num.simul){
    sum.u<-sum.u+u;
    if(use.v){
      sum.v<-sum.v+v;
    }
    if(use.cofactors){
      sum.c.val<-sum.c.val+c.val;
    }
    sum.w<-sum.w+w;
    sum.y<-sum.y+y;
    if(use.insp){
      sum.beta<-sum.beta+beta;
    }
    yp<-as.integer(y>0)
    sum.yp<-sum.yp+yp;
  }
  if(use.v){
    LLHv<-sum(dnorm(v,0,Kv,log=TRUE));
  }
  if(use.cofactors){
    LLHc<-sum(dnorm(c.val,0,Kc,log=TRUE));
  }
  if(use.insp){
    LLHb<-sum(dbeta(beta,abeta,bbeta,log=TRUE));
  }

  # adapt sampling
  if((num.simul)%% 20 ==0 && !adaptOK){
    adaptOK<-TRUE
    if(fit.spatstruct){
      # adapt sampling of T
      if(use.streets){
	logsdTprop<-adaptSDProp(logsdTprop,acceptT)
	adaptOK<-adaptOK && attributes(logsdTprop)$noupdate
      }
      # adapt sampling of f
      logsdfprop<-adaptSDProp(logsdfprop,acceptf)
      adaptOK<-adaptOK && attributes(logsdfprop)$noupdate
    }
    # adapt sampling of c.val
    if(use.cofactors){
      sdc.val<-adaptSDProp(sdc.val,acceptc.val)
      adaptOK<-adaptOK && attributes(sdc.val)$noupdate
    }
    if(adaptOK){
      cat("\n Adaptation of sampling variances OK\n");
      if(final.run){
	lastAdaptProp<-num.simul;
	nbiterations<-lastAdaptProp+ItTestNum;
	sampled<-resized(sampled,nr=nbiterations+1);
	if(use.cofactors){
	  c.valsamp<-resized(c.valsamp,nr=nbiterations+1)
	}
      }
    }else{
      # adapt sampled size if needed
      if(nbiterations<num.simul+1 && final.run){
	nbiterations<-num.simul+ItTestNum
	sampled<-resized(sampled,nr=nbiterations+1)
	if(use.cofactors){
	  c.valsamp<-resized(c.valsamp,nr=nbiterations+1)
	}
      }
    }
    cat(file="convergence_tests.txt","num.simul:",num.simul,"AdaptOK:",adaptOK,append=TRUE)
    cat(file="convergence_tests.txt","lastAdaptProp",lastAdaptProp,"nbiterations",nbiterations,"\n",append=TRUE)
  }
  LLHyw<-llh.ygivw(y,w);
  LLHy<-llh.zgivy(y,zpos,zneg,bivect);
  LLH<-llh.zgivw(w,zpos,zneg,bivect);
  LLHTotal<-LLH;
  llhf<-lik.f(f,mf,sdlf);
  llhT<-lik.T(T,mT,sdlT);
  llhKu<-dgamma(Ku,shape=Kushape,scale=Kuscale,log=TRUE)
  llhKv<-dgamma(Kv,shape=Kvshape,scale=Kvscale,log=TRUE)

  LLHTotal<-LLH+LLHu+LLHv+LLHc+llhf+llhT+llhKu+llhKv+LLHb

  # cat("LLHtotal:",LLHTotal,"(",LLH,"+",LLHu,"+",LLHv,"+",LLHc,"+",llhf,"+",llhT,"+",llhKu,"+",llhKv,"+",LLHb,")\n");
  # cat("mu:",mean(u),"sdu:",sd(u),"T:",T,"f:",f);

  # Same block weight in spatial component
  if(use.streets){
  kernSB<-kern(1,SBdist,f)
  kernAS<-kern(T,ASdist,f)
  kernASnoT<-kern(1,ASdist,f)
  SBweights<-apply_by_row_not_null.spam(kernSB,sum,na.rm=TRUE)
  ASweights<-apply_by_row_not_null.spam(kernAS,sum,na.rm=TRUE)
  ASweightsNoT<-apply_by_row_not_null.spam(kernASnoT,sum,na.rm=TRUE)
  SBshares<-SBweights/(ASweights+SBweights)
  SBsharesNoT<-SBweights/(ASweightsNoT+SBweights)
  meanSBshare<-mean(SBshares,na.rm=TRUE)
  meanSBshareNoT<-mean(SBsharesNoT,na.rm=TRUE)
  sumSBshares<-sumSBshares+SBshares
  cat("mean SB share:",meanSBshare,"w/o barriers would be:",meanSBshareNoT)
  }else{
  }

  ## monitored variables
  if(fit.spatstruct){
    sampled[num.simul+1,1]<-T;
    sampled[num.simul+1,2]<-LLHu;
    sampled[num.simul+1,3]<-f;
    sampled[num.simul+1,4]<-LLHu;
    sampled[num.simul+1,5]<-K[1];
    sampled[num.simul+1,6]<-LLHy;
    sampled[num.simul+1,7]<-LLH;
    sampled[num.simul+1,8]<-LLHyw;
    sampled[num.simul+1,9]<-num.simul;
    sampled[num.simul+1,10]<-K[2];
    sampled[num.simul+1,11]<-mean(u);
    sampled[num.simul+1,12]<-Kc;
    sampled[num.simul+1,13]<-LLHv;
    sampled[num.simul+1,14]<-LLHc;
    sampled[num.simul+1,15]<-LLHb;
    sampled[num.simul+1,16]<-LLHTotal;
    sampled[num.simul+1,17]<-meanSBshare
    sampled[num.simul+1,18]<-meanSBshareNoT
    sampled[num.simul+1,19]<-mean(beta);
  }else{
    sampled[num.simul+1,1]<-mean(u)
    sampled[num.simul+1,2]<-sd(u)
    sampled[num.simul+1,3:spacer]<-u[grid.stab]
    sampled[num.simul+1,spacer+1]<-mean(w)
    sampled[num.simul+1,spacer+2]<-sd(u)
    sampled[num.simul+1,(spacer+3):(2*spacer)]<-w[grid.stab]
    sampled[num.simul+1,(2*spacer)+1]<-LLHu
    sampled[num.simul+1,(2*spacer)+2]<-LLHyw;
    sampled[num.simul+1,(2*spacer)+3]<-LLHy;
    sampled[num.simul+1,(2*spacer)+4]<-LLH;
  }
  if(use.cofactors){
    c.valsamp[num.simul+1,]<-c.val
  }
  ## auto stopping
  if((num.simul==ItTestNum+lastAdaptProp && num.simul>lastAdaptProp +1 && final.run)){
    # accommodate for limited length acceptable by gibbsit
    KthinAcc<-ceiling(num.simul/45000)
    cb<-cb.diag(sampled[(1+lastAdaptProp):num.simul,-9],logfile="convergence_tests.txt",KthinInit=KthinAcc);
    Kthin<-cb$Kthin*KthinAcc
    lastBurnIn<-cb$burnIn+lastAdaptProp
    ItTestNum<-min(cb$newNbIt,num.simul*3);
    if(use.autostop){
      nbiterations<-ItTestNum+lastAdaptProp;
      if(!cb$ok){
	sampled<-resized(sampled,nr=nbiterations+1);
	if(use.cofactors){
	  c.valsamp<-resized(c.valsamp,nr=nbiterations+1)
	}
      }
    }

    cat("cb: num.simul \t begEst \t ok \t newNbIt \t new nbiterations \t nbItEff \t burnIn \t Kthin \t Kind \t min(resG)\n",file="convergence_tests.txt",append=TRUE)
    cat("\t",num.simul,"\t",lastAdaptProp,"\t",cb$ok,"\t",cb$newNbIt,"\t",nbiterations,"\t",cb$nbItEff,"\t",cb$burnIn,"\t",cb$Kthin,"\t",cb$Kind,"\t",cb$"min(resG)","\n",file="convergence_tests.txt",append=TRUE)
  }

  ### save to files every freqsave or before stopping
  if(num.simul%%freqsave==0 || num.simul==(nbiterations)){ 
    cat("\nSaving at",num.simul,", out of:",nbiterations," ");
    write.table(t(u), "usamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(w), "wsamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(as.numeric(yprime)), "ypsamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)

    toWrite<-sampled[(lastsaved+1):(num.simul+1),]
    if(is.vector(toWrite)){
      toWrite<-t(toWrite)
    }
    write.table(toWrite, monitorfile,append=TRUE, sep="\t",col.names=FALSE,row.names=FALSE)
    if(visu.progression){
      plot_reel(db$X,db$Y,y,main="y")
      plot_reel(db$X,db$Y,u,main="u")
    }
    if(use.cofactors){
      if(visu.progression){
	  plot(c.val)
      }
      toWrite<-c.valsamp[(lastsaved+1):(num.simul+1),]
    if(is.vector(toWrite)){
      toWrite<-t(toWrite)
    }
    write.table(toWrite, coffile, sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    }
    if(use.v){
      v<-x[dimension+(1:dimension)]-u;
      if(visu.progression){
	plot_reel(db$X,db$Y,v,main="v")
      }
    }
    if(use.insp){
      write.table(t(beta), betafile, sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      if(visu.progression){
	plot_reel(db$X,db$Y,bivect*2-1,main="insp")
      }
    }

    lastsaved<-num.simul+1
    ## manual stop
    try(source("manual.stop.r",local=TRUE))
  }
  
  num.simul<-num.simul+1
}

mainLoopStop<-proc.time()
cat("Time main loop:",mainLoopStop-mainLoopStart,"\n")

nbiterations<-(num.simul-1)
est.u<-sum.u/(nbiterations-lastAdaptProp);
dump("est.u",file=paste("estimated.txt",sep=""),append=TRUE)
est.y<-sum.y/(nbiterations-lastAdaptProp);
dump("est.y",file=paste("estimated.txt",sep=""),append=TRUE)
est.w<-sum.w/(nbiterations-lastAdaptProp);
dump("est.w",file=paste("estimated.txt",sep=""),append=TRUE)
est.yp<-sum.yp/(nbiterations-lastAdaptProp);
dump("est.yp",file=paste("estimated.txt",sep=""),append=TRUE)
db$p.i <- est.yp
db$est.u <- est.u
if(use.v){
  est.v<-sum.v/(nbiterations-lastAdaptProp);
  dump("est.v",file=paste("estimated.txt",sep=""),append=TRUE)
  db$est.v <- est.v
}
if(use.insp){
  est.beta<-sum.beta/(nbiterations-lastAdaptProp);
  dump("est.beta",file=paste("estimated.txt",sep=""),append=TRUE)
  db$est.obs<-inspector%*%est.beta
}
if(use.cofactors){
  est.c.val<-sum.c.val/(nbiterations-lastAdaptProp);
  dump("est.c.val",file=paste("estimated.txt",sep=""),append=TRUE)
  db$est.c <- c.map%*%est.c.val
}
sampled<-sampled[1:(nbiterations+1),];

cat("\n")
attributes(db)$lastAdapt<-lastAdaptProp
attributes(db)$lastBurnIn<-lastBurnIn
attributes(db)$Kthin<-Kthin
attributes(db)$freqsave<-freqsave
attributes(db)$nbiterations<-nbiterations

save(list=ls(),file="EndSampleImage.img") # allow to examin the environment later

return(db)
}

extrapol.spatautocorel<-function(...){
  return(fit.spatautocorel(...,fit.spatstruct=FALSE))
}

cut.burnIn<-function(dbFit=NULL,samplesDB){
  if(!is.null(dbFit) && length(samplesDB)>0){
    # check that sampled not already chopped 
    sampleLength<-dim(samplesDB)[1]
    expectedLength<-attributes(dbFit)$nbiterations

    if(sampleLength == expectedLength){
      cat("Using only after burn in\n")
      range<-seq(attributes(dbFit)$lastBurnIn+1,sampleLength)
      samplesDB<-samplesDB[range,]
    }
  } 
  return(samplesDB)
}
get.sampled<-function(samples=NULL,file=monitorfile,dbFit=NULL){
  if(is.null(samples$sampled)){
    # for some reason this is much faster than read.table()
    cat("Importing spatial autocorrelation values\n")
    noms<-read.table(file,nrow=1,header=TRUE)
    sampled<-as.data.frame(matrix(scan(skip=1,file=file,sep="\t"),ncol=dim(noms)[2],byrow=TRUE))
    names(sampled)<-names(noms)
  }else{
    sampled<-samples$sampled
  }
  sampled<-cut.burnIn(dbFit,sampled)

  return(sampled)
}
get.cofactors<-function(samples=NULL,file=coffile,dbFit=NULL){
  if(is.null(samples$c.vals)){
    cat("Import cofactors\n")
    if(file.exists(file)){
      noms<-read.table(file,nrow=1,header=TRUE)
      c.vals<-as.data.frame(matrix(scan(skip=1,file=file,sep="\t"),ncol=dim(noms)[2],byrow=TRUE))
      names(c.vals)<-names(noms)
    }else{
      cat("Import cofactors impossible, no",file,"file \n")
      c.vals<-NULL
    }
  }else{
    c.vals<-samples$c.vals
  }

  c.vals<-cut.burnIn(dbFit,c.vals)

  return(c.vals)
}
get.betas<-function(samples=NULL,file=betafile,dbFit=NULL){
  if(is.null(samples$betas)){
    if(file.exists(file)){
      cat("Import observers \n")
      noms<-read.table(file,nrow=1,header=FALSE)
      betas<-as.data.frame(matrix(scan(file=file,sep="\t"),ncol=dim(noms)[2],byrow=TRUE))
      names(betas)<-paste("obs",1:length(names(betas)),sep="")
    }else{
      cat("Import observers failed: No",file,"file \n")
      betas<-NULL
    }

  }else{
    betas<-samples$betas
  }

  betas<-cut.burnIn(dbFit,betas)

  return(betas)
}
traces<-function(db,nl=3,nc=4){
  db<-as.data.frame(db)
  if(dim(db)[2]>100){
    type="p"
    pch<-"."
  }else{
    type="l"
  }
  for(num in 1:length(names(db))){
    if(num %% (nl*nc) ==1){ 
      dev.new()
      par(mfrow=c(nl,nc))
    }
    plot(db[[num]],main=names(db)[num],pch=pch,type=type)
  }
}

trace.mcmc<-function(samples=NULL,dbFit=NULL){
  # import autocorrelation parameters and likelihoods

  ## monitored variables/parameters
  sampled<-get.sampled(dbFit=dbFit,samples=samples)
  niterations<-dim(sampled)[1]

  if(!is.null(sampled$T)){
  par(mfrow=c(2,3))
  plot(1/sampled$T,main="Barrier effect (1/T)",pch=".")
  plot(sampled$f,main="Characteristic distance",pch=".")
  plot(c(1,niterations),c(0,1),type="n",main="Mean SB share",xlab="Index")
  lines(sampled$meanSBshare,pch=".")
  lines(sampled$meanSBshareNoT,main="Mean SB share",pch=".",col="blue",type="p")
  legend("bottomright",c("Actual","If no barrier effect"),pch=20,col=c("black","blue"),pt.cex=0.5)
  plot(sampled$meanBeta,main="Mean observer quality",pch=".")
  plot(sampled$LLH,main="LLH",pch=".")
  }else{
    par(mfrow=c(3,6))
    for( var in 1:length(names(sampled))){
      plot(sampled[,var],main=names(sampled)[var],pch=".")
    }
  }


  ## cofactors
  c.vals<-get.cofactors(dbFit=dbFit,samples=samples)
  if(!is.null(c.vals)){
    # traces it
    dev.new()
    par(mfrow=c(3,2))
    for(cofactor in names(c.vals)){
      plot(c.vals[,cofactor],main=cofactor,pch=".")
    }
  }

  ## observers
  betas<-get.betas(samples=samples,dbFit=dbFit,)
  if(!is.null(betas)){
    dev.new()
    plot(c(1,dim(betas)[1]),c(0,1),type="n",main="Individual quality of observers",xlab="Iterations",ylab="Detection rate of positive")
    for(insp in 1: dim(betas)[2]){
      lines(betas[,insp],col=sample(colors(),1),pch=".",type="p")
    }
  }

  return(invisible(list(sampled=sampled,c.vals=c.vals,betas=betas)))
}
get.estimate<-function(C,name="",visu=TRUE,leg=TRUE,true.val=NULL){
  C<-C[which(!is.infinite(C))]
  if(length(which(!is.na(C)))>1){
    estimate<-c(mean(C),quantile(C,probs=c(0.025,0.5,0.975)))
    names(estimate)[1]<-"Mean"
    names(estimate)[3]<-"Median"
  }else{
    estimate<-rep(NA,4)
  }
  if(length(levels(as.factor(C)))>1){ # avoid to estimate unvarying
    # fit the density of C, avoiding errors with overspread C 
    Nthin<-1
    while(class(densfit<-try(locfit(~C[seq(1,length(C),Nthin)])))=="try-error" && Nthin<length(C)/100){
      Nthin<-Nthin*10
    }
    attributes(densfit)$Nthin<-Nthin
    if(Nthin>1){
      name<-paste(name," thinned x",Nthin,sep="")
    }

    vals<-predict(densfit,estimate)
    if(visu){
      plot(densfit,xlab=name)
      lines(rep(estimate[1],2),c(0,vals[1]),col="black")
      for(q in 2:4){
	lines(rep(estimate[q],2),c(0,vals[q]),col="blue")
      }
      if(!is.null(true.val)){
	abline(v=true.val,col="green")
      }

      if(leg){ # legend
	# legend
	if(mean(densfit$box)>estimate[3]){
	  loc<-"topright"
	}else{
	  loc<-"topleft"
	}
	leg.text<-c(paste("CrI/med.",round(estimate[3],2)),paste("Mean",round(estimate[1],2)))
	leg.col<-c("blue","black")
	if(!is.null(true.val)){
	  leg.text<-c(leg.text,paste("True val.(",true.val,")",sep=""))
	  leg.col<-c(leg.col,"green")
	}
	legend(loc,leg.text,col=leg.col,lty=1)
      }
    }
  }else{
    densfit<-NULL
    vals<-NULL
  }
  attributes(estimate)$densfit<-densfit
  attributes(estimate)$vals<-vals

  return(estimate)
}
group.posteriors<-function(db,main=NULL,visu=TRUE,leg=NULL,true.vals=NULL,pal=strongColors){
  estimates<-list()
  if(!is.null(db)){
    db<-as.data.frame(db)
    if(!is.null(db)){
      maxdens<-0
      for(monitored in names(db)){
	estimates[[monitored]]<-get.estimate(db[,monitored],visu=FALSE)
	if(visu){
	  densfit<-attributes(estimates[[monitored]])$densfit
	  r<-densfit$box
	  vals<-seq(r[1],r[2],(r[2]-r[1])/100)
	  maxDensMon<-max(predict(densfit,vals))
	  maxdens<-max(maxdens,maxDensMon)
	}
      }
      if(visu){
	palette(pal)
	if(is.null(main)){
	  if(length(names(db))<length(palette())){
	    main<-paste(names(db))
	  }else{
	    main<-paste(names(db)[1],"...",tail(names(db),1))
	  }
	}
	plot(range(db),c(0,maxdens),type="n",main=main,xlab="Value",ylab="Density")
	for(monitored in names(db)){
	  col<-which(names(db)==monitored)
	  lines(attributes(estimates[[monitored]])$densfit,col=col)
	  if(!is.null(true.vals[monitored])){
	    abline(v=true.vals[monitored],col=col)
	  }
	}
	if(is.null(leg)){
	  leg<-(length(names(db))<length(palette()))
	}
	if(leg){
	  legend("topleft",names(db),col=(1:length(names(db))),lty=1)
	}
	palette("default")
      }
    }
  }
  return(invisible(estimates))
}

posteriors.mcmc<-function(samples=NULL,dbFit=NULL,visu=TRUE){
  ### plot posteriors and get estimates

  estimates<-list()
  ## monitored parameters and variables
  sampled<-get.sampled(dbFit=dbFit,samples=samples)

  if(!is.null(sampled$T)){
    if(visu){
      dev.new()
      par(mfrow=c(2,2))
    }
    estimates$T<-get.estimate(sampled$T,name="Across barrier factor (sigma)",leg=FALSE,visu=visu)
    estimates$f<-get.estimate(sampled$f,name="Charac. Dist. ",leg=FALSE,visu=visu)
    estimates$meanBeta<-get.estimate(sampled$meanBeta,name="Mean rate obs.",leg=FALSE,visu=visu)
    estimates$meanSBshare<-get.estimate(sampled$meanSBshare,name="Mean share SG",leg=TRUE,visu=visu)
  }else{
    if(visu){
      dev.new()
      par(mfrow=c(3,6))
    }
    for(var in 1: length(names(sampled))){
      nameVar<-names(sampled)[var]
      estimates[[nameVar]]<-get.estimate(sampled[,nameVar],name=nameVar,leg=FALSE,visu=visu)
    }
  }

  ## cofactors
  c.vals<-get.cofactors(dbFit=dbFit,samples=samples)
  estimates<-c(estimates,group.posteriors(c.vals,main="Cofactors' posteriors",visu=visu))
  ## inspectors
  betas<-get.betas(dbFit=dbFit,samples=samples)
  estimates<-c(estimates,group.posteriors(betas,main="Observers' posteriors",visu=visu))
  attributes(estimates)$sizeEstimate<-dim(sampled)[1]

  return(invisible(estimates))
}

summary.spatcontrol<-function(estimates=NULL,samples=NULL,dbFit=NULL){
  # compute it if needed
  if(is.null(estimates)){
    estimates<-posteriors.mcmc(visu=FALSE,samples=samples,dbFit=dbFit)
  } 
  # format it for easy printing
  toPrint<-t(as.data.frame(estimates))

  # keep only the necessary
  keepVar<-function(varVals){
    # do not keep NA filled
    NotGood1<-all(is.na(varVals))
    # do not keep invariant
    NotGood2<-all(varVals==varVals[1])
    return(!NotGood1 && !NotGood2)
  }
  toPrint<-toPrint[apply(toPrint,1,keepVar),]

  # reorder for better visibility
  toPrint<-toPrint[,c(1,3,2,4)]
  # names(toPrint)<-c("Mean","Median","2.5%","97.5%")
  print(toPrint,digits=2)

  return(invisible(toPrint))
}


