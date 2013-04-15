# source("/home/cbarbu/Documents/outils/programmation/R/zoom.r")
source("zoom.r")
library(spam);
# sum on rows of a matrix
rowsum.spam <- function(A){
	return(A%*%rep(1,A@dimension[2]));
}
mean.spam<-function(A){
	sumMat<-sum(A@entries);
	nEntries<-A@dimension[1]*A@dimension[2];
	return(sumMat/nEntries)
}
quantile.spam<-function(A,...){
	A<-as.matrix(A)
	return(quantile(A,...))
}
# generic function that send back to the rowsum system
rsum<-function(x,...)
UseMethod("rowsum")

# get only the mean part of xprop for a given y
best.mv.canonical <- function(b,Q,Rstruct=NULL,...){
	N = dim(Q)[1]
	if (is(Rstruct, "spam.chol.NgPeyton")) 
		R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
	else R <- chol(Q, ...)
	if (is(R, "spam.chol.NgPeyton")) {
		mu <- drop(solve.spam(R, b))
	}
	else {
		mu <- backsolve(R, forwardsolve(t(R), b))
	}
	nu <- backsolve(R, array(rnorm(n * N), c(N, n)))
	return(drop(mu))
}

# check best.mv : check that each point is the mean of the neighbourghs
get_quality_of_spatial_mean<-function(u,Q){
	MeanMat<- -Q
	diag(MeanMat)<-0
	MeanVect<-MeanMat%*%u;
	sumSqDev<-t(MeanVect)%*%MeanVect;
	return(drop(sumSqDev));
}

# # basic test
# muu<-best.mv.canonical(y.r,Q)
# get_quality_of_spatial_mean(muu,Q)
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
	z_plot<-((z-min(z,asmin))/(max(z,asmax)-min(z,asmin)+0.00001))# ^(1/2) # +0.01 avoid NAN
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
rmvnorm.canonical.pseudo <- function (n, b, Q, Rstruct = NULL,given_rnorm=NULL, ...) 
{
	N = dim(Q)[1]
	if (is(Rstruct, "spam.chol.NgPeyton")) 
		R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
	else R <- chol(Q, ...)
	if (is(R, "spam.chol.NgPeyton")) {
		mu <- drop(solve.spam(R, b))
	}
	else {
		mu <- backsolve(R, forwardsolve(t(R), b))
	}
	if(is.null(given_rnorm)){
		given_rnorm<-rnorm(n * N)
	}
	nu <- backsolve(R, array(given_rnorm, c(N, n)))
	return(t(nu + mu))
}
rmvnorm.prec.pseudo <- function (n, mu = rep(0, nrow(Q)), Q, Rstruct = NULL,given_rnorm=NULL, ...) 
{
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
moran.spam<-function(dmt,values){
	# dmt a spam matrix with the weights between houses
	# values the values for each location

	# correction of dmt for correct calculus of I
	dim_mp<-length(data_non_NA[,1])
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
	names(data_points)<-c("easting","northing","block_num")
	
	return(data_points)
}


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

