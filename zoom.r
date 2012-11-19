# allow to to zoom on last existing plot using zm()
# optionally, graphs can be plotted with zplot and lines with zlines
# allowing to go back and forth with bf()

# WARNING: when you zoom on a multiplot window, you can only zoom on the last plot and the other plots will be rescaled
# accordingly, can be very usefull or anoying but I don't plan to change this soon

library(Hmisc); # for errbar

full.screen<-function(...){
	dev.new(width=13.5,height=9.5,...)
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

zoom_loaded<-1

multipanc<-function(ancien,fact,lim){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
	meanRange<-mean(ancien)
	newRange<-ancien+(1/fact-1)*(ancien-meanRange)
	return(newRange);
}
keepanc<-function(ancien,fact,lim){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
	return(ancien);
}
usenew<-function(ancien,fact,lim){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
	return(lim);
}
zoomplot.zoom <- function (xlim=NULL, ylim = NULL,fact=NULL,rp=NULL,...) 
{
	# rp is a recorded plot
	# fact is a factor of manification/outzoom
	# fact has priority on xlim/ylim
	if(! is.null(fact)&& is.numeric(fact) && length(fact)==1){
		xlimfn <-multipanc;
		ylimfn <-multipanc;
	}else{
		if(!is.null(xlim)&& is.numeric(xlim) && length(xlim)==2){
			xlimfn <- usenew;
		}else{
			xlimfn <-keepanc;
		}
		if(!is.null(ylim) && is.numeric(ylim) && length(ylim)==2){
			ylimfn <-usenew;
		}else{
			ylimfn <-keepanc;
		}
	}

	if(is.null(rp)){
		tmp <- recordPlot()[[1]]
	}else{
		tmp<-rp[[1]]
	}

	for (i in seq(along = tmp)) {
		fn <- tmp[[i]][[1]]
		alst <- as.list(tmp[[i]][[2]])
		tmp2 <- all.equal(.Primitive("locator"), fn)
		if (is.logical(tmp2) && tmp2) {
			next
		}
		tmp2 <- all.equal(.Primitive("plot.window"), fn)
		if (is.logical(tmp2) && tmp2) {
			# print(alst)
			# cat("alst orig:",alst[[1]],alst[[2]],"\n")
			alst[[1]] <- xlimfn(alst[[1]],fact,xlim)
			alst[[2]] <- ylimfn(alst[[2]],fact,ylim)
		}
		do.call(fn, alst)
	}
}
# get the limits of the plot in argument
getalst<-function(tmp=recordPlot()[[1]]){
	for (i in seq(along = tmp)) {
		fn <- tmp[[i]][[1]]
		alst <- as.list(tmp[[i]][[2]])
		tmp2 <- all.equal(.Primitive("locator"), fn)
		if (is.logical(tmp2) && tmp2) {
			next
		}
		tmp2 <- all.equal(.Primitive("plot.window"), fn)
		if (is.logical(tmp2) && tmp2) {
			alstwin<-alst
		}
		# do.call(fn, alst)
	}
	return(alstwin);
}

inout.zoom<-function(...){
	cat("Enter magnification factor: \n")
	f<-scan(n=1)
	zoomplot.zoom(fact=f,...);
	return()
}
in.zoom<-function(...){
	# Ideally later:
	# center<-locator(1)
	# cat("Click on center: \n")
	zoomplot.zoom(fact=2,...);
	return()
}
out.zoom<-function(...){
	# Ideally later:
	# center<-locator(1)
	# cat("Click on center: \n")
	zoomplot.zoom(fact=0.5,...);
	return()
}
orig.zoom<-function(orig){
	replayPlot(orig)
	return()
}
sq.zoom<-function(...){
	# use locator to zoom with the mouse (two left clicks)
	# specially, ... can be used to pass a recorded plot rp
	cat("Click two times on the plot to define a square to zoom in\n");
	square<-locator(2)
	xmin<-min(square$x)
	xmax<-max(square$x)
	ymin<-min(square$y)
	ymax<-max(square$y)
	zoomplot.zoom(xlim=c(xmin,xmax),ylim=c(ymin,ymax),...)
}
session.zoom<-function(...){
	orig <- recordPlot()
	go_on<-TRUE
	while(go_on){
		cat("Interactive plot, you can\n")
		cat("    continue: Enter\n")
		cat("     zoom in: 1\n")
		cat("    zoom out: 2\n")
		cat(" zoom square: 3\n")
		cat(" set magnif.: 4\n")
		cat("back to init: 9\n")
		sel<-scan(n=1)
		if(length(sel)==0){
			go_on<-FALSE;
		}else{
			if(exists("exec.zoom")){
				rm(exec.zoom)
			}
			exec.zoom<-switch(sel,
				in.zoom,out.zoom,
				sq.zoom,inout.zoom,
				NULL,NULL,NULL,NULL,
				orig.zoom)
			if(!is.null(exec.zoom)){
				exec.zoom(orig=orig,...);
			}else{
				cat("Say it again?\n")
			}
		}
	}
	return(recordPlot());
}
zm<-session.zoom


# # example:
# plot(runif(1000))
# abline(0,1)
# original<-session.zoom()
# 
# # and to go back to the original
# replayPlot(original);

# in addition, launch any plot with a zoom session
ord.zoom<-list()
num.zoom<-0
max.zoom<-0
plot.zoom<-function(zs=1,...){
	out<-try(plot(...))
	if(length(out)==0){ # success of the try
		ord.zoom<<-list()
		num.zoom<<-0
		max.zoom<<-0
		if(zs==1){
			original<-session.zoom()
		}else{
			original<-recordPlot()
		}
		ord.zoom[[1]] <<- original
		num.zoom <<-1# do it after not to change it if failure
		max.zoom <<-1
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
	return(original)
}
zplot<-plot.zoom

# # initial prototype
# lines.zoom<-function(...){
# 	out<-try(lines(...))
# 	if(length(out)==0){ # success of the try
# 		recplot<-session.zoom()
# 		ord.zoom[[num.zoom+1]] <<- recplot
# 		num.zoom<<-num.zoom+1 
# 		max.zoom <<-num.zoom
# 	}
# 	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
# 	return(recplot)
# }
fn.zoom<-function(fn,zs=1,...){
	out<-try(fn(...))
	if(length(out)==0){ # success of the try
		if(zs==1){
			recplot<-session.zoom()
		}else{
			recplot<-recordPlot();
		}
		ord.zoom[[num.zoom+1]] <<- recplot
		num.zoom<<-num.zoom+1 
		max.zoom <<-num.zoom
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
	return(recplot)
}
zabline<-function(...){
	fn.zoom(abline,...)
}
zlines<-function(...){
	fn.zoom(lines,...)
}
zplot<-function(...){
	fn.zoom(plot,...)
}
zerrbar<-function(...){
	fn.zoom(errbar,...)
}
# function that plot the recorded plot num_to_plot 
# with the same zoom than the plot num_w_zoom
# or not according to keepzoom (1 yes, other no)
plotsamezoom.zoom<-function(num_to_plot,num_w_zoom,keepzoom){
	rp<-ord.zoom[[num_to_plot]] # get recorded plot
	if(keepzoom==T){ # put zoom of now to recorded plot
		if(! is.null(num_w_zoom)){
			wz<-ord.zoom[[num_w_zoom]]  # get with zoom plot
		}else{
			wz<-recordPlot();
		}
		lim<-getalst(tmp=wz[[1]]);
		zoomplot.zoom(xlim=lim[[1]],ylim=lim[[2]],rp=rp);
	}else{
		replayPlot(rp)
	}
}
back.zoom<-function(keepzoom=F,...){
	if(num.zoom-1<1){
		cat("Sorry, already at first item\n");
	}else{
		plotsamezoom.zoom(num.zoom-1,NULL,keepzoom)
		num.zoom <<- num.zoom-1;
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
}
forward.zoom<-function(keepzoom=F,...){
	if(num.zoom+1>max.zoom){
		cat("Sorry, already at final item\n");
	}else{
		plotsamezoom.zoom(num.zoom+1,NULL,keepzoom)
		num.zoom <<- num.zoom+1;
		
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
}
bkeepzoom.zoom<-function(...){
	back.zoom(keepzoom=T)
}
fkeepzoom.zoom<-function(...){
	forward.zoom(keepzoom=T)
}

bf.zoom<-function(){
	orig <- recordPlot()
	go_on<-T
	cat("You can go back an forth in saved graphs:\n");
	cat("1: backward\n");
	cat("2: forward\n");
	cat("4: backward same zoom\n");
	cat("5: forward same zoom\n");
	cat("9: original\n");
	while(go_on){
		sel<-scan(n=1)
		if(length(sel)==0){
			go_on<-FALSE;
		}else{
			if(exists("exec.zoom")){
				rm(exec.zoom)
			}
			exec.zoom<-switch(sel,
				back.zoom,
				forward.zoom,
				NULL,
				bkeepzoom.zoom,
				fkeepzoom.zoom,
				NULL,
				NULL,
				NULL,
				orig.zoom
				)
			if(! is.null(exec.zoom)){
				exec.zoom(orig=orig);
			}else{
				cat("Say it again?\n")
			}
		}
	}
	return(recordPlot())
}
bf<-bf.zoom


# # example:
# zplot(runif(1000))
# zlines(c(100,300),c(0.2,0.4))
# zabline(v=600)
#
# # then to go back 
# back.zoom()
# # or to go forward 
# forward.zoom()

