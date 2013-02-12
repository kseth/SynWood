library(testthat)
# source("functions_sampling.R")
source("SynWood/functions_sampling.R",chdir=TRUE)
graphics.off()
test_that("getParamBeta OK",{
	  a<-9
	  b<-1
	  musigback<-getMeanSdBetaDis(a,b)
	  expect_equal(0.9,musigback$mu)

	  mu<-0.1
	  sig<-0.05
	  out<-getParamBeta(mu,sig)
	  musigback<-getMeanSdBetaDis(out$a,out$b)
	  expect_equal(mu,musigback$mu)
	  expect_equal(sig,musigback$sig)
})

test_that("omniSample OK for norm/lnorm",{
	  # simple model to test omniSample
	  Model<-function(theta,Data){
	    names(theta)<-Data$parm.names
	    LL<-sum(dnorm(Data$y,mean=theta["mean"],sd=theta["sd"],log=TRUE))
	    LP<-LL
	    yhat<-rnorm(length(Data$y),mean=theta["mean"],sd=exp(theta["sd"]))
	    return(list(LP=LP,
			Dev=-2*LL, # deviance, probably not to be changed
			Monitor=c(LP,theta), # to be monitored/ploted
			yhat=yhat, # data generated for that set of parameter
			# will be used for posterior check
			parm=theta # the parameters, possibly constrained by the model
			))
	  }

	  # simple Data to test normSample
	  nbit<-1000
	  upFreq<- 0 # display state every upFreq, 0 to never display
	  set.seed(777)
	  MyData<-list()
	  MyData$parm.names<-c("mean","sd")
	  MyData$sampling<-c("norm","lnorm")
	  MyData$mon.names<-c("LP",MyData$parm.names)
	  MyData$realMean<-0 # to be guessed
	  MyData$realSd<-1   # to be guessed
	  MyData$y<-rnorm(100,mean=MyData$realMean,sd=MyData$realSd)

	  # basic initialization
	  theta<-c(10,20)

	  #--- Necessary Machinery---#
	  # init of Data:
	  names(MyData$sampling)<-MyData$parm.names
	  nparams<-length(theta)

	  # init of theta attributes and saving scheme
	  outModel<-Model(theta,MyData)
	  Monitor<-mat.or.vec(nbit+1,length(MyData$mon.names))
	  Monitor[1,]<-outModel$Monitor
	  attributes(theta)$outModel<-outModel

	  accepts<-as.data.frame(matrix(rep(0,nparams*nbit),ncol=nparams))
	  names(accepts)<-MyData$parm.names

	  # simple MCMC chain
	  Rprof()
	  for(numit in 1:nbit){
		  if(numit%%upFreq==0 && upFreq!=0 ){
			  cat("it:",numit,"of",nbit,"current theta:",theta,"\n");
		  }
	    for(paramName in MyData$parm.names){
	      theta<-omniSample(Model,MyData,theta,paramName,0.4)
	      accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    }
	    Monitor[numit+1,]<-attributes(theta)$outModel$Monitor
	  }
	  Rprof(NULL)

	  # post treatment
	  Monitor<-as.data.frame(Monitor)
	  names(Monitor)<-MyData$mon.names
	  burn.in<-ceiling(nbit/10)
	  estMean<-mean(Monitor[-(1:burn.in),"mean"])
	  estSd<-mean(Monitor[-(1:burn.in),"sd"])
	  yMean<-mean(MyData$y)
	  ySd<-sd(MyData$y)

	  # cat("rateAccept:",apply(accepts,2,mean),"\n")
	  # cat("estimate(mean)",estMean,"estimate(sd)",estSd,"\n")

	  # dev.new()
	  # par(mfrow=c(1,3))
	  # plot(Monitor[,"LP"])
	  # plot(Monitor[,"mean"])
	  # plot(Monitor[,"sd"])

	  # dev.new()
	  # par(mfrow=c(1,2))
	  # hist(Monitor[,"mean"])
	  # abline(v=estMean,col="black")
	  # abline(v=MyData$realMean,col="blue")
	  # abline(v=yMean,col="red")
	  # hist(Monitor[,"sd"])
	  # abline(v=estSd,col="black")
	  # abline(v=MyData$realSd,col="blue")
	  # abline(v=ySd,col="red")

	  expect_true(abs(estMean-yMean)<0.1)
	  expect_true(abs(estSd-ySd)<0.1)
})

# test_that("omniSample OK for beta sampling",{
	  # guess a binomial rate
	  Model<-function(theta,Data){
	    names(theta)<-Data$parm.names
	    LL<-sum(dbinom(Data$draws,Data$nByDraw,theta[["rate"]],log=TRUE))
	    LP<-LL
	    yhat<-0
	    return(list(LP=LP,
			Dev=-2*LL, # deviance, probably not to be changed
			Monitor=c(LP,theta), # to be monitored/ploted
			yhat=yhat, # data generated for that set of parameter
			# will be used for posterior check
			parm=theta # the parameters, possibly constrained by the model
			))
	  }

	  # simple Data to test normSample
	  nbit<-1000
	  upFreq<- 0 # display state every upFreq, 0 to never display
	  set.seed(777)
	  MyData<-list()
	  MyData$parm.names<-c("rate")
	  MyData$sampling<-c("beta")
	  MyData$mon.names<-c("LP",MyData$parm.names)
	  MyData$realProba<-0 # to be guessed
	  MyData$nByDraw<-100
	  MyData$draws<-rbinom(100,MyData$nByDraw,MyData$realProba)

	  # basic initialization
	  theta<-c(0.5)

	  #--- Necessary Machinery---#
	  # init of Data:
	  names(MyData$sampling)<-MyData$parm.names
	  nparams<-length(theta)

	  # init of theta attributes and saving scheme
	  outModel<-Model(theta,MyData)
	  Monitor<-mat.or.vec(nbit+1,length(MyData$mon.names))
	  Monitor[1,]<-outModel$Monitor
	  attributes(theta)$outModel<-outModel

	  accepts<-as.data.frame(matrix(rep(0,nparams*nbit),ncol=nparams))
	  names(accepts)<-MyData$parm.names

	  # simple MCMC chain
	  Rprof()
	  for(numit in 1:nbit){
		  if(numit%%upFreq==0 && upFreq!=0 ){
			  cat("it:",numit,"of",nbit,"current theta:",theta,"\n");
		  }
	    for(paramName in MyData$parm.names){
	      theta<-omniSample(Model,MyData,theta,paramName,0.4)
	      accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    }
	    Monitor[numit+1,]<-attributes(theta)$outModel$Monitor
	  }
	  Rprof(NULL)

	  # post treatment
	  Monitor<-as.data.frame(Monitor)
	  names(Monitor)<-MyData$mon.names
	  burn.in<-ceiling(nbit/10)
	  estMean<-mean(Monitor[-(1:burn.in),"mean"])
	  estSd<-mean(Monitor[-(1:burn.in),"sd"])
	  yMean<-mean(MyData$y)
	  ySd<-sd(MyData$y)

	  # cat("rateAccept:",apply(accepts,2,mean),"\n")
	  # cat("estimate(mean)",estMean,"estimate(sd)",estSd,"\n")

	  # dev.new()
	  # par(mfrow=c(1,3))
	  # plot(Monitor[,"LP"])
	  # plot(Monitor[,"mean"])
	  # plot(Monitor[,"sd"])

	  # dev.new()
	  # par(mfrow=c(1,2))
	  # hist(Monitor[,"mean"])
	  # abline(v=estMean,col="black")
	  # abline(v=MyData$realMean,col="blue")
	  # abline(v=yMean,col="red")
	  # hist(Monitor[,"sd"])
	  # abline(v=estSd,col="black")
	  # abline(v=MyData$realSd,col="blue")
	  # abline(v=ySd,col="red")

	  expect_true(abs(estMean-yMean)<0.1)
	  expect_true(abs(estSd-ySd)<0.1)
# })


