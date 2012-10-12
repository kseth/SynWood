library(testthat)
# test_that("normSample OK",{
	  # simple model to test normSample
	  Model<-function(theta,Data){
	    names(theta)<-Data$parm.names
	    LL<-sum(dnorm(Data$y,mean=theta["mean"],sd=exp(theta["sd"]),log=TRUE))
	    LP<-LL
	    yhat<-rnorm(length(Data$y),mean=theta["mean"],sd=exp(theta["logsd"]))
	    return(list(LP=LP,
			Dev=-2*LL, # deviance, probably not to be changed
			Monitor=c(LP,theta), # to be monitored/ploted
			yhat=yhat, # data generated for that set of parameter
			# will be used for posterior check
			parm=theta # the parameters, possibly constrained by the model
			))
	  }

	  # simple Data to test normSample
	  nbit<-10000
	  upFreq<-10
	  set.seed(777)
	  MyData<-list()
	  MyData$parm.names<-c("mean","logsd")
	  MyData$mon.names<-c("LP",MyData$parm.names)
	  MyData$realMean<-0
	  MyData$realSd<-1
	  MyData$y<-rnorm(100,mean=MyData$realMean,sd=MyData$realSd)

	  # basic initialization
	  theta<-c(1,log(2))
	  nparams<-length(theta)

	  outModel<-Model(theta,MyData)
	  Monitor<-mat.or.vec(nbit+1,length(MyData$mon.names))
	  Monitor[1,]<-outModel$Monitor
	  LLH<-outModel$LP

	  accepts<-as.data.frame(matrix(rep(0,nparams*nbit),ncol=nparams))
	  names(accepts)<-MyData$parm.names

	  # basic MCMC
	  Rprof()
	  for(numit in 1:nbit){
		  if(numit%%upFreq==0){
			  cat("it:",numit,"of",nbit,"current theta:",theta,"\n");
		  }
		  
	    for(paramName in MyData$parm.names){
	      theta<-normSample(Model,MyData,theta,paramName,0.4,LLH)
	      accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    }
	    Monitor[numit+1,]<-attributes(theta)$outModel$Monitor
	  }
	  Rprof(NULL)

	  # post treatment
	  Monitor<-as.data.frame(Monitor)
	  names(Monitor)<-MyData$mon.names
	  Monitor$sd<-exp(Monitor[,"logsd"])
	  estMean<-mean(Monitor[,"mean"])
	  estSd<-mean(Monitor[,"sd"])
	  yMean<-mean(MyData$y)
	  ySd<-sd(MyData$y)

	  cat("rateAccept:",apply(accepts,2,mean),"\n")
	  cat("estimate(mean)",estMean,"estimate(sd)",estSd,"\n")

	  par(mfrow=c(1,3))
	  plot(Monitor[,"LP"])
	  plot(Monitor[,"mean"])
	  plot(Monitor[,"sd"])

	  par(mfrow=c(1,2))
	  hist(Monitor[,"mean"])
	  abline(v=estMean,col="black")
	  abline(v=MyData$realMean,col="blue")
	  abline(v=yMean,col="red")
	  hist(Monitor[,"sd"])
	  abline(v=estSd,col="black")
	  abline(v=MyData$realSd,col="blue")
	  abline(v=ySd,col="red")

	  expect_true(abs(estMean-yMean)<0.1)
	  expect_true(abs(estSd-ySd)<0.1)

# })


