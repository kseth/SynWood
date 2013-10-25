library(testthat)
# source("functions_sampling.R")
source("functions_sampling.R",chdir=TRUE)
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

start1 <- Sys.time()
test_that("omniSample OK for norm/lnorm",{
	  # simple model to test omniSample
	  Model<-function(theta,Data){
	    names(theta)<-Data$parmNames
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
	  MyData$parmNames<-c("mean","sd")
	  MyData$sampling<-c("norm","lnorm")
	  MyData$monNames<-c("LP",MyData$parmNames)
	  MyData$realMean<-0 # to be guessed
	  MyData$realSd<-1   # to be guessed
	  MyData$y<-rnorm(100,mean=MyData$realMean,sd=MyData$realSd)

	  # basic initialization
	  theta<-c(10,20)

	  #--- Necessary Machinery---#
	  # init of Data:
	  names(MyData$sampling)<-MyData$parmNames
	  nparams<-length(theta)

	  # init of theta attributes and saving scheme
	  outModel<-Model(theta,MyData)
	  Monitor<-mat.or.vec(nbit+1,length(MyData$monNames))
	  Monitor[1,]<-outModel$Monitor
	  attributes(theta)$outModel<-outModel

	  accepts<-as.data.frame(matrix(rep(0,nparams*nbit),ncol=nparams))
	  names(accepts)<-MyData$parmNames

	  # simple MCMC chain
	  # Rprof()
	  for(numit in 1:nbit){
		  if(numit%%upFreq==0 && upFreq!=0 ){
			  cat("it:",numit,"of",nbit,"current theta:",theta,"\n");
		  }
	    for(paramName in MyData$parmNames){
	      theta<-omniSample(Model,MyData,theta,paramName,0.4)
	      accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    }
	    Monitor[numit+1,]<-attributes(theta)$outModel$Monitor
	  }
	  # Rprof(NULL)

	  # post treatment
	  Monitor<-as.data.frame(Monitor)
	  names(Monitor)<-MyData$monNames
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

start2 <- Sys.time()
test_that("omniSample OK for [0, 1] sampling with boundednorm",{
	  # guess a binomial rate
	  Model<-function(theta,Data){
	    names(theta)<-Data$parmNames
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
	  MyData$parmNames<-c("rate")
	  MyData$sampling<-c("boundednorm")
	  MyData$monNames<-c("LP",MyData$parmNames)
	  MyData$realProba<-0 # to be guessed
	  MyData$nByDraw<-100
	  MyData$draws<-rbinom(100,MyData$nByDraw,MyData$realProba)

	  # basic initialization
	  theta<-c(0.5)

	  #--- Necessary Machinery---#
	  # init of Data:
	  names(MyData$sampling)<-MyData$parmNames
	  nparams<-length(theta)

	  # init of theta attributes and saving scheme
	  outModel<-Model(theta,MyData)
	  Monitor<-mat.or.vec(nbit+1,length(MyData$monNames))
	  Monitor[1,]<-outModel$Monitor
	  attributes(theta)$outModel<-outModel

	  accepts<-as.data.frame(matrix(rep(0,nparams*nbit),ncol=nparams))
	  names(accepts)<-MyData$parmNames

	  # simple MCMC chain
	  # Rprof()
	  for(numit in 1:nbit){
		  if(numit%%upFreq==0 && upFreq!=0 ){
			  cat("it:",numit,"of",nbit,"current theta:",theta,"\n");
		  }
	    for(paramName in MyData$parmNames){
	      theta<-omniSample(Model,MyData,theta,paramName,0.4)
	      accepts[numit,paramName]<-as.numeric(attributes(theta)$new)
	    }
	    Monitor[numit+1,]<-attributes(theta)$outModel$Monitor
	  }
	  # Rprof(NULL)

	  # post treatment
	  Monitor<-as.data.frame(Monitor)
	  names(Monitor)<-MyData$monNames
	  burn.in<-ceiling(nbit/10)
	  estProba<-mean(Monitor[-(1:burn.in),"rate"])

	  # cat("rateAccept:",apply(accepts,2,mean),"\n")
	  # cat("estimated rate",estProba,"\n")

	  # dev.new()
	  # par(mfrow=c(1,2))
	  # plot(Monitor[,"LP"], pch = ".")
	  # plot(Monitor[,"rate"], pch = ".")

	  expect_true(abs(estProba-MyData$realProba)<0.1)
})

end <- Sys.time()
# cat("norm MCMC ")
# print(start2-start1)
# cat("boundednorm MCMC ")
# print(end-start2)

#===============================
## intent to normalize any stat
#===============================
library(cobs)
z<-runif(1000)^2
## wood procedure
trans<-get.trans(rbind(z,z))
w<-trans.stat(rbind(z,z),trans)
hist(w[1,])
shapiro.test(w)
ks.test(w,"pnorm")
# # => shapiro no but ks ok and imply defining by hand, painful

### inverse of quantile based transform to normality
## prototype
par(mfrow=c(1,2))
z<-runif(1000)^2 # something really not normal
hist(z)
cz<-sort(z) # quantiles
invfn<-function(x) sqrt(x) # define an inverse of the quantile
	# it should always be possible as the quantile is injective
w<-qnorm(invfn(z)) # transform to normality
shapiro.test(w)
hist(w)
ks.test(w,"pnorm")
# => perfect

## real stuff
# initial state
par(mfrow=c(1, 4))
# z<-runif(1000,min=0,max=1)^2 # something really not normal
z<-rpois(1000, 50)
hist(z)

# fit a polynomial fn to the inverse of the quantiles
xs<-sort(z)
ys<-1:length(z)/(length(z)+1)

# # intent with polynomial fitting
# model<-lm(ys ~ poly(xs, 10, raw=TRUE))

# # intent with specific fitting
# model<-nls(ys ~ 1/(1+alpha*exp(-beta*xs)),start=list(alpha=1,beta=2))

# intent with spline fitting
model<-smooth.spline(xs,ys)
model2<-cobs(xs, ys, constraint = "increase", lambda=-1)
model3<-nls(ys ~ , start=list(beta=1,chi=1), lower=c(0, 0), algorithm="port")

#lines(xs,predict(model,xs)$y,col=4)
#lines(xs,predict(model2,xs)[, "fit"],col=5)
plot(xs,predict(model3,xs),type="l" col="grey")
points(xs,ys,pch=".")

invz<-predict(model,xs)$y
invz2<-predict(model2,xs)[, "fit"]
invz3<-predict(model3,xs)
ks.test(invz,"punif") # test uniformity
ks.test(invz2,"punif") # test uniformity
ks.test(invz3,"punif") # test uniformity
	
# uniform to normality
w<-qnorm(invz)
w2<-qnorm(invz2)
w3<-qnorm(invz3)
#hist(w)
#hist(w2)
hist(w3)
#hist(qnorm(sqrt(z)))
shapiro.test(w) 
shapiro.test(w2)
shapiro.test(w3)
# not quite as good as the real transform


## boxcox transforms
library(car)
