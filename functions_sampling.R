library(msm)

#============================
# Functions for MH
#===========================
# linking sd/mean to parameters of beta distribution
getParamBeta<-function(mu,sig){
	K <- (1/mu -1)
	a<- mu *(mu^2*K/(sig^2) - 1)
	b<-K*a
	return(list(a=a,b=b))
}
getMeanSdBetaDis<-function(a,b){
	mu<- a/(a+b)
	sig<-sqrt(a*b/((a+b)^2 *(a+b+1)))
	return(list(mu=mu,sig=sig))
}
# function sampling a parameter fed to Model
normSample<-function(Model,Data,oldTheta,nameParam,sdprop){

  # identify param to sample
  names(oldTheta)<-Data$parm.names
  old<-oldTheta[nameParam]
  # cat("parm:",nameParam,"old:",old,"oldTheta:",oldTheta,"\n")

  # sample proposal
  prop<-rnorm(1,mean=old,sd=sdprop);

  # include proposal in theta
  propTheta<-oldTheta
  attributes(propTheta)<-NULL # important to avoid growing thetas
  names(propTheta)<-Data$parm.names
  propTheta[nameParam]<-prop

  # get LLH for proposal
  outModel<-Model(propTheta,Data)
  LLHprop<-outModel$LP
  LLHold<-attributes(oldTheta)$outModel$LP

  # accept/reject
  # always 0 for symmetric distribution, only for record
  hasting_term<-dnorm(old,prop,sdprop,log=TRUE)-dnorm(prop,old,sdprop,log=TRUE);

  lnr <- LLHprop-LLHold+hasting_term;
  # cat("otheta:",oldTheta,"ptheta",propTheta,"lnr:",lnr,"(",LLHprop,"-",LLHold,"+",hasting_term);
 
  if(lnr>=log(runif(1))) {
    newTheta <- propTheta;
    attributes(newTheta)$new<-TRUE
    attributes(newTheta)$outModel<-outModel
    # cat(nameParam," accept 1\n");
  }else{
    newTheta<-oldTheta
    attributes(newTheta)$new<-FALSE
    # cat(nameParam," accept 0\n");
  }

  return(newTheta)
}

# generic function for sampling, see test-functions_sampling.R for use
omniSample<-function(Model,Data,oldTheta,nameParam,sdprop, recompLLHold = TRUE){
  # identify param to sample
  names(oldTheta)<-Data$parm.names
  old<-oldTheta[nameParam]
  # cat("parm:",nameParam,"old:",old,"oldTheta:",oldTheta,"\n")

  # init rprop and dprop according to Data$sampling
  if(Data$sampling[nameParam]=="norm"){
	  rprop<-function(center,disp){
		  return(rnorm(1,mean=center,sd=disp))
	  }
	  dprop<-function(val,center,disp){
		  return(dnorm(val,mean=center,sd=disp,log=TRUE))
	  }
  }else if(Data$sampling[nameParam]=="lnorm"){
	  rprop<-function(center,disp){
		  return(rlnorm(1,meanlog=log(center),sdlog=disp))
	  }
	  dprop<-function(val,center,disp){
		  return(dlnorm(val,meanlog=log(center),sdlog=disp,log=TRUE))
	  }
  }else if(Data$sampling[nameParam]=="beta"){
	  rprop<-function(center,disp){
		  paramBeta<-getParamBeta(center,disp)
		  return(rbeta(1,paramBeta$a,paramBeta$b))
	  }
	  dprop<-function(val,center,disp){
		  paramBeta<-getParamBeta(center,disp)
		  return(dbeta(val,paramBeta$a,paramBeta$b,log=TRUE))
	  }
  }else if(Data$sampling[nameParam]=="boundednorm"){
	  rprop<-function(center,disp){
		  return(rtnorm(1,mean=center,sd=disp,lower=0,upper=1))
	  }
	  dprop<-function(val,center,disp){
		  return(dtnorm(val,mean=center,sd=disp,lower=0,upper=1,log=TRUE))
	  }
  }else{
	  stop("unknown sampling method for ",nameParam)
  }

  # sample proposal
  prop<-rprop(old,sdprop);

  # include proposal in theta
  propTheta<-oldTheta
  attributes(propTheta)<-NULL # important to avoid growing thetas
  names(propTheta)<-Data$parm.names
  propTheta[nameParam]<-prop

  # get LLH for proposal
  outModel<-Model(propTheta,Data)
  LLHprop<-outModel$LP
  
  #if recomputing or not recomputing LLHold
  if(!recompLLHold){
	
	#old LL stays constant	
	LLHold<-attributes(oldTheta)$outModel$LP
 
  }else{
	
	# set up dummy theta
	recompTheta<-oldTheta
	attributes(recompTheta)<-NULL
	names(recompTheta)<-Data$parm.names
	
	#recompute LLH for old
	recompModel<-Model(recompTheta, Data)
	LLHold<-recompModel$LP
		
	#redirect oldTheta
  	attributes(oldTheta)$outModel <- recompModel
  	attributes(oldTheta)$LLH <- LLHold
  }

  # accept/reject
  # always 0 for symmetric distribution, only for record
  hasting_term<-dprop(old,prop,sdprop)-dprop(prop,old,sdprop);

  lnr <- LLHprop-LLHold+hasting_term

  # if the sum is nan b/c of Infs, check them out
  if(is.nan(lnr))
	if(is.nan(LLHprop))
		lnr <- -Inf
	else if(LLHprop == -Inf)
		lnr <- -Inf
	else if(is.nan(LLHold))
		lnr <- hasting_term
	else if(LLHold == -Inf)
		lnr <- hasting_term

  rand<-log(runif(1))
  # cat("otheta: ",oldTheta[nameParam]," ptheta: ", propTheta[nameParam]," lnr: ",lnr," (",LLHprop,"-",LLHold,"+",hasting_term, ")", " rand: ", rand,"\n");
  
  if(lnr>=rand){
    newTheta <- propTheta;
    attributes(newTheta)$new<-TRUE
    attributes(newTheta)$LLH<-LLHold
    attributes(newTheta)$outModel<-outModel
    # cat(nameParam," accept 1\n");
  }else{
    newTheta<-oldTheta
    attributes(newTheta)$new<-FALSE
    # cat(nameParam," accept 0\n");
  }

  return(newTheta)
}

# adapt the standard deviation of the proposal
adaptSDProp <- function(sdprop, accepts, lowAcceptRate=0.15, highAcceptRate=0.4,tailLength=20){
  # according to Gelman1996
  # adapt the sampling deviation so that acceptance rate fall within:
  # 0.15 and 0.4 (careful to begin the sampling after that)

  acceptRate <- mean(tail(accepts,tailLength))
  cat("accept rate:",acceptRate);

  attributes(sdprop)$noupdate<-FALSE
  if(acceptRate < lowAcceptRate){
    newsdprop<-sdprop*0.9
    cat("update sdprop",sdprop,"to",newsdprop);
    return(newsdprop)
  }else if(acceptRate > highAcceptRate){
    newsdprop<-sdprop*1.1
    cat("update sdprop",sdprop,"to",newsdprop);
    return(newsdprop)
  }else{
    attributes(sdprop)$noupdate<-TRUE
    return(sdprop)
  }
}




