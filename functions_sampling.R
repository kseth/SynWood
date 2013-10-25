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
  names(oldTheta)<-Data$parmNames
  old<-oldTheta[nameParam]
  # cat("parm:",nameParam,"old:",old,"oldTheta:",oldTheta,"\n")

  # sample proposal
  prop<-rnorm(1,mean=old,sd=sdprop);

  # include proposal in theta
  propTheta<-oldTheta
  attributes(propTheta)<-NULL # important to avoid growing thetas
  names(propTheta)<-Data$parmNames
  propTheta[nameParam]<-prop

  # get LLH for proposal
  outModel<-Model(propTheta,Data,...)
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
omniSample<-function(Model,Data,oldTheta,nameParam,sdprop,recompLLHold=TRUE){
  # identify param to sample
  names(oldTheta)<-Data$parmNames
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
  }else if(Data$sampling[nameParam]=="poisson"){
    rprop<-function(center,disp){
      return(rpois(1,center))
    }
    dprop<-function(val,center,disp){
      return(dpois(val,center,log=TRUE))
    }
  }else{
    stop("unknown sampling method for ",nameParam)
  }

  # sample proposal
  prop<-rprop(old,sdprop);

  # include proposal in theta
  propTheta<-oldTheta
  attributes(propTheta)<-NULL # important to avoid growing thetas
  names(propTheta)<-Data$parmNames
  propTheta[nameParam]<-prop

  # get LLH for proposal
  outModel<-Model(propTheta,Data)
  LLHprop<-outModel$LP
  LLHprop_prioronly <- outModel$Lprioronly
  LLHprop_statsonly <- outModel$LL
  
  #if recomputing or not recomputing LLHold
  if(!recompLLHold){
	#old LL stays constant	
	LLHold<-attributes(oldTheta)$outModel$LP
	LLHold_prioronly <- attributes(oldTheta)$outModel$Lprioronly
	LLHold_statsonly <- attributes(oldTheta)$outModel$LL
 
  }else{
	
	# set up dummy theta
	recompTheta<-oldTheta
	attributes(recompTheta)<-NULL
	names(recompTheta)<-Data$parmNames
	
	#recompute LLH for old
	recompModel<-Model(recompTheta, Data)
	LLHold<-recompModel$LP
	LLHold_prioronly <- recompModel$Lprioronly
	LLHold_statsonly <- recompModel$LL
		
	#redirect oldTheta
  	attributes(oldTheta)$outModel <- recompModel
  	attributes(oldTheta)$LLH <- LLHold
  }
  # avoid stumbling on a null prior
  if(is.null(LLHprop_prioronly)){
	  LLHprop_prioronly<-0
  }
  if(is.null(LLHold_prioronly)){
	  LLHold_prioronly<-0
  }

  # accept/reject
  # always 0 for symmetric distribution, only for record
  hasting_term<-dprop(old,prop,sdprop)-dprop(prop,old,sdprop)

  # compute the log accept/reject ratio
  if(is.finite(LLHprop) && is.finite(LLHold)){ # deal with NAN or Inf in synlik
	# to avoid problems with LLH_statsonly >> priors, nullifying the impact of the prior
	if(LLHprop_statsonly == LLHold_statsonly){ # if both stats only LL equal, consider only priors
		#browser()
		lnr <- LLHprop_prioronly-LLHold_prioronly+hasting_term
	}else{
  		lnr <- LLHprop-LLHold+hasting_term
	}
  }else{
	if(!is.na(LLHprop) && !is.na(LLHold)){ # if only infinites, no NA
		if(LLHprop > LLHold)
			lnr <- Inf
		else if(LLHprop < LLHold)
			lnr <- -Inf
		else 
			lnr <- LLHprop_prioronly-LLHold_prioronly+hasting_term
	}else{
		warning(paste("NAN from synLik: LLHprop",LLHprop,"LLHold",LLHold,"\n"))
		if(is.na(LLHprop) && is.na(LLHold))
			lnr <- LLHprop_prioronly-LLHold_prioronly+hasting_term
		else if(is.na(LLHold))
			lnr <- Inf
		else if(is.na(LLHprop))
			lnr <- -Inf
	}
  }
 
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




