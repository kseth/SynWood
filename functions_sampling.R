
#============================
# Perso MH
#===========================
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
omniSample<-function(Model,Data,oldTheta,nameParam,sdprop){
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
  }else{
	  stop("unknown sampling method for",nameParam)
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
  LLHold<-attributes(oldTheta)$outModel$LP

  # accept/reject
  # always 0 for symmetric distribution, only for record
  hasting_term<-dprop(old,prop,sdprop)-dprop(prop,old,sdprop);

  lnr <- LLHprop-LLHold+hasting_term;
  rand<-log(runif(1))
  # cat("otheta:",oldTheta,"ptheta",propTheta,"lnr:",lnr,"(",LLHprop,"-",LLHold,"+",hasting_term,"rand:",rand,"\n");
  
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


