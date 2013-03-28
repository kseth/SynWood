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

setClass("cb.diag.out",representation(colAnalysed = "vector", geweke= "vector",limitGeweke="numeric",raftery="matrix"),contains="data.frame")

gibbsit <- function(data, varnames = NULL, q = 1/40, r =
	  1/80, s = 19/20, epsilon = 1/1000, spacing = 1, resfile = "",
	  warning=TRUE,NminOnly=FALSE){

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

library(LaplacesDemon)

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
	burnIn<-resRafMax[2];
	Nmin<-resRafMax[4];
	if(nbItMin<=nbItBrut){
		cat("Raftery positive (",nbItBrut,">",nbItMin,")\n",file=logfile,append=TRUE);
		## if enough iterations according to the raftery, test if Geweke ok on non burnin 
		discard<-burnIn;
		result<-TRUE;
		nbRep<-1
		nMaxRep<-5
		limitGeweke<-1-(1-baseLimitGeweke)^(nMaxRep/length(sampThined[1,]))
		while( nbRep<=nMaxRep && min(resG)<limitGeweke && class(result)!="try-error"){
			cat("test Geweke with",discard,"discarded and limit=",limitGeweke,"\n",file=logfile,append=TRUE);
			resGbrut<-1+limitGeweke # somegthing bigger than limitGeweke
			result<-try(resGbrut<-Geweke.Diagnostic(sampThined[-(1:discard),]),silent=TRUE) # return p-values, 
			resG<-pnorm(-abs(resGbrut))
			# resG<-pnorm(-abs(Geweke.Diagnostic(sampThined[-(1:discard),]))) # return p-values, 
			discard<-ceiling((nbIt-burnIn)*nbRep*0.1);
			nbRep<-nbRep+1
		}
		if(min(resG)<limitGeweke){
			## if the convergence is not ok according to Geweke
			cat("Geweke not ok:",resG,"(threshold:",limitGeweke,")\n",file=logfile,append=TRUE);
			newNbIt<-ceiling(nbItBrut*1.5);
		}else{
			## if the convergence is ok test if we have enough samples
			cat("Geweke ok :",resG,"(threshold:",limitGeweke,"), burnIn=",discard*KthinInit,"\n",file=logfile,append=TRUE);
			resESS<-ESS(sampThined[-(1:discard),])
			nbItEff<-ceiling(max(resESS))
			if(nbItEff<Nmin){
				cat("Need some more iterations, only",nbItEff,"efficient iterations, we need:",Nmin,")\n",file=logfile,append=TRUE);
				## not yet enough iterations according to get the precision we need
				newNbIt<-as.integer(ceiling(Nmin*nbItBrut/nbItEff));
			}else{
				cat("All test passed, enough iterations (",nbItEff," it. eff. >",Nmin,")\n",file=logfile,append=TRUE);
				enoughIt<-TRUE
				newNbIt<-nbItBrut
			}
		}
	}else{
		cat("Raftery negative (",nbIt,"<",nbItMin,")\n",file=logfile,append=TRUE);
		newNbIt<-nbItMin;
	}
	output<-as.data.frame(t(c(enoughIt,newNbIt,nbItEff,KthinInit*max(discard,burnIn),Kthin,Kind,min(resG))))
	names(output)<-c("ok","newNbIt","nbItEff","burnIn","Kthin","Kind","min(resG)")

	return(new("cb.diag.out",output,colAnalysed=colAnalysed,geweke=resG,limitGeweke=limitGeweke,raftery=resRaf));
}


