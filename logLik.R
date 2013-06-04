#=====================
# Transformation Functions
#=====================
## routine's for transformation of statistics to better meet normality assumptions,
## and for checking the MVN approximation

## S is an ns by n.reps matrix of statistics. This routine works through
## its rows finding piecewise linear transformations to normality, by 
## interactive use of `locator'.
get.trans <- function(S){
  op <- par(mfrow=c(1,1))
  if (!is.matrix(S)) S <- matrix(S,1,length(S))
  ns <- nrow(S)    ## the number of statistics
  n.rep <- ncol(S) ## the number of replicates
  z <- qnorm((1:n.rep-.5)/n.rep) ## standard normal quantiles
  trans <- list()
  for (i in 1:ns) { ## loop through stats...
    plot(sort(S[i,]),z)
    tr <- locator(,type="l",col=2)
    if (length(tr$x)<2) { 
      warning("no transform");
      trans[[i]] <- NULL  
    } else
    ## now extend end segments
    if (length(tr$x)==2) { ## single segment --- extend both ends
      xr <- tr$x[2]-tr$x[1]
      slope <- (tr$y[2]-tr$y[1])/xr
      tr$x[1] <- tr$x[1] - 1000*xr 
      tr$y[1] <- tr$y[1] - slope*1000*xr
      tr$x[2] <- tr$x[2] + 1000*xr
      tr$y[2] <- tr$y[2] + slope*1000*xr
      trans[[i]] <- tr      
    } else { ## extend end segments
      xr <- max(tr$x) - min(tr$x)
      slope <- (tr$y[2]-tr$y[1])/(tr$x[2]-tr$x[1])
      tr$x[1] <- tr$x[1] - 1000*xr 
      tr$y[1] <- tr$y[1] - slope*1000*xr
      nt <- length(tr$x)
      slope <- (tr$y[nt]-tr$y[nt-1])/(tr$x[nt]-tr$x[nt-1])
      tr$x[nt] <- tr$x[nt] + 1000*xr
      tr$y[nt] <- tr$y[nt] + slope*1000*xr
      trans[[i]] <- tr     
    }
  }
  trans[[ns+1]] <- NA
  par(op)
  trans ## the transformation object
}

## apply a piecewise linear `trans' object to the rows 
## of a statistics matrix, S
trans.stat <- function(S,trans) {
  if (!is.matrix(S)) S <- matrix(S,length(S),1)
  for (i in 1:nrow(S)) {
    if (!is.null(trans[[i]]))
    S[i,] <- approx(trans[[i]]$x,trans[[i]]$y,S[i,],rule=2)$y
  }
  if (ncol(S)==1) S <- as.numeric(S)
  S
}

#=================
# Visualization Function
#=================
## graphical check of multivariate normality, from
## Krzanowski (1988) 7.5
MVN.check <- function(S,s=NULL,cex.axis=1,cex.lab=1) {
  p <- nrow(S)
  n <- ncol(S)
  if (n<10*p) warning("You don't really have enough reps for this approach")
  
  ps <- s
  for (i in 1:nrow(S)) ps[i] <- sum(S[i,]<s[i])/ncol(S)

  um <- robust.vcov(S)
  ms <- as.numeric(um$mY) # means
  S <- S-ms # centering to 0
  ## Malahanobis for each column of S
  z <- colSums(S*(t(um$E)%*%um$E%*%S)) # standardize the variance
  
  q <- log(qchisq((1:n-.5)/n,df=p)) # expected quantiles of a multivariate standardized normal
  z <- log(sort(z)) # observed quantiles for the stats
  
  # plot the q-q multivariate plot
  dev.new()
  plot(q,z,type="l",col="grey",
       xlab="log theoretical quantiles",ylab="log observed quantiles",
       cex.axis=cex.axis,cex.lab=cex.lab)
  points(q,z,pch=".")
  abline(0,1,col=2)
  cat("\nproportion |log(z)-log(q)|>.25 = ",sum(abs(z-q)>.25)/n,"\n")
  
  if (!is.null(s)) { ## QQ plot for observed stats
    z <- um$E%*%(s-ms)
    q <- sum(z^2)
    abline(h=log(q),lty=2)
  }


  ## Marginal stats

  for (i in 1:nrow(S)) S[i,] <- S[i,]/um$sd[i]
  n <- ncol(S)
  z <- qnorm((1:n-.5)/n)
  rz <- range(z)
  rz <- c(rz[1]-2,rz[2]+2)
  
  # marginal qq-plot, one line per stat
  dev.new()
  plot(z,sort(S[1,]),type="l",col="grey",ylim=rz,
       xlab="N(0,1) quantiles",ylab="marginal quantiles",
       cex.axis=cex.axis,cex.lab=cex.lab)
  points(z,sort(S[1,]),pch=".")
  for (i in 2:nrow(S)) { 
    lines(z,sort(S[i,]),col="grey")
    points(z,sort(S[i,]),pch=".")
  }
  abline(0,1,col=2)  

  if (!is.null(s)) { ## QQ plot for observed stats
    z <- um$E%*%(s-ms)
    qqnorm(z,cex.axis=cex.axis,cex.lab=cex.lab);qqline(z,col=2)
  }

  ps
 
}

#=====================
# Syn Lik Computation Functions
#=====================
## marginal normalization routine. Each statistic gets replaced by a 
## normal quantile.... doesn't work - only fixed transforms are allowed
marginal.norm <- function(sy,sY) {
  ns <- nrow(sY)
  n <- ncol(sY)+1
  for (i in 1:ns) {
    x <- c(sy[i],sY[i,]) ## the data
    x <- qnorm((rank(x)-.5)/n)
    sy[i] <- x[1]
    sY[i,] <- x[-1]
  }
  list(sy=sy,sY=sY)
}

## trim the smallest and largest p of data from each row
trim.stat <- function(sY,p=.01) {
  n <- ncol(sY)
  for (i in 1:nrow(sY)) {
    r <- rank(sY[i,])
    sY[i,r<=p*n|r>(n-p*n)] <- NA
  }
  sY
}

## Uses Campbell's robust approach as described on p 231 of Krzanowski 1988
## But adds pre-conditioning for stable computation....
robust.vcov <- function(sY,alpha=2,beta=1.25) {
 
  mY <- rowMeans(sY)
  sY1 <- sY - mY 
  ## use pre-conditioning to stabilize computation
  D <- rowMeans(sY1*sY1)^.5 
  Di <- 1/D  ## diagonal pre-conditioner
  
  sY1 <- Di*sY1 ## pre-conditioned for better scaling
  R <- qr.R(qr(t(sY1)))/sqrt(ncol(sY1)-1) ## Va = DR'RD - initial estimate
  zz <- forwardsolve(t(R),sY1)
  d <- sqrt(colSums((zz)^2)) ## Mahalonobis distance for each column

  ## create Campbell weight vector...
  d0 <- sqrt(nrow(sY)) + alpha/sqrt(2)
  w <- d*0 + 1
  ind <- d>d0
  w[ind] <- exp(-.5*(d[ind]-d0)^2/beta)*d0/d[ind] 
  mY <- colSums(w*t(sY))/sum(w)
  sY1 <- sY - mY

  ## preconditioning...
  D <- rowMeans(sY1*sY1)^.5
  Di <- 1/D  ## diagonal pre-conditioner
  sY1 <- Di*sY1 ## pre-conditioned for better scaling
  
  R <- qr.R(qr(w*t(sY1)))/sqrt(sum(w*w)-1) ## Va = DR'RD
  sd <- rowSums((D*t(R))^2)^.5
  E <- t(Di*backsolve(R,diag(nrow(R))))                   ## V^{-1} = E'E 
  half.ldet.V <- sum(log(abs(diag(R)))) + sum(log(D))

  list(E=E,half.ldet.V=half.ldet.V,mY=mY,sd=sd)

}

## Uses Campbell's robust approach as described on p 231 of Krzanowski 1988
## But adds pre-conditioning for stable computation....
## modified by KS
#### avoid the QR decomposition used for preconditioning in robust.vcov which seems to be acting up
#### instead, use the eigendecomposition provided in robust.vcov.old
robust.vcov.modified <- function(sY,alpha=2,beta=1.25) {
 
  mY <- rowMeans(sY)
  sY1 <- sY - mY 
  ## use pre-conditioning to stabilize computation
  D <- rowMeans(sY1*sY1)^.5 
  Di <- 1/D  ## diagonal pre-conditioner
  
  sY1 <- Di*sY1 ## pre-conditioned for better scaling
  Va <- sY1%*%t(sY1)/(ncol(sY1)-1) #var-cov matrix
  ev <- eigen(Va,symmetric=TRUE) #eigendecomp
  zz <- t(ev$vectors)%*%sY1
  ival <- ev$values
  ind <- ival > ival[1]*.Machine$double.eps^.9
  ival[ind] <- 1/sqrt(ival[ind])
  ival[!ind] <- 0
  d <- sqrt(colSums((ival*zz)^2)) ## mahalanobis distance (arrived at without using QR decomposition) 

  ## create Campbell weight vector...
  d0 <- sqrt(nrow(sY)) + alpha/sqrt(2)
  w <- d*0 + 1
  ind <- d>d0
  w[ind] <- exp(-.5*(d[ind]-d0)^2/beta)*d0/d[ind] 
  mY <- colSums(w*t(sY))/sum(w)
  sY1 <- sY - mY

  ## preconditioning...
  D <- rowMeans(sY1*sY1)^.5
  Di <- 1/D  ## diagonal pre-conditioner
  sY1 <- Di*sY1 ## pre-conditioned for better scaling

  R <- qr.R(qr(w*t(sY1)))/sqrt(sum(w*w)-1) ## Va = DR'RD
  sd <- rowSums((D*t(R))^2)^.5
  E <- t(Di*backsolve(R,diag(nrow(R))))                   ## V^{-1} = E'E 
  half.ldet.V <- sum(log(abs(diag(R)))) + sum(log(D))

  list(E=E,half.ldet.V=half.ldet.V,mY=mY,sd=sd)
}

## Uses Campbell's robust approach as described on p 231 of Krzanowski 1988
robust.vcov.old <- function(sY,alpha=2,beta=1.25) {

  mY <- rowMeans(sY)
  sY1 <- sY - mY 
 
  Va <- sY1%*%t(sY1)/(ncol(sY1)-1)
  if (FALSE) {
    R <- chol(Va)
    M <- forwardsolve(t(R),sY1)
    d <- colSums(M*M)^.5 ## Mahalanobis distance for cols of sY
  } else {
    ev <- eigen(Va,symmetric=TRUE)
    zz <- t(ev$vectors)%*%sY1
    ival <- ev$values
    ind <- ival > ival[1]*.Machine$double.eps^.9
    ival[ind] <- 1/sqrt(ival[ind])
    ival[!ind] <- 0
    d <- sqrt(colSums((ival*zz)^2)) 
  }
  d0 <- sqrt(nrow(sY)) + alpha/sqrt(2)
  w <- d*0 + 1
  ind <- d>d0
  w[ind] <- exp(-.5*(d[ind]-d0)^2/beta)*d0/d[ind] 
  mY <- colSums(w*t(sY))/sum(w)
  sY1 <- sY - mY
  w2 <- w*w
  Va <- sY1%*%(w2*t(sY1))/(sum(w2)-1)
  list(mY = mY, Va=Va)
}

## get the log synthetic likelihood
  # sY: matrix with stats for theta (not necessary if er given)
  # sy: vector with stats in data
  # trans: a result of call get.trans
  #	   contains piecewise transform to normality (is interactive)
  # er: if given, use this vcov object instead of vcov of sY
  #	if given, will override sY
  #	if given, must be same format as object from robust.vcov (or alternately robust.vcov.modified)
synLik<-function(sY=NULL, sy, trans=NULL, er=NULL){

  if(is.null(sY) && is.null(er))
  	stop("either sY (stats) or er (vcov) must be passed")	  

  ## extreme transform to normality
  if (!is.null(trans)){
    if(!is.null(sY))
	    sY <- trans.stat(sY,trans)
    
    sy <- trans.stat(sy,trans)
  }

  if(is.null(er)){ ## if er is not passed, calc vcov matrix

        ## sY <- trim.stat(sY) ## trimming the marginal extremes to robustify
  	## commented out by Wood, KS

	# only keep the observation for which all stats are finite
  	sY <- sY[,is.finite(colSums(sY))] 	

	er <- try(robust.vcov(sY), silent=TRUE)

	if(class(er)=="try-error" || is.na(er)){ ## cannot find vcov matrix
		ll<-NA
		attr(ll,"rss") <- NA
		attr(ll,"er") <- NA
		attr(ll,"sy") <- sy
		attr(ll,"sY") <- sY
		return(ll)
  	}
  }

  rss <- sum((er$E%*%(sy-er$mY))^2)
  ll <- -rss/2 - er$half.ldet.V

  attr(ll,"rss") <- rss
  attr(ll,"sy") <- sy
  attr(ll,"sY") <- sY 
  attr(ll,"er") <- er

  return(ll)
}

## fits quadratic regression to chain output in th.
## rows of th contain chain output, and must include a 
## row of log likelihood obs. 
## para is an array of variable names.
## ll is the name of the log likelihood field 
chain2ll <- function(th,para=NULL,ll="ll",start=2000) {

  ## get default predictor names
  if (is.null(para)) para <- rownames(th)
  para <- para[para!=ll]
  
  ## separate the ll and parameter information
  llr <- th[rownames(th)==ll,] 
  th <- th[rownames(th)!=ll,]
  
  ## discard burn-in and centre
  n.mc <- ncol(th)
  th <- th[,start:n.mc] ## discard burn in
  llr <- llr[start:n.mc] 
  thm <- rowMeans(th)
  th <- th - thm ## centre variables
  
  ## start constructing the regression formula
  form <- paste(ll,"~",paste(para,collapse=" + "))
  m <- length(para) 
  for (i in 1:m) form <- paste(form," + I(",para[i],"^2)",sep="")
  for (i in 1:m) if (i<m) for (j in (i+1):m) form <- 
     paste(form," + I(",para[i]," * ",para[j],")",sep="")

  df <- rbind(llr,th)
  rownames(df)[1] <- ll
  df <- as.data.frame(t(df))

  model <- lm(form,data=df)
  b <- coef(model)

  ## extract the Hessian of the loglik

  H <- matrix(0,m,m) 
  k <- 2 * m + 1
  for (i in 1:m) if (i<m) for (j in (i+1):m) 
  { k <- k + 1
    H[i,j] <- H[j,i] <- b[k]
  }
  k <- m + 1
  for (i in 1:m) {k<- k+1;H[i,i] <- 2*b[k]}  
  rownames(H) <- colnames(H) <- para
  
  eh <- eigen(-H,symmetric=TRUE)
  ev <- eh$values;
  ind <- abs(ev)>max(abs(ev))*.Machine$double.eps^.9
  ev[ind] <- 1/ev[ind]
  ev[!ind] <- 0
  Hi <- eh$vectors%*%(ev*t(eh$vectors))

  th.mle <- as.numeric(Hi%*%b[2:(m+1)]) + thm
  th.se <- diag(Hi)^.5

  ml <- predict(model,newdata=as.list(th.mle-thm),se=TRUE)
  
  list(logLik=ml$fit,logLik.se=ml$se.fit,mle=th.mle,se=th.se,Hi=Hi,H=H)
}
