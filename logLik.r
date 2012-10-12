




## naive likelihoods and related....

ricker.fey <-  function(theta,y,e) {
## naive joint density of data y and random effects e for ricker model
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  n.t <- nrow(y)
  n.reps <- 1
  n <- y*0
  sig.e <- theta[2]
  oo <- .C("ricker",n=as.double(n),as.double(theta),as.double(e),as.integer(0),as.integer(n.t),
                    as.integer(n.reps),as.double(log(max(y[1],1))),PACKAGE="sl")
  mu <- exp(oo$n) ## E(y|e)
  sum (dpois(y,mu,log=TRUE)) + sum(dnorm(e,sd=sig.e,log=TRUE)) ## log f(y,e)

}


## default log synthetic likelihoods for models 


ricker.ll <- function(theta,y,e,u,burn.in,trans=NULL,stats=FALSE) {
## function to obtain log synthetic likelihood for the Ricker model
## as coded in `ricker'

## simulate from model with current parameter values

  Y <- ricker(theta,e,u,burn.in)

## Now assemble the relevant statistics
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  acf.Y <- sl.acf(Y,max.lag=5)
  acf.y <- sl.acf(y,max.lag=5)

  b0.Y <- nlar(Y^.3,lag=c(1,1),power=c(1,2))
  b0.y <- nlar(y^.3,lag=c(1,1),power=c(1,2))

  b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   

## combine the statistics...

  sy <- c(as.numeric(acf.y),
          as.numeric(b0.y),
          as.numeric(b1.y),
          mean(y),sum(y==0)
         )
  sY <- rbind(acf.Y,
              b0.Y, 
              b1.Y,
              colMeans(Y),
              colSums(Y==0)
             )


  if (!is.null(trans)) {
    sy <- trans.stat(sy,trans)
    sY <- trans.stat(sY,trans)
  }

## get the log synthetic likelihood

  sY <- sY[,is.finite(colSums(sY))]
  if (stats) {  ## return statistics
   attr(sY,"observed") <- sy
   return(sY) 
  }

  er <- robust.vcov(sY)

  rss <- sum((er$E%*%(sy-er$mY))^2)
  ll <- -rss/2 - er$half.ldet.V
  
} ## end of ricker.ll


marginal.norm <- function(sy,sY) {
## marginal normalization routine. Each statistic gets replaced by a 
## normal quantile.... doesn't work - only fixed transforms are allowed
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

trim.stat <- function(sY,p=.01) {
## trim the smallest and largest p of data from each row
  n <- ncol(sY)
  for (i in 1:nrow(sY)) {
    r <- rank(sY[i,])
    sY[i,r<=p*n|r>(n-p*n)] <- NA
  }
  sY
}

robust.vcov <- function(sY,alpha=2,beta=1.25) {
## Uses Campbell's robust approach as described on p 231 of Krzanowski 1988
## But adds pre-conditioning for stable computation....
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


robust.vcov.old <- function(sY,alpha=2,beta=1.25) {
## Uses Campbell's robust approach as described on p 231 of Krzanowski 1988
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

blowfly.ll <- function(theta,y,lu,lu1,burn.in,trans=NULL,stats=FALSE,step=2) {
## function to obtain log synthetic likelihood for the Nisbet and Gurney
## blowfly model as coded in `ng.bf' 
## also check the help for this funtion in the blowfly package
# theta: vector of current blowfly parameters (see ng.bf in sl.r)
#        c(P, N0, sig.p, tau, sig.d)
# y : vector of the observed abundances
# lu/lu1: noise matrices (see ng.bf in transform.r)
# burn.in: length of burn.in
# trans: if not null use trans to perform a "radical normalization of the summary statistic"
# stats: if false return only the stats and not the log(synthetic likelihood)
# step: weird, multiply the size of y?!?

## simulate from model with current parameter values

  n.y <- length(y)
  Y <- ng.bf(theta,lu,lu1,burn.in)[1:n.y*step,]

###################################
## Now assemble the relevant statistics
###################################
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  acf.Y <- sl.acf(Y,max.lag=11)
  acf.y <- sl.acf(y,max.lag=11)

  b0.Y <- nlar(Y,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))
  b0.y <- nlar(y,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))

  b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   

## combine the statistics for old and current proposed 

  sy <- c(as.numeric(acf.y),
          as.numeric(b0.y),
          as.numeric(b1.y),
          mean(y),
          mean(y)-median(y)#, ## possibly mean would be better here?
          #sum(abs(diff(sign(diff(y)))))/2 ## count turning points
         )
  sY <- rbind(acf.Y,
              b0.Y, 
              b1.Y,
              colMeans(Y),
              colMeans(Y)-apply(Y,2,median)#,
              #colSums(abs(diff(sign(diff(Y)))))/2
             )

###################################
## get the log synthetic likelihood
###################################
  # ## extreme transform to normality.... CB:what the heck?!?

  # if (!is.null(trans)) {
  #   sY <- trans.stat(sY,trans)
  #   sy <- trans.stat(sy,trans)
  # }



  # sY <- sY[,is.finite(colSums(sY))] # only keep the observation for which all stats are finite?

###   sY <- trim.stat(sY) ## trimming the marginal extremes to robustify

  # if (stats) { 
  #   attr(sY,"observed") <- sy 
  #   return(sY)   ## just return statistics
  # }

  # er <- robust.vcov(sY)

  # rss <- sum((er$E%*%(sy-er$mY))^2)
  # ll <- -rss/2 - er$half.ldet.V

  # attr(ll,"rss") <- rss
  # ll # the log synthetic likelihood
    ll<-synLik(sY,sy,trans)
  if(stats){
    sY<-attributes(ll)$sY
    sy<-attributes(ll)$sy
    attr(sY,"observed") <- sy 
    return(sY)   ## just return statistics
  }else{
    return(ll)
  }
  
} ## end of blowfly.ll

synLik<-function(sY,sy,trans=NULL){
## get the log synthetic likelihood
  # sY: matrix with stats for theta
  # sy: vector with stats in data

  ## extreme transform to normality.... CB:what the heck?!?
  if (!is.null(trans)) {
    sY <- trans.stat(sY,trans)
    sy <- trans.stat(sy,trans)
  }

  sY <- sY[,is.finite(colSums(sY))] # only keep the observation for which all stats are finite?

  ##  sY <- trim.stat(sY) ## trimming the marginal extremes to robustify

  er <- robust.vcov(sY)

  rss <- sum((er$E%*%(sy-er$mY))^2)
  ll <- -rss/2 - er$half.ldet.V

  attr(ll,"rss") <- rss
  attr(ll,"sy") <- sy
  attr(ll,"sY") <- sY

  return(ll)
}

dsbf.ll <- function(theta,y,burn.in,n.rep=500,trans=NULL,stats=FALSE,step=2) {
## function to obtain log synthetic likelihood for the Nisbet and Gurney
## blowfly model as coded in `ds.bf' --- i.e. for the model where *all*
## stochasticity is demographic.
## theta contains delta,P,N0 and tau, in that order.

## simulate from model with current parameter values

  n.y <- length(y)
  Y <- ds.bf(theta,burn.in=burn.in,n.t=n.y*step,n.rep=n.rep)[1:n.y*step,]

## Now assemble the relevant statistics
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  acf.Y <- sl.acf(Y,max.lag=11)
  acf.y <- sl.acf(y,max.lag=11)

  b0.Y <- nlar(Y,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))
  b0.y <- nlar(y,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))

  b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   

## combine the statistics...

  sy <- c(as.numeric(acf.y),
          as.numeric(b0.y),
          as.numeric(b1.y),
          mean(y),
          mean(y)-median(y) ## possibly mean would be better here?
          ##sum(abs(diff(sign(diff(y)))))/2 ## count turning points
         )
  sY <- rbind(acf.Y,
              b0.Y, 
              b1.Y,
              colMeans(Y),
              colMeans(Y)-apply(Y,2,median)
              ##colSums(abs(diff(sign(diff(Y)))))/2
             )

## extreme transform to normality....

  if (!is.null(trans)) {
    sY <- trans.stat(sY,trans)
    sy <- trans.stat(sy,trans)
  }


## get the log synthetic likelihood
  sY <- sY[,is.finite(colSums(sY))]

##  sY <- trim.stat(sY) ## trimming the marginal extremes to robustify

  if (stats) { 
    attr(sY,"observed") <- sy 
    return(sY)   ## just return statistics
  }

  er <- robust.vcov(sY)

  ## robustify the likelihood...
 
  rss <- sum((er$E%*%(sy-er$mY))^2)
  
  ll0 <- -rss/2 - er$half.ldet.V ## true l_s

  d0 <- qchisq(.99,nrow(sY))^.5

  rss <- not.sq(sqrt(rss),alpha=.1,d0=d0)
 
  ll <- -rss/2 - er$half.ldet.V ## robustified l_s
  attr(ll,"true") <- ll0 ## append the true l_s
  ll

} ## end of dsbf.ll

desbf.ll <- function(theta,y,burn.in,n.rep=500,trans=NULL,stats=FALSE,step=2) {
## function to obtain log synthetic likelihood for the Nisbet and Gurney
## blowfly model as coded in `des.bf' --- i.e. for the model where 
## stochasticity is demographic + environmental. Given demographic contribution
## it trivial this is not really useful --- might as well use bf.ll
## theta contains delta,P,N0,tau,sig2.p,sig2.d in that order.

## simulate from model with current parameter values

  n.y <- length(y)
  Y <- des.bf(theta,burn.in=burn.in,n.t=n.y*step,n.rep=n.rep)[1:n.y*step,]

## Now assemble the relevant statistics
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  acf.Y <- sl.acf(Y,max.lag=11)
  acf.y <- sl.acf(y,max.lag=11)

  b0.Y <- nlar(Y,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))
  b0.y <- nlar(y,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))

  b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   

## combine the statistics...

  sy <- c(as.numeric(acf.y),
          as.numeric(b0.y),
          as.numeric(b1.y),
          mean(y),
          mean(y)-median(y), ## possibly mean would be better here?
          sum(abs(diff(sign(diff(y)))))/2 ## count turning points
         )
  sY <- rbind(acf.Y,
              b0.Y, 
              b1.Y,
              colMeans(Y),
              colMeans(Y)-apply(Y,2,median),
              colSums(abs(diff(sign(diff(Y)))))/2
             )

## transform to normality....

  if (!is.null(trans)) {
    sY <- trans.stat(sY,trans)
    sy <- trans.stat(sy,trans)
  }


## get the log synthetic likelihood
  sY <- sY[,is.finite(colSums(sY))]

##  sY <- trim.stat(sY) ## trimming the marginal extremes to robustify

  if (stats) { 
    attr(sY,"observed") <- sy 
    return(sY)   ## just return statistics
  }

  er <- robust.vcov(sY)

  ## robustify the likelihood...
 
  rss <- sum((er$E%*%(sy-er$mY))^2)
  
  ll0 <- -rss/2 - er$half.ldet.V ## true l_s

  d0 <- qchisq(.99,nrow(sY))^.5

  rss <- not.sq(sqrt(rss),alpha=.1,d0=d0)
 
  ll <- -rss/2 - er$half.ldet.V ## robustified l_s
  attr(ll,"true") <- ll0 ## append the true l_s
  ll

} ## end of desbf.ll



bupalus.ll <- function(theta,y,e,u,burn.in,trans=NULL,stats=FALSE) {
## function to obtain log synthetic likelihood for the host parasite  model
## as coded in `bup.para'

## simulate from model with current parameter values

  Y <- bup.para(theta,e,u,burn.in)

## Now assemble the relevant statistics
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  acf.Y <- sl.acf(Y,max.lag=15)
  acf.y <- sl.acf(y,max.lag=15)

#  b0.Y <- nlar(Y,lag=c(1,1,2,2),power=c(1,2,1,2))
#  b0.y <- nlar(y,lag=c(1,1,2,2),power=c(1,2,1,2))

  b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   

## combine the statistics...

  sy <- c(as.numeric(acf.y),
        #  as.numeric(b0.y),
          as.numeric(b1.y),
          mean(y)
         )
  sY <- rbind(acf.Y,
         #     b0.Y, 
              b1.Y,
              colMeans(Y)
             )

## get the log synthetic likelihood
  sY <- sY[,is.finite(colSums(sY))]

  if (!is.null(trans)) {
    sy <- trans.stat(sy,trans)
    sY <- trans.stat(sY,trans)
  }

  if (stats) { 
    attr(sY,"observed") <- sy 
    return(sY)
  }

  er <- robust.vcov(sY)

  rss <- sum((er$E%*%(sy-er$mY))^2)
  ll <- -rss/2 - er$half.ldet.V

  ll

} ## end of bupalus.ll


chain2ll <- function(th,para=NULL,ll="ll",start=2000) {
## fits quadratic regression to chain output in th.
## rows of th contain chain output, and must include a 
## row of log likelihood obs. 
## para is an array of variable names.
## ll is the name of the log likelihood field
  
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
