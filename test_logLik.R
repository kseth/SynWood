source("logLik.R")

set.seed(13)
n <- 5    # dimension
N <- 1000  # sample size
Sigmainv <- .25^abs(outer(1:n,1:n,"-"))
Sigmainv <- as.spam( Sigmainv, eps=1e-4)
Mean<- rnorm(n,sd=3)
mvsample<-rmvnorm.prec(N,mu=Mean,Q=Sigmainv)

sl1<-synLik(sY=t(mvsample),sy=mvsample[1,])
er<-attr(sl1,"er")
attributes(sl1)<-NULL
sl2<-synLik(sY=t(mvsample),sy=mvsample[2,])
attributes(sl2)<-NULL

sls<-synLik(sY=t(mvsample),sy=t(mvsample[1:2,]))
attributes(sls)<-NULL

expect_equal(c(sl1,sl2),sls)
