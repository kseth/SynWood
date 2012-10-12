halfDistJ<-50
halfDistH<-10
rateMove<-0.0375
useDelta<-TRUE # if TRUE, weightSkipInMove<-0
delta<- 0.25
set.seed(777)

weightHopInMove<-1 # should not be changed, the ref for the two others
weightSkipInMove<-1 # ignored if useDelta == TRUE
weightJumpInMove<-0.07

# set corresponding rates
if (useDelta)
{
	weightSkipInMove<-0 
}

totalWeight<-(weightHopInMove+weightSkipInMove+weightJumpInMove)
rateHopInMove<-weightHopInMove/totalWeight
rateSkipInMove<-weightSkipInMove/totalWeight
rateJumpInMove<-weightJumpInMove/totalWeight # 1 - rateHopInMove - rateSkipInMove


################
# definition of the grid
################
library(splancs)
factSize<-1
L<-num_obs<-1600*factSize^2  #number of households
mx<-500*factSize   #length, in meters, of population boundary--squared is area
# starting point
infestHouse<-rep(0,L);
infestHouse[(L+sqrt(L))/2] <- 100

poly<- array(c(1,1,mx,mx,1,1,mx,mx,1,1), c(5,2))
gridpts(poly,L)->sp

# affect to blocks
nbblocks<- L/(16)
nbBlockHoriz<-(sqrt(L)/8)
nbBlockVert<-(sqrt(L)/2)

numBlock<-0
blockIndex<-rep(0,L)
for(jblock in 0:(nbBlockVert-1)){
	for(iblock in 0:(nbBlockHoriz-1)){
		numBlock<-numBlock+1
		xint<-((iblock*8):(iblock*8+7))+1
		yint<-((jblock*2):(jblock*2+1))+1
		# cat("x:",min(xint),max(xint),"y:",min(yint),max(yint),"BN:",numBlock,"\n")
		blockIndex[xint+(yint[1]-1)*sqrt(L)]<-numBlock
		blockIndex[xint+(yint[2]-1)*sqrt(L)]<-numBlock
	}
}
# plot(sp,col=blockIndex)
library(spam)

# define same block neighbors
spam.options(nearestdistnnz=c(13764100,400))
SB <- nearest.dist(x=cbind(blockIndex,rep(0,length(blockIndex))), method="euclidian", upper=NULL,delta=0.1)
SB@entries<-rep(1,length(SB@entries))
SB<-as.spam(SB);

source("spam_complement.r")

