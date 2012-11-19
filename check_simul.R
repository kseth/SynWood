source("extrapol_field.R")
source("param.r")
Monitored<-read.table("thetasamples_all.txt",header=TRUE)
# obtain values from param.r
true.val<-c(rateMove, rateJumpInMove, delta, halfDistH, halfDistJ)
names(true.val)<-names(Monitored)[3:length(names(Monitored))]
traces(Monitored)

dev.new()
par(mfrow=c(2,3))
for(param in names(true.val)){
get.estimate(Monitored[[param]],true.val=true.val[param],name=param)
}
