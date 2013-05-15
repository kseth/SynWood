load("../whitehall_lat_long_pos.img")
whitehall <- out[, -which(names(out) %in% c("lat", "long"))]
whitehall <- cbind(whitehall, lat = out$long, lon = out$lat)
rm(out)

#convert lat, long to utms
conv <- LL.to.our.utms(whitehall[, c("lon", "lat")])
whitehall <- cbind(whitehall, conv)

maps <- unique(whitehall[, c("X", "Y", "name")])
startingInfested_names <- whitehall[which(whitehall$year == 2010), "name"]
startingInfested <- match(startingInfested_names, maps$name)
endInfestedHouses_names <- whitehall[which(whitehall$year == 2011), "name"]
endInfestedHouses <- match(endInfestedHouses_names, maps$name)
binomialEndInfested <- rep(0, length(maps$X))
binomialEndInfested[endInfestedHouses] <- 1

maps <- maps[, -3]

dev.new()
par(mfrow=c(2, 1))
st01 <- rep(0, length(maps$X))
st01[startingInfested] <- 1
plot_reel(maps$X, maps$Y, st01,  base=0)
plot_reel(maps$X, maps$Y, binomialEndInfested, base=0)
dev.off()
