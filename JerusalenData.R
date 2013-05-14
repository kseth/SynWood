## csv of Jerusalen
## "Encuesta_melgar.csv" contains data from the Encuesta for all of Mariana Melgar ~ 2008
## "jerusalen_encuesta_2011_f2_collapsed_KHL.csv" contains the 2011 data for just jerusalen
encuesta_melgar <- read.csv(file = "../Encuesta_Melgar.csv")
jerusalen_encuesta <- read.csv(file = "../jerusalen_encuesta_2011_full_f2_collapsed_KHL.csv")

keep <- which(encuesta_melgar$unicode %in% jerusalen_encuesta$Unicode)
match <- match(jerusalen_encuesta$Unicode, encuesta_melgar[keep, "unicode"])
jer_dat <- cbind(jerusalen_encuesta[, -which(names(jerusalen_encuesta) %in% c("EASTING", "NORTHING"))], OLDSTATUS = encuesta_melgar[keep, "status"][match])
jer_dat <- cbind(jer_dat, X = jerusalen_encuesta[, "EASTING"], Y = jerusalen_encuesta[, "NORTHING"])
maps <- jer_dat[, c("X", "Y")]
## fields now contained in jer_dat:
## "Unicode","POINT_X","POINT_Y", "STATUS","D.x","L.x","V.x","BLOCK_NUM","TOTAL_C","TOTAL_P", "OLDSTATUS", "X", "Y"

startInfestH <- which(maps$OLDSTATUS == 1)
endInfestH <- which(maps$STATUS == 1)

## create a dummy timeH
timeH <- rep(-2, length(startInfestH))

### the vec of stats for the second timepoint data
stats <- noKernelMultiGilStat(stratHopSkipJump = stratHopSkipJump, blockIndex = blockIndex, infestH = endInfestH, timeH=timeH, endTime = nbit, rateMove = rateMove, weightSkipInMove = weightSkipInMove, weightJumpInMove = weightJumpInMove, Nrep = 1, coords = maps, breaksGenVar = genIntervals, simul=FALSE, getStats = TRUE, seed = seedSimul, dist_out = bin_dist_out, typeStat = useStats, map.partitions = map.partitions, conc.circs = circles, rateIntro = rateIntro)

if(!is.vector(stats$statsTable)){
	statsData <- stats$statsTable[, 1]
}else{
	statsData <- stats$statsTable
}
