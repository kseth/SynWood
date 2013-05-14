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

startingInfested <- which(maps$OLDSTATUS == 1)
endInfestedHouses <- which(maps$STATUS == 1)
binomialEndInfested <- rep(0, length(maps$X))
binomialEndInfested[endInfestedHouses] <- 1
