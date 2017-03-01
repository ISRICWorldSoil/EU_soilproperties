## LUCAS points cleaned up by Lesschen, Jan Peter <janpeter.lesschen@wur.nl>
## about 20,000 points

setwd("/data/Integrator/points")
library(foreign)
library(sp)
library(gstat)
library(rgdal)
library(lattice)
library(maptools)
library(RColorBrewer)
library(plotKML)

LUCAS = read.csv("Lucas_soil_and_land_use_data_2009.csv", na.strings = c("","#NA","#N/A","NA"))
str(LUCAS)

## ID column:
plyr:::nunique(LUCAS$POINT_ID)
nrow(LUCAS)
summary(LUCAS$LU1)

## rename columns:
LUCAS <- plyr::rename(LUCAS, c("POINT_ID"="SOURCEID", "sand"="SNDPPT", "clay"="CLYPPT", "silt"="SLTPPT", "pH_in_H2O"="PHIHOX", "pH_in_CaCl"="PHICAL", "OC"="ORCDRC", "GPS_LAT"="LATWGS84", "GPS_LONG"="LONWGS84", "CEC"="CECSUM"))
str(LUCAS)
summary(LUCAS$ORCDRC)
summary(LUCAS$PHICAL)

LUCAS$SAMPLEID <- make.unique(paste0("LUCAS_LC", LUCAS$SOURCEID))
LUCAS$SOURCEDB = "LUCAS_LC"
LUCAS$TIMESTRR <- as.Date("2009", format="%Y")
LUCAS$UHDICM = 0
LUCAS$LHDICM = 20
LUCAS$DEPTH <- LUCAS$UHDICM + (LUCAS$LHDICM - LUCAS$UHDICM)/2
SPROPS.LUCAS_LC <- LUCAS[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIHOX","PHICAL","ORCDRC","CECSUM","LC1")]
str(SPROPS.LUCAS_LC)
## 19,860
save(SPROPS.LUCAS_LC, file="SPROPS.LUCAS_LC.rda")
plot(SPROPS.LUCAS_LC$LONWGS84, SPROPS.LUCAS_LC$LATWGS84, pch="+")

