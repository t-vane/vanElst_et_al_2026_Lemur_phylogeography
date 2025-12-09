### Script to generate Figures S1-S2 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(ggplot2)
library(mapdata)
library(mapplots)
library(maps)
library(maptools)
library(raster)
library(rgdal)
library(rgeos)
library(scales)
library(sp)

## Load raster data
madarivers_major_shapes_t = spTransform(readOGR("River_Mada_1"), "+init=epsg:4326")
madarivers_major_df_t <- fortify(madarivers_major_shapes_t)
DEM <- raster("DEM_JAXA_WGS84.tif")


##########################
## Fig. S1 - Microcebus ##
##########################
## Read sample data
samples <- read.table("microcebus_sampling.txt", header=TRUE)

## Assign colors and shapes
cols <- samples$spec
cols[cols == "jonahi"] <- "#F8766D"
cols[cols == "macarthurii"] <- "orange"
cols[cols == "lehilahytsara"] <- "#00BA38"
cols[cols == "simmonsi"] <- "#619CFF"

shps <- samples$status
shps[shps == "new"] <- 21
shps[shps == "old"] <- 22
shps <- as.double(shps)

## Jitter coordinates
samples$lon <- jitter(samples$lon, 100)
samples$lat <- jitter(samples$lat, 100)

## Plot (labels were added subsequently)
svg("Fig_S1.svg", 6, 8)
map("worldHires","Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), col="#ffffbf", fill=TRUE)
axis(1, at=c(48,48.5,49,49.5,50,50.5), labels=c("48", "","49","","50",""))
axis(2, at=c(-18.5,-18,-17.5,-17,-16.5,-16,-15.5,-15,-14.5), labels=c("","-18","","-17","","-16","","-15",""))
col <- c("#ffefb0","#fadb9d","#ebc291","#cfa288","#b58381")
brk <- c(250, 400, 600, 800, 1000, 4000)
plot(DEM, col=col, breaks=brk, main="", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), legend=FALSE, border=FALSE, add=TRUE, axes=FALSE)
plot(madarivers_major_shapes_t, xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), col=alpha("royalblue",1), border=FALSE, add=TRUE)
points(samples$lon,samples$lat, col="black", bg=cols, pch = shps, cex=1.5)
dev.off()


#####################
## Fig. S2 - Avahi ##
#####################
## Read sample data
samples <- read.table("avahi_sampling.txt", header=TRUE)

## Assign colors and shapes
cols[cols == "laniger"] <- "#F8766D"
cols[cols == "mooreorum"] <- "#00BA38"

shps <- samples$status
shps[shps == "new"] <- 21
shps[shps == "old"] <- 22
shps <- as.double(shps)

## Jitter coordinates
samples$lon <- jitter(samples$lon, 100)
samples$lat <- jitter(samples$lat, 100)

## Plot (labels were added subsequently)
svg("Fig_S2.svg", 6, 8)
map("worldHires","Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), col="#ffffbf", fill=TRUE)
axis(1, at=c(48,48.5,49,49.5,50,50.5), labels=c("48","","49","","50",""))
axis(2, at=c(-18.5,-18,-17.5,-17,-16.5,-16,-15.5,-15,-14.5), labels=c("","-18","","-17","","-16","","-15",""))
col <- c("#ffefb0","#fadb9d","#ebc291","#cfa288","#b58381")
brk <- c(250, 400, 600, 800, 1000, 4000)
plot(DEM, col=col, breaks=brk, main="",xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), legend=FALSE, border=FALSE, add=TRUE,axes=FALSE)
plot(madarivers_major_shapes_t, xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), col=alpha("royalblue",1), border=FALSE, add=TRUE)
points(samples$lon,samples$lat, col="black", bg=cols, pch = shps, cex=1.5)
dev.off()