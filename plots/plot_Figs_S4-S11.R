### Script to generate Figures S4-S11 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(mapdata)
library(maps)
library(raster)
library(rgdal)
library(sp)
library(terra)

## Read raster data
elevation <- raster("Elevation_GM_Res150m.tif")
rivers <- raster("rivers_padded.tif")
landscape_het <- raster("Landscape_H_GM_Res150m.tif")
enm_jonahi <- raster("ENM_jonahi_Res150m.tif")
enm_lehilahytsara <- raster("ENM_lehilahytsara_Res150m.tif")
enm_simmonsi <- raster("ENM_simmonsi_Res150m.tif")
enm_laniger <- raster("ENM_laniger_Res150m.tif")
forest <- raster("Forest1953_Res150m.tif")
forest[is.na(forest)] <- 0
forest <- ratify(forest)
RAT <- levels(forest)[[1]]
RAT$VALUE <- c("blank", "feature")
levels(forest) <- RAT

#########################
## Fig. S4 - Elevation ##
#########################
svg("Fig_S4.svg", height=8, width=6)
plot(extent(elevation), col=NA)
plot(elevation, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

#######################################
## Fig. S5 - Landscape heterogeneity ##
#######################################
svg("Fig_S5.svg", height=8, width=6)
plot(extent(landscape_het), col=NA)
plot(landscape_het, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

#################################
## Fig. S6 - Flow accumulation ##
#################################
svg("Fig_S6.svg", height=8, width=6)
plot(extent(rivers), col=NA)
## Increase river size artifically for better visualization
threshold <- 10000
river_mask <- rivers > threshold
# Create a buffer by applying a focal() on the mask
kernel <- matrix(1, 9, 9)
river_buffer <- focal(river_mask, w = kernel, fun = max)
# Fill buffer with nearby flow values (using nearest river cell)
# Distance to nearest river pixel
river_mask[river_mask == 0] <- NA
dist_to_river <- distance(river_mask)
# Assign river values only to buffered cells (interpolated using nearest river cell)
river_values <- mask(rivers, river_mask) 
buffer_fill <- focal(river_values, w = kernel, fun = max, na.policy="omit")
buffered_rivers <- ifel(river_buffer == 1, buffer_fill, NA)
# Overlay with original raster (preserve real values where available)
rivers_vis <- cover(rivers, buffered_rivers)
# Plot
plot(rivers_vis, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

###############################
## Fig. S7 - ENM (M. jonahi) ##
###############################
svg("Fig_S7.svg", height=8, width=6)
plot(extent(enm_jonahi), col=NA)
plot(enm_jonahi, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

######################################
## Fig. S8 - ENM (M. lehilahytsara) ##
######################################
svg("Fig_S8.svg", height=8, width=6)
plot(extent(enm_simmonsi), col=NA)
plot(enm_simmonsi, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

#################################
## Fig. S9 - ENM (M. simmonsi) ##
#################################
svg("Fig_S9.svg", height=8, width=6)
plot(extent(enm_lehilahytsara), col=NA)
plot(enm_lehilahytsara, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

#################################
## Fig. S10 - ENM (M. laniger) ##
#################################
svg("Fig_S10.svg", height=8, width=6)
plot(extent(enm_laniger), col=NA)
plot(enm_laniger, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

#######################
## Fig. S11 - Forest ##
#######################
svg("Fig_S11.svg", height=8, width=6)
plot(extent(forest), col=NA)
plot(forest, add=T)
map("worldHires", "Madagascar", xlim=c(47.71056,50.70222), ylim=c(-18.5,-14.10417), fill=FALSE, add=TRUE)
dev.off()

