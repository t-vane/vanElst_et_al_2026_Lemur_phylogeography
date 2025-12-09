### Script to generate Figures S12 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(corrplot)
library(gridExtra)
library(ppcor)
library(raster)
library(sf)
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
forest[is.na(forest)]<-0
forest <- ratify(forest)
RAT <- levels(forest)[[1]]
RAT$VALUE <- c("blank", "feature")
levels(forest) <- RAT

## Create cropped raster stacks for each species and plot
species_list <- c("jonahi", "simmonsi", "lehilahytsara", "laniger")
plots <- list()

for (species in species_list) {
    # Load coordinates
    coords_df <- read.table(paste0(species, "_coordinates.txt"), header = TRUE)
    coords <- coordinates(SpatialPoints(coords_df, proj4string = CRS("+init=epsg:4326")))
    # Build convex hull
    hull_idx <- chull(coords) #creates the smallest convex polygon that contains all the points
    hull_poly <- Polygon(coords[hull_idx,])
    hull <- SpatialPolygons(list(Polygons(list(hull_poly), ID="hull")))
    proj4string(hull) <- "EPSG:4326"
    # Buffer convex hull
    buffer <- raster::buffer(hull, width = 35000)
    # Crop rasters
    elevation_sp <- crop(elevation, buffer)
    rivers_sp <- crop(rivers, buffer)
    landscape_het_sp <- crop(landscape_het, buffer)
    enm_sp <- crop(get(paste0("enm_", species)), buffer)
    forest_sp <- crop(forest, buffer)
    # Stack rasters
    cov_stack <- raster::stack(list(
        elevation = raster::scale(elevation_sp),
        landscape_het = raster::scale(landscape_het_sp),
        rivers = raster::scale(rivers_sp),
        enm = raster::scale(enm_sp),
        forest = forest_sp
    ))
    # Estimate correlations
    cor <- pcor(sampleRandom(cov_stack, size= ncell(elevation) * 0.05 ), method = "spearman")
    # Generate and capture plot
    p <- recordPlot({corrplot(cor$estimate, method = "number")})
    plots[[species]] <- p
}

## Plot full
svg("Fig_S12.svg", 12, 9)
grid.arrange(grobs = plots, ncol = 2)
dev.off()
