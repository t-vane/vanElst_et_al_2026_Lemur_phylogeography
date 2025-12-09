#!/bin/sh
#!/usr/bin/env Rscript

library(radish)
library(raster)
library(terra)
library(sf)
library(sp)
library(ggplot2)
library(corrplot)
library(ppcor)

args <- commandArgs(trailingOnly=TRUE)
wd <- args[1]
prefix <- args[2]
model <- args[3]

setwd(wd)

## Report
cat("Running script with the following parameters:\n")
cat("Working directory: ", wd, "\n", sep = "")
cat("Prefix: ", prefix, "\n", sep = "")
cat("Model: ", model, "\n\n", sep = "")

## Set function to buffer river data
update_matrix <- function(mat) {
    rows <- nrow(mat)
    cols <- ncol(mat)
    updated_mat <- mat
  
    for (i in 1:rows) {
        cat(".........Processing row ", i, " of total ", rows, " (", i*100/rows, "% done)\n", sep = "")
        for (j in 1:cols) {
            a=0
            if (mat[i, j] > 1) {
                # Check for presence of top-right diagonal neighbor and absence of top and right neighbor; then update top position
                if  (i-1!=0 && j+1<=cols) {
                    if (mat[i - 1, j + 1] > 1 && mat[i - 1, j] < mat[i,j] && mat[i, j + 1] < mat[i,j] && a==0) {
                        updated_mat[i - 1, j] <- mat[i,j]
                        a=1
                    }
                }
                # Check for presence of top-left diagonal neighbor and absence of top and left neighbor; then update top position
                if (i-1!=0 && j-1!=0) {
                    if (mat[i - 1, j - 1] > 1 && mat[i - 1, j] < mat[i,j] && mat[i, j - 1] < mat[i,j] && a==0) {
                        updated_mat[i - 1, j] <- mat[i,j]
                        a=1
                    }
                }
                # Check for presence of bottom-right diagonal neighbor and absence of bottom and right neighbor; then update bottom position
                if (i+1<=rows && j+1<=cols) {
                    if (mat[i + 1, j + 1] > 1 && mat[i + 1, j] < mat[i,j] && mat[i, j + 1] < mat[i,j] && a==0) {
                        updated_mat[i + 1, j] <- mat[i,j]
                        a=1
                    }
                }
                # Check for presence of bottom-left diagonal neighbor and absence of bottom and left neighbor; then update bottom position
                if (j-1!=0 && i+1<=rows) {
                    if (mat[i + 1, j - 1] > 1 && mat[i + 1, j] < mat[i,j] && mat[i, j - 1] < mat[i,j] && a==0) {
                        updated_mat[i + 1, j] <- mat[i,j]
                        a=1
                    }
                }
            }
        }
    }
    return(updated_mat)
}

## Load geographic distances as spatial points class
cat("Creating spatial points class from coordinates\n")
coords <- paste0(prefix, "_coordinates.txt")
coordinates <- SpatialPoints(read.table(coords, header=TRUE), proj4string = CRS("EPSG:4326"))

## Load raster data
cat("Loading raster data\n")
cat("...Elevation\n")
elevation <- raster("Elevation_GM_Res150m.tif")
cat("...Rivers\n")
if (file.exists("rivers_padded.tif")) {
    rivers <- raster("rivers_padded.tif")
} else {
    rivers_raw <- raster("FlowAcc_GM_Res150m.tif")
    cat("......Buffering river raster\n")
    rivers <- update_matrix(rivers_raw)
    writeRaster(rivers, "rivers_padded.tif")
    cat("......Buffering river raster done\n")
}
cat("...Landscape heterogeneity\n")
landscape_het <- raster("Landscape_H_GM_Res150m.tif")
cat("...ENM\n")
enm <- raster(paste0("ENM_", prefix, "_Res150m.tif"))
cat("...Forest\n")
forest <- raster("Forest1953_Res150m.tif")

## Estimate extent (with buffer)
cat("Estimating extent\n")
coords <- coordinates(coordinates)
hull_idx <- chull(coords)
hull_poly <- Polygon(coords[hull_idx,])
hull <- SpatialPolygons(list(Polygons(list(hull_poly), ID = "hull")))
proj4string(hull) <- "EPSG:4326"
buffer <- raster::buffer(hull, width=35000)
coordinates@bbox <- buffer@bbox

## Crop to estimated extent and modify if necessary
cat("Cropping rasters to estimated extent\n")
elevation.c <- crop(elevation, buffer)
rivers.c <- crop(rivers, buffer)
landscape_het.c <- crop(landscape_het, buffer)
enm.c <- crop(enm, buffer)
forest.c <- crop(forest, buffer)
forest.c[is.na(forest.c)] <- 0
forest.c <- ratify(forest.c)
RAT <- levels(forest.c)[[1]]
RAT$VALUE <- c("blank", "feature")
levels(forest.c) <- RAT

## Create raster stack
cat("Creating raster stack\n")
covariates <- raster::stack(list(elevation = raster::scale(elevation.c), rivers = raster::scale(rivers.c), enm = raster::scale(enm.c),landscape_het = raster::scale(landscape_het.c), forest = forest.c))

## Load genetic distance as normal matrix
cat("Creating genetic distance matrix\n")
genetic <- read.table(paste0(prefix, "_genetic_distances_Dps.txt"), header=TRUE)
genetic <- as.matrix(genetic)

## Create conductance surface
cat("Creating conductance surface\n")
surface <- conductance_surface(covariates, coordinates, directions = 8)

## Fit radish model
cat("Fitting radish model\n")
models <- list(
    "0"   = list(formula = genetic ~ 1, desc = "null model"),
    "1"   = list(formula = genetic ~ elevation, desc = "elevation"),
    "2"   = list(formula = genetic ~ rivers, desc = "rivers"),
    "3"   = list(formula = genetic ~ enm, desc = "ENM"),
    "4"   = list(formula = genetic ~ landscape_het, desc = "landscape heterogeneity"),
    "5"   = list(formula = genetic ~ forest, desc = "forest"),
    "12"  = list(formula = genetic ~ elevation + rivers, desc = "elevation + rivers"),
    "13"  = list(formula = genetic ~ elevation + enm, desc = "elevation + ENM"),
    "14"  = list(formula = genetic ~ elevation + landscape_het, desc = "elevation + landscape heterogeneity"),
    "15"  = list(formula = genetic ~ elevation + forest, desc = "elevation + forest"),
    "23"  = list(formula = genetic ~ rivers + enm, desc = "rivers + ENM"),
    "24"  = list(formula = genetic ~ rivers + landscape_het, desc = "rivers + landscape heterogeneity"),
    "25"  = list(formula = genetic ~ rivers + forest, desc = "rivers + forest"),
    "34" = list(formula = genetic ~ enm + landscape_het, desc = "ENM + landscape heterogeneity"),
    "35" = list(formula = genetic ~ enm + forest, desc = "ENM + forest"),
    "45" = list(formula = genetic ~ landscape_het + forest, desc = "landscape heterogeneity + forest"),
    "123" = list(formula = genetic ~ elevation + rivers + enm, desc = "elevation + rivers + ENM"),
    "124" = list(formula = genetic ~ elevation + rivers + landscape_het, desc = "elevation + rivers + landscape heterogeneity"),
    "125" = list(formula = genetic ~ elevation + rivers + forest, desc = "elevation + rivers + forest"),
    "134" = list(formula = genetic ~ elevation + enm + landscape_het, desc = "ENM + landscape heterogeneity"),
    "135" = list(formula = genetic ~ elevation + enm + forest, desc = "ENM + forest"),
    "145" = list(formula = genetic ~ elevation + landscape_het + forest, desc = "elevation + landscape heterogeneity + forest"),
    "234" = list(formula = genetic ~ rivers + enm + landscape_het, desc = "rivers + ENM + landscape heterogeneity"),
    "235" = list(formula = genetic ~ rivers + enm + forest, desc = "rivers + ENM + forest"),
    "245" = list(formula = genetic ~ rivers + landscape_het + forest, desc = "rivers + landscape heterogeneity + forest"),
    "345" = list(formula = genetic ~ enm + landscape_het + forest, desc = "ENM + landscape heterogeneity + forest")
)
model_info <- models[[as.character(model)]]
desc <- model_info$desc
formula <- model_info$formula
cat("...Fitting model with", desc, "\n")
fit <- radish(
    formula,
    data = surface,
    conductance_model = radish::loglinear_conductance,
    measurement_model = radish::mlpe,
    control = NewtonRaphsonControl(maxit=2000, verbose=TRUE)
)
summary(fit)

## Save workspace
cat("Saving workspace\n")
save.image(file=paste0(prefix, "_model", model, ".RData"))

