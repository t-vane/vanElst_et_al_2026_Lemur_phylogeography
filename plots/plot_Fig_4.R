### Script to generate Figure 4 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(raster)
library(sf)
library(sp)
library(terra)

## Read raster data
rivers <- raster("rivers_padded.tif")
landscape_het <- raster("Landscape_H_GM_Res150m.tif")
enm_jonahi <- raster("ENM_jonahi_Res150m.tif")
enm_lehilahytsara <- raster("ENM_lehilahytsara_Res150m.tif")

## Create cropped raster stacks
species_list <- c("jonahi", "lehilahytsara")
cov_stacks <- list()

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
    rivers_sp <- crop(rivers, buffer)
    landscape_het_sp <- crop(landscape_het, buffer)
    enm_sp <- crop(get(paste0("enm_", species)), buffer)
    # Stack rasters
    cov_stack <- raster::stack(list(
        landscape_het = raster::scale(landscape_het_sp),
        rivers = raster::scale(rivers_sp),
        enm = raster::scale(enm_sp)
    ))
    # Store
    cov_stacks[[species]] <- cov_stack
}

## Set function to generate plot data frame
make_plot_df <- function(raster_layer, coef, ste, ci, varname) {
    rng <- seq(cellStats(raster_layer, min), cellStats(raster_layer, max), length = 100)
    
    z <- qnorm((1 + ci) / 2)
    coef_lower <- coef - z * ste
    coef_upper <- coef + z * ste

    data.frame(
        mean  = exp(coef * rng),
        range = rng,
        lower = exp(coef_lower * rng),
        upper = exp(coef_upper * rng),
        variable = varname
    )
}

## Generate plot data frames for M. jonahi
params_jonahi <- list(
  rivers = list(layer = cov_stacks[["jonahi"]]$rivers, coef = -1.376e+01, ste  = 4.461e-05, name = "Flow accumulation"),
  enm = list(layer = cov_stacks[["jonahi"]]$enm, coef = 3.354e-01, ste  = 8.779e-06, name = "Climatic niche suitability"),
  landhet = list(layer = cov_stacks[["jonahi"]]$landscape_het, coef = 6.404e-01, ste  = 6.116e-05, name = "Landscape heterogeneity")
)
plot_df_list_jonahi <- lapply(params_jonahi, function(x) {
  make_plot_df(raster_layer = x$layer, coef = x$coef, ste = x$ste, ci = 0.95, varname = x$name)
})
plot_df_jonahi <- dplyr::bind_rows(plot_df_list_jonahi)

## Generate plot data frames for M. lehilahytsara
params_lehilahytsara <- list(
  rivers = list(layer = cov_stacks[["lehilahytsara"]]$rivers, coef = -1.376e+01, ste  = 4.461e-05, name = "Flow accumulation"),
  landhet = list(layer = cov_stacks[["lehilahytsara"]]$landscape_het, coef = 6.404e-01, ste  = 6.116e-05, name = "Landscape heterogeneity")
)
plot_df_list_lehilahytsara <- lapply(params_lehilahytsara, function(x) {
  make_plot_df(raster_layer = x$layer, coef = x$coef, ste = x$ste, ci = 0.95, varname = x$name)
})
plot_df_lehilahytsara <- dplyr::bind_rows(plot_df_list_lehilahytsara)

## Final plot
svg("Fig_4.svg", 12, 4.5)               
p_jonahi <- ggplot(plot_df_jonahi, aes(range, mean, lty = variable)) +
  geom_line(size = 1, color = "#F8766D") + 
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha = .1, linetype = 0) +
  scale_linetype_manual(values = c(3, 1, 5)) +
  xlab("Scaled variable") + 
  ylab("Estimated conductance") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 3)) +
  theme(axis.text=element_text(size = 12), axis.title=element_text(size = 14), panel.border = element_rect(colour = "black", fill = NA))
p_lehilyhtsara <- ggplot(plot_df_lehilahytsara, aes(range, mean, lty = variable)) +
  geom_line(size = 1, color = "#00BA38") + 
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha = .1, linetype = 0) +
  scale_linetype_manual(values = c(1, 5)) +
  xlab("Scaled variable") + 
  ylab("Estimated conductance") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 3)) +
  theme(axis.text=element_text(size = 12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill = NA))
grid.arrange(p_jonahi, p_lehilyhtsara, ncol=2)
dev.off()


