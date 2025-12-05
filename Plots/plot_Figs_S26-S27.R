### Script to generate Figures S26-S27 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(adegenet)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pegas)
library(vcfR)
library(vegan)

## Set function to do mantel test and create data frame
run_pipeline <- function(species, set) {
    ## Read data
    geographic <- read.table(paste0(species, "_geographic_distances_", set, ".txt"), header = TRUE)
    genetic <- read.table(paste0(species, "_genetic_distances_Dps_", set, ".txt"), header = TRUE)
    genetic[is.na(genetic)] <- 0

    ## Convert to distance objects
    dist.geographic <- as.dist(as(geographic, "matrix"))
    dist.genetic    <- as.dist(as(genetic, "matrix"))

    ## Mantel test
    mantel_res <- mantel(dist.geographic, dist.genetic, method = "spearman", permutations = 9999, na.rm = TRUE)
    capture.output(mantel_res, file = paste0(species, "_", set, "mantel.txt"))

    ## Create and output data frame
    populationPairs <- c()
    geneticDistances <- c()
    geographicDistances <- c()
    for(row in 1:nrow(genetic)) {
        for(col in 1:ncol(genetic)) {
            if (row < col) {
                populationPairs <- c(populationPairs, paste0(rownames(genetic)[row], "-", colnames(genetic)[col]))
                geneticDistances <- c(geneticDistances, genetic[row, col])
                geographicDistances <- c(geographicDistances, geographic[row, col])
            }
        }
    }
    combined <- data.frame(
        geographicDistances = geographicDistances,
        geneticDistances = geneticDistances,
        populationPair = populationPairs
    )

    return(combined)
}

types <- c("ind", "pop")

##############################
## M. jonahi/M. macarthurii ##
##############################
## Get data frame and annotate for individual level
combined_ind <- run_pipeline("jonahiMacarthurii", "ind")
combined_ind[18246:18336,4] <- "macarthurii"
combined_ind$category[grepl("Mmac", combined_ind$populationPair) & grepl("Mjon", combined_ind$populationPair)] <- "between"
## Get data frame and annotate for population level
combined_pop <- run_pipeline("jonahiMacarthurii", "pop")
combined_pop[1,4] <- "macarthurii"
combined_pop[2:43,4] <- "between"

## Plot
plots_jonahiMacarthurii <- list()
for(type in types) {
    df <- get(paste0("combined_", type))
    p <- ggplot(df, aes(x = geographicDistances/1000, y = geneticDistances, fontface=1)) + 
        geom_point(aes(fill=category), size=3, pch=21) +
        xlab("Geographic distance (km)") + 
        ylab(expression(Genetic~distance~(italic(D[PS])))) + 
        theme_bw() + 
        scale_fill_manual(values=c("black", "#F8766D", "#fff00a")) +
        scale_y_continuous(limits=c(0, 0.35)) +
        scale_x_continuous(limits=c(0, 235)) +
        theme(legend.position = c(0.14, 0.897),legend.direction = "vertical",legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
            panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"),
            axis.text = element_text(size=12), axis.title=element_text(size=14)) 
    plots_jonahiMacarthurii[[type]] <- p
}

#################
## M. simmonsi ##
#################
## Get data frame and annotate for individual level
combined_ind <- run_pipeline("simmonsi", "ind")
combined_ind$category <- rep("all", nrow(combined_ind))
## Get data frame and annotate for population level
combined_pop <- run_pipeline("simmonsi", "pop")
combined_pop$category <- rep("all", nrow(combined_pop))

## Plot
plots_simmonsi <- list()
for(type in types) {
    df <- get(paste0("combined_", type))
    p <- ggplot(df, aes(x = geographicDistances/1000, y = geneticDistances, fontface=1)) + 
        geom_point(fill= "#619CFF", size=3, pch=21) + 
        xlab("Geographic distance (km)") + 
        ylab(expression(Genetic~distance~(italic(D[PS])))) + 
        theme_bw() + 
        scale_y_continuous(limits=c(0, 0.35)) +
        scale_x_continuous(limits=c(0, 235)) +
        theme(legend.position = c(0.14, 0.897),legend.direction = "vertical",legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
            panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"),
            axis.text = element_text(size=12), axis.title=element_text(size=14)) 
    plots_simmonsi[[type]] <- p
}

######################
## M. lehilahytsara ##
######################
## Get data frame and annotate for individual level
combined_ind <- run_pipeline("lehilahytsara", "ind")
combined_ind$category <- rep("all", nrow(combined_ind))
## Get data frame and annotate for population level
combined_pop <- run_pipeline("lehilahytsara", "pop")
combined_pop$category <- rep("all", nrow(combined_pop))

## Plot
plots_lehilahytsara <- list()
for(type in types) {
    df <- get(paste0("combined_", type))
    p <- ggplot(df, aes(x = geographicDistances/1000, y = geneticDistances, fontface=1)) + 
        geom_point(fill= "#00BA38", size=3, pch=21) +
        xlab("Geographic distance (km)") + ylab(expression(Genetic~distance~(italic(D[PS])))) + 
    theme_bw() + 
    scale_y_continuous(limits=c(0, 0.35)) +
    scale_x_continuous(limits=c(0, 235)) +
    theme(legend.position = c(0.14, 0.897),legend.direction = "vertical",legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
            panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"),
            axis.text = element_text(size=12), axis.title=element_text(size=14))
    plots_lehilahytsara[[type]] <- p
}


#############################
## A. laniger/A. mooreorum ##
#############################
## Get data frame and annotate for individual level
combined_ind <- run_pipeline("lanigerMooreorum", "ind")
combined_ind[494:496,4] <- "mooreorum"
combined_ind$category[grepl("Alan", combined_ind$populationPair) & grepl("Amoo", combined_ind$populationPair)] <- "between"
## Get data frame and annotate for population level
combined_pop <- run_pipeline("lanigerMooreorum", "pop")
combined_pop[1,4] <- "mooreorum"
combined_pop[2:17,4] <- "between"

## Plot
plots_lanigerMooreorum <- list()
for(type in types) {
    df <- get(paste0("combined_", type))
    p <- ggplot(df, aes(x = geographicDistances/1000, y = geneticDistances, fontface=1)) + 
        geom_point(aes(fill=category), size=3, pch=21) +
        xlab("Geographic distance (km)") + ylab(expression(Genetic~distance~(italic(D[PS])))) + 
        theme_bw() + 
        scale_fill_manual(values=c("black", "#F8766D", "#fff00a")) +
        scale_y_continuous(limits=c(0, 0.35)) +
        scale_x_continuous(limits=c(0, 235)) +
        theme(legend.position = c(0.14, 0.897),legend.direction = "vertical",legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
                panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"),
                axis.text = element_text(size=12), axis.title=element_text(size=14))
    plots_lanigerMooreorum[[type]] <- p
}

###################################
## Figure S26 - individual level ##
###################################
svg("Fig_S26.svg", 8, 8)
grid.arrange(plots_jonahiMacarthurii$ind, plots_lehilahytsara$ind, plots_simmmonsi$ind, plots_lanigerMooreorum$ind, ncol = 2)
dev.off()

###################################
## Figure S27 - population level ##
###################################
svg("Fig_S27.svg", 8, 8)
grid.arrange(plots_jonahiMacarthurii$pop, plots_lehilahytsara$pop, plots_simmmonsi$pop, plots_lanigerMooreorum$pop, ncol = 2)
dev.off()