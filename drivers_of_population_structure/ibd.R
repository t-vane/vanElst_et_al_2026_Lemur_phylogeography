### Script to test for isolation-by-distance and generate Figures S26-S27 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(adegenet)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pegas)
library(vcfR)
library(vegan)

## Set function to estimate 1 - proportion of shared alleles (D_PS)
cal.pairwise.Dps <- function(gi){
    if(is.null(pop(gi))){
        stop("Please assign genotypes to corresponding populations.")
    }
    message("transforming data...")
    pops <- seppop(gi)
    g.list <- lapply(pops, FUN = function(x){
        apply(x@tab, 2, function(y){
            if(all(is.na(y))){
                return(NA)
            }else{
                out <- mean(y, na.rm = TRUE)
                return(out)
            }
        })
    })
    g.matrix <- matrix(unlist(g.list), ncol = length(g.list))
    colnames(g.matrix) <- names(g.list)
    pairs.comb <- combn(1:ncol(g.matrix), 2)
    rm(pops,g.list)
    message("calculating distance...")
    out.matrix <- matrix(NA, ncol = ncol(g.matrix), nrow = ncol(g.matrix))
    rownames(out.matrix) <- colnames(out.matrix) <- colnames(g.matrix)
    for(i in 1:ncol(pairs.comb)){
        sel1 <- pairs.comb[1,i]
        sel2 <- pairs.comb[2,i]
        p <- g.matrix[,c(sel1, sel2)]
        min.af <- apply(na.omit(p), 1, function(x){
            return(min(x))
        })
        d <- mean(min.af, na.rm = T)
        out.matrix[sel1, sel2] <- d
        out.matrix[sel2, sel1] <- d
    }
    diag(out.matrix) <- 1
    out.matrix <- 1 - out.matrix
    return(as.dist(out.matrix))
}

## Set function to do mantel test from D_PS and create data frame
run_pipeline <- function(prefix, set) {
    setwd()
    vcf_file <- paste0("/scratch/projects/nib00015/northeastProject/microcebus/", prefix, "/stacks/vcfFiltering/populations.snps.filtered06.mac3.vcf")
    sample_file <- paste0(prefix, "_pops.txt")
    
    # Read input data
    genind <- vcfR2genind(read.vcfR(vcf_file))
    samples <- read.table(sample_file, stringsAsFactors = FALSE)
    geographic <- read.table(paste0(prefix, "_geographic_distances_", set, ".txt"), header = TRUE)
    
    # D_PS estimation
    if (set == "pop") {
        pops <- samples[, 2]
    } else {
        pops <- samples[, 1]
    }
    genind@pop <- as.factor(pops)
    genetic <- as.matrix(cal.pairwise.Dps(gi = genind))
    genetic[is.na(genetic)] <- 0
    write.table(genetic, file = paste0(prefix, "_genetic_distances_Dps_", set, ".txt"), quote = FALSE, sep = "\t")

    ## Convert to distance objects
    dist.geographic <- as.dist(as(geographic, "matrix"))
    dist.genetic    <- as.dist(as(genetic, "matrix"))

    ## Mantel test
    mantel_res <- mantel(dist.geographic, dist.genetic, method = "spearman", permutations = 9999, na.rm = TRUE)
    capture.output(mantel_res, file = paste0(prefix, "_", set, "mantel.txt"))

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
combined_ind <- run_pipeline("jonahi", "ind")
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
    p <- ggplot(df, aes(x = geographicDistances/1000, y = geneticDistances, fontface = 1)) + 
        geom_point(aes(fill=category), size=3, pch=21) +
        xlab("Geographic distance (km)") + 
        ylab(expression(Genetic~distance~(italic(D[PS])))) + 
        theme_bw() + 
        scale_fill_manual(values=c("black", "#F8766D", "#fff00a")) +
        scale_y_continuous(limits=c(0, 0.35)) +
        scale_x_continuous(limits=c(0, 235)) +
        theme(legend.position = c(0.14, 0.897),legend.direction = "vertical", legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
            panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
            axis.text = element_text(size=12), axis.title=element_text(size = 14)) 
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
        theme(legend.position = c(0.14, 0.897),legend.direction = "vertical", legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
            panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
            axis.text = element_text(size=12), axis.title=element_text(size = 14)) 
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
    theme(legend.position = c(0.14, 0.897),legend.direction = "vertical", legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
            panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
            axis.text = element_text(size=12), axis.title=element_text(size = 14))
    plots_lehilahytsara[[type]] <- p
}

#############################
## A. laniger/A. mooreorum ##
#############################
## Get data frame and annotate for individual level
combined_ind <- run_pipeline("laniger", "ind")
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
        theme(legend.position = c(0.14, 0.897),legend.direction = "vertical", legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
                panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
                axis.text = element_text(size=12), axis.title=element_text(size = 14))
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