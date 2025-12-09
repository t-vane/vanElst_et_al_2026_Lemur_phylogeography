#!/usr/bin/env Rscript
### Script adapted and modified from Gabriele M. Sgarlata
require(data.table)
require(ggplot2)

## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_to_file <- args[1]
out_path <- args[2]

## Plot SNPs per locus over n with fixed m=3, M=4 and r=0.8
snps_per_loc <- read.delim(paste0(path_to_file, '/n_snps_per_locus.tsv')) 
snps_per_loc <- subset(snps_per_loc, m == 3 & M == 2 & r == 8) 
colnames(snps_per_loc)[1] <- 'par_set' 
snps_per_loc$n <- as.factor(snps_per_loc$n) 
final.plot.n <- ggplot(snps_per_loc, aes(x = n_snps,y = n_loci)) +
    geom_point(aes(colour = n)) +
    geom_line(aes(colour = n, group = n)) +
    scale_colour_discrete() +
    theme_bw() +
    labs(x = "N° of SNPs", y = "N° of loci")
ggsave(filename = paste0(out_path, "/SNPs_distribution_by_locus.png"), plot = final.plot.n, width = 7, height = 5)

## Plot SNPs per locus for all parameter combinations with fixed r=0.8
snps_per_loc <- read.delim(paste0(path_to_file, '/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, r == 8)
colnames(snps_per_loc)[1] <- 'par_set'
final.plot <- ggplot(snps_per_loc, aes(x = n_snps, y = n_loci)) +
    geom_point(aes(colour = par_set)) +
    geom_line(aes(colour = par_set, group = par_set)) +
    scale_colour_discrete() +
    theme_bw() +
    labs(x = "N° of SNPs", y = "N° of loci") +
    guides(colour = guide_legend(title = "Parameter combinations"))
ggsave(filename = paste0(out_path, "/SNPs_distribution_by_locus_all_param_comb.png"), plot = final.plot, width = 10, height = 5)
