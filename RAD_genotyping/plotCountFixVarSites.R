#!/usr/bin/env Rscript
### Script adapted and modified from Gabriele M. Sgarlata
require(data.table)
require(forcats)
require(ggplot2)
require(reshape2)

## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_to_file <- args[1]
out_path <- args[2]

## Read data
count_fix_var <- read.delim(paste0(path_to_file, '/count_fixed_variant.tsv'))

## Plot the number of variant and fixed sites while iterating over n while keeping m=3 and M=2 fixed
count_fix_var <- subset(count_fix_var, m==3 & M==2)
count.df <- reshape2::melt(count_fix_var, id.vars = c("par_set", "m","M","n","r"), variable.name = "Sites")
final.plot.n <- ggplot(count.df, aes(x = factor(n), y = log10(value))) + 
    geom_point(aes(colour=Sites, shape=Sites)) + 
    theme_bw() + 
    labs(x="N° mismatches between Sample loci for catalog", y="log10(N° of SNPs)") +
    scale_colour_discrete(name = "", breaks = c("TotalSites", "VariantSites","FixedSites"), labels = c("All Sites", "Variant Sites","Fixed Sites")) +
    scale_shape_discrete(name = "", breaks = c("TotalSites", "VariantSites","FixedSites"), labels = c("All Sites", "Variant Sites","Fixed Sites"))
ggsave(filename = paste0(out_path, "/Count_fix_var_n_values.png"), plot = final.plot.n, width = 7, height = 5)

## Plot the number of variant sites fixed sites for all parameter combinations
count_fix_var <- read.delim(paste0(path_to_file, '/count_fixed_variant.tsv'))
colnames(count_fix_var)[1] <- 'par_set'
count_fix_var$m <- as.numeric(count_fix_var$m) 
count_fix_var$M <- as.numeric(count_fix_var$M)
count_fix_var$n <- as.numeric(count_fix_var$n)
d <- count_fix_var[order(count_fix_var$m, count_fix_var$M, count_fix_var$n),]
count.df <- reshape2::melt(d, id.vars = c("par_set", "m", "M", "n", "r"), variable.name = "Sites")  
final.plot <- ggplot(count.df, aes(x = par_set, y = log10(value))) + 
    geom_point(aes(colour = Sites, shape = Sites)) +
    theme_bw() +
    aes(x = fct_inorder(par_set)) +
    labs(x = "m-M-n combinations", y = "log10(N° of SNPs)") +
    scale_colour_discrete(name = "", breaks = c("TotalSites", "VariantSites", "FixedSites"), labels = c("All Sites", "Variant Sites","Fixed Sites")) +
    scale_shape_discrete(name = "", breaks = c("TotalSites", "VariantSites", "FixedSites"), labels = c("All Sites", "Variant Sites","Fixed Sites")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(out_path, "/Count_fix_var_m&M&n_comb.png"), plot = final.plot, width = 7, height = 5)
  