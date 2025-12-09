#!/usr/bin/env Rscript
### Script adapted and modified from Gabriele M. Sgarlata
require(data.table)
require(ggplot2)
require(forcats)

## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_to_file <- args[1]
out_path <- args[2]

## Set function for plotting number of loci, number of polymorphic loci and number of SNPs
plot_iterating<-function(snps_per_loc){
    # Rename column 1
    colnames(snps_per_loc)[1] <- 'par_set'
    
    # Create data frame to contain number of loci and polymorphic loci
    d = snps_per_loc[,c('par_set', 'm', 'M', 'n', 'r')]
    d = d[!duplicated(d),]
  
    # Compute these numbers for each parameter set, using the par_set column as an ID
    rownames(d) = d$par_set #set rownames
    for(p in rownames(d)) {
        s = subset(snps_per_loc, par_set == p)
        d[p, 'n_loci'] = sum(s$n_loci)
        s2 = subset(s, n_snps > 0)
        d[p, 'n_loci_poly'] = sum(s2$n_loci)
        d[p, 'n_snps'] = sum(s2$n_loci*s2$n_snps)
    }
  
    # Ensure the table is ordered and plot
    d <- d[order(d$M),]
    d$r <- factor(d$r)
    plot_d <- reshape2::melt(d, id = c("par_set", "m", "M", "n", "r"))
    levels(plot_d$variable) <- c("n° assembled loci (K)", "n° polymorphic loci (K)", "n° SNPs (K)")
    return(plot_d)
}

## Plot while iterating over M with fixed n=1 and m=3
snps_per_loc <- read.delim(paste0(path_to_file, '/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, n == 1 & m == 3)
plot_M <- plot_iterating(snps_per_loc)
final.plot.M <- ggplot(plot_M, aes(x = M, y = value/1000)) +
    geom_point(aes(colour = r, shape = r)) + 
    facet_wrap(~variable, scales = "free", strip.position = "left") +
    theme_bw() +
    ylab(NULL) + 
    theme(strip.background = element_blank(), strip.placement = "outside")+
    scale_colour_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals","80% of individuals")) + 
    scale_shape_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals","80% of individuals"))
ggsave(filename = paste0(out_path, "/iterating_M_values.png"), plot = final.plot.M, width = 14, height = 5)

## Plot while iterating over m with fixed n=1 and M=2
snps_per_loc <- read.delim(paste0(path_to_file, '/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, n == 1 & M == 2)
plot_m <- plot_iterating(snps_per_loc)
final.plot.m <- ggplot(plot_m, aes(x = m, y = value/1000))+
    geom_point(aes(colour = r, shape = r)) +
    facet_wrap(~variable, scales = "free", strip.position = "left") +
    theme_bw() +
    ylab(NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside") + 
    scale_colour_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals", "80% of individuals")) +
    scale_shape_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals", "80% of individuals"))
ggsave(filename = paste0(out_path, "/iterating__m_values.png"), plot = final.plot.m, width = 14, height = 5)

## Plot while iterating over n with fixed m=3 and M=2
snps_per_loc <- read.delim(paste0(path_to_file, '/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, m == 3 & M == 2)
plot_n <- plot_iterating(snps_per_loc)
final.plot.n <- ggplot(plot_n, aes(x = n, y = value/1000)) +
    geom_point(aes(colour = r, shape = r)) + 
    facet_wrap(~variable,scales = "free", strip.position = "left") +
    theme_bw() +
    ylab(NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    scale_colour_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals","80% of individuals")) +
    scale_shape_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals","80% of individuals"))
ggsave(filename = paste0(out_path, "/iterating_n_values.png"), plot = final.plot.n, width = 14, height = 5)

## Plot all parameter combinations 
snps_per_loc <- read.delim(paste0(path_to_file, '/n_snps_per_locus.tsv'))
plot_M <- plot_iterating(snps_per_loc)
d <- plot_M[order(plot_M$m, plot_M$M, plot_M$n, plot_M$r),]
d$par_set <- gsub("-r[0-9]+", "", d$par_set)
final.plot <- gplot(d, aes(x = par_set, y = log10(value))) +
    geom_point(aes(colour = r, shape = r)) + 
    facet_wrap(~variable, scales = "free", strip.position = "left") +
    theme_bw() + 
    ylab(NULL) + 
    xlab("m & M combination") +
    aes(x = fct_inorder(par_set)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(strip.background = element_blank(), strip.placement = "outside") + 
    scale_colour_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals","80% of individuals"))+
    scale_shape_discrete(name = "", breaks = c("40", "60","80"), labels = c("40% of individuals", "60% of individuals","80% of individuals"))
ggsave(filename = paste0(out_path, "/iterating_m_M_n_comb.png"), plot = final.plot, width = 14, height = 5)
