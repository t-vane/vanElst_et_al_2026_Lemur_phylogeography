### Script to conduct principal componant analyses and generate Figures S19 and S24 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(SNPRelate)

## Assign prefixes
prefixes <- c("jonahiMacarthurii", "jonahi", "lehilahytsara", "simmonsi", "lanigerMooreorum", "laniger")

## Loop over prefixes
for(prefix in prefixes) {
    # Read biallelic SNPs from VCF file
    vcf.fn <- paste0("/scratch/projects/nib00015/northeastProject/microcebus/", species, "/stacks/vcfFiltering/populations.snps.filtered06.mac3.vcf")
    snpgdsVCF2GDS(vcf.fn, paste0(prefix,".gds"),  method="biallelic.only")

    # Load generate geno file
    genofile <- openfn.gds(paste0(prefix,".gds"))

    # Load populations file

    # Load population codes and bind to sample IDs
    sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    pop_code <- read.table(paste0(prefix,"_pops.txt"), header = FALSE, colClasses=c("character", "character"))
    cbind(sample.id, pop_code[,2])

    # Run PCA
    pca <- snpgdsPCA(genofile, autosome.only=FALSE, remove.monosnp=FALSE)
    tab <- data.frame(sample.id = pca$sample.id,pop = factor(pop_code[,2])[match(pca$sample.id, sample.id)],
        EV1 = pca$eigenvect[,1], 
        EV2 = pca$eigenvect[,2],
        EV3 = pca$eigenvect[,3],
        EV4 = pca$eigenvect[,4],
        stringsAsFactors = FALSE)
    
    ## Assign variables
    assign(paste0("pca_", prefix), pca)
    assign(paste0("tab_", prefix), tab)
}

##############
## Fig. S19 ##
##############
## Plot (legends were added manually)
svg("Fig_S19.svg", 12, 24)
par(mfrow = c(2, 4))
# M. jonahi/M. macarthurii
plot(tab_jonahiMacarthurii$EV1, tab_jonahiMacarthurii$EV2, 
     xlab=paste0("PC1; ",round(pca$eigenval[1]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC2; ",round(pca$eigenval[2]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
plot(tab_jonahiMacarthurii$EV3, tab_jonahiMacarthurii$EV4, 
     xlab=paste0("PC3; ",round(pca$eigenval[3]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC4; ",round(pca$eigenval[4]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
# M. jonahi
plot(tab_jonahi$EV1, tab_jonahi$EV2, 
     xlab=paste0("PC1; ",round(pca$eigenval[1]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC2; ",round(pca$eigenval[2]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
plot(tab_jonahi$EV3, tab_jonahi$EV4, 
     xlab=paste0("PC3; ",round(pca$eigenval[3]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC4; ",round(pca$eigenval[4]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
# M. lehilahytsara
plot(tab_lehilahytsara$EV1, tab_lehilahytsara$EV2, 
     xlab=paste0("PC1; ",round(pca$eigenval[1]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC2; ",round(pca$eigenval[2]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
plot(tab_lehilahytsara$EV3, tab_lehilahytsara$EV4, 
     xlab=paste0("PC3; ",round(pca$eigenval[3]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC4; ",round(pca$eigenval[4]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
# M. simmonsi
plot(tab_simmonsi$EV1, tab_simmonsi$EV2, 
     xlab=paste0("PC1; ",round(pca$eigenval[1]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC2; ",round(pca$eigenval[2]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
plot(tab_simmonsi$EV3, tab_simmonsi$EV4, 
     xlab=paste0("PC3; ",round(pca$eigenval[3]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC4; ",round(pca$eigenval[4]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
dev.off()


##############
## Fig. S24 ##
##############
## Plot (legends were added manually)
svg("Fig_S24.svg", 12, 12)
par(mfrow = c(2, 2))
# A. laniger/A. mooreorum
plot(tab_lanigerMooreorum$EV1, tab_lanigerMooreorum$EV2, 
     xlab=paste0("PC1; ",round(pca$eigenval[1]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC2; ",round(pca$eigenval[2]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
plot(tab_lanigerMooreorum$EV3, tab_lanigerMooreorum$EV4, 
     xlab=paste0("PC3; ",round(pca$eigenval[3]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC4; ",round(pca$eigenval[4]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
# A. laniger
plot(tab_laniger$EV1, tab_laniger$EV2, 
     xlab=paste0("PC1; ",round(pca$eigenval[1]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC2; ",round(pca$eigenval[2]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
plot(tab_laniger$EV3, tab_laniger$EV4, 
     xlab=paste0("PC3; ",round(pca$eigenval[3]/sum(pca$eigenval[0:32])*100,2),"%"),ylab=paste0("PC4; ",round(pca$eigenval[4]/sum(pca$eigenval[0:32])*100,2),"%"),
     lwd=2, cex = 1.5, col=as.integer(tab$pop),pch=as.integer(tab$pop), cex.lab=1.5, cex.axis=1.5)
dev.off()
