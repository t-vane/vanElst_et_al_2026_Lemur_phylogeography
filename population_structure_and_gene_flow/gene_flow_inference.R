### Script to infer rates of gene flow with the R package gene.flow.inference (https://github.com/elundgre/gene-flow-inference)

library(adegenet)
library(gene.flow.inference)
library(Matrix)
library(pegas)
library(scales)
library(vcfR)

## Set function for estimation of 1 - proportion of shared alleles (DPS)
cal.pairwise.Dps <- function(gi){
    if(is.null(pop(gi))){stop("Please assign genotypes to corresponding populations.")}
    message("Transforming data...")
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
    message("Calculating distance...")
    out.matrix <- matrix(NA, ncol = ncol(g.matrix), nrow = ncol(g.matrix))
    rownames(out.matrix) <- colnames(out.matrix) <- colnames(g.matrix)
    for(i in 1:ncol(pairs.comb)){
        print(paste0("Estimating combination ", i, " of ", length(pairs.comb)/2, " total combinations (", round(100*2*i/length(pairs.comb), digits=2),"% done...)"))
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

## Set function to prepare input data for gene flow inference
prepare_data <- function(wd, info_file, centroid_file) {
    # Set working directory
    setwd(wd)
    # Read info about individuals (format: accession, species, latitude, longitude, group)
    info <- read.table(info_file, header=TRUE)
    # Read VCF and convert to genind
    vcf <- read.vcfR("populations.snps.filtered06.mac3.vcf")
    genind <- vcfR2genind(vcf)
    # Read centroid file
    centroids <- read.table(centroid_file, header=TRUE)
    # Assign populations
    pops <- info[,1]
    genind@pop <- as.factor(pops)
    # Compute and save genetic distances
    dist <- cal.pairwise.Dps(gi = genind)
    dist <- as.matrix(dist)
    write.table(dist, file="dps.dist.txt")
    # Get group mean genetic distances and standard errors
    indiv <- lapply(1:n, FUN=function(x){which(info$groups==centroids$group[x])})
    Hs <- matrix(nrow=n,ncol=n)
    Hs_se <- matrix(nrow=n,ncol=n)
    dist_ij <- list()
    n <- length(centroids$group)
    k <- 0
    for(i in 1:n){ 
        for(j in 1:i){ 
            k <- k+1
            dist_ij[[k]] <- dist[indiv[[i]],indiv[[j]]]
            Hs[i,j] <- mean(dist_ij[[k]],na.rm=TRUE)
            Hs_se[i,j] <- sd(dist_ij[[k]],na.rm=TRUE)/sqrt(min(length(indiv[[i]]),length(indiv[[j]])))
        }
    }
    Hs[upper.tri(Hs)] <- t(Hs)[upper.tri(Hs)]
    Hs_se[upper.tri(Hs_se)] <- t(Hs_se)[upper.tri(Hs_se)]
    hs <- as.vector(Hs[upper.tri(Hs,diag=TRUE)])
    hs_se <- as.vector(Hs_se[upper.tri(Hs_se,diag=TRUE)])
    # Make adjacency matrix (sparse) and return
    G_adj <- Matrix(0, nrow=length(centroids$group), ncol=length(centroids$group), sparse=TRUE)
    G_adj <- as(as(G_adj, "generalMatrix"), "CsparseMatrix")
    # Return output
    # Return multiple objects in a list
    return(list(
        G_adj = G_adj,
        H = Hs,
        h_se = hs_se
    ))
}

## Set seed
set.seed(sample(1000000000, 1))

## Set (pre-)burn-in and iterations
preburn_iter <- 1e6
burn_iter <- 3e6
iter <- 4e6


##############################
## M. jonahi/M. macarthurii ##
##############################
## Set working directory
wd <- "/scratch/projects/nib00015/northeastProject/microcebus/jonahiMacarthurii/geneFlow"
## Get data
data <- prepare_data(wd, "jonahiMacarthurii_inds.txt", "jonahiMacarthurii_coordinates.txt")
## Assign edges of adjacency matrix
data$G_adj[1,c(2)] <- 1 # Anjiahely/Marovovonana (5) <-> Antsiatsiaka/Beanana (6)
data$G_adj[2,c(1,3,4)] <- 1 # Antsiatsiaka/Beanana (6) <-> Anjiahely/Marovovonana (5), Ambodimanga/Beanana (7), Ambolozatsy/Behovana/Antanetiambo (7)
data$G_adj[3,c(2,4,5,7)] <- 1 #A mbodimanga/Beanana (7) <-> Antsiatsiaka/Beanana (6), Ambolozatsy/Behovana/Antanetiambo (7), Antanetiambo (8), Ankoetrika (9)
data$G_adj[4,c(2,3,5)] <- 1 # Ambolozatsy/Behovana/Antanetiambo (7) <-> Antsiatsiaka/Beanana (6), Ambodimanga/Beanana (7), Antanetiambo (8)
data$G_adj[5,c(3,4,6,7)] <- 1 # Antanetiambo (8) <-> Ambodimanga/Beanana (7), Ambolozatsy/Behovana/Antanetiambo (7), Ambavala/Madera (9), Ankoetrika (9)
data$G_adj[6,c(5,7,8,9)] <- 1 # Ambavala/Madera (9) <-> Antanetiambo (8), Ankoetrika (9), Antanambe/Mananara (10), Sahafafana/Soatanana (10/11)
data$G_adj[7,c(3,5,6,9)] <- 1 # Ankoetrika (9) <-> Ambodimanga/Beanana (7), Antanetiambo (8), Ambavala/Madera (9), Sahafafana/Soatanana (10/11)
data$G_adj[8,c(6,9,10)] <- 1 # Antanambe/Mananara (10) <-> Ambavala/Madera(9), Sahafafana/Soatanana (10/11), Sasomanga (11)
data$G_adj[9,c(6,7,8,10,11)] <- 1 # Sahafafana/Soatanana (10/11) <-> Ambavala/Madera (9), Ankoetrika (9), Antanambe/Mananara (10), Sasomanga (11), Antara/Sasomanga (12)
data$G_adj[10,c(8,9,11)] <- 1 # Sasomanga (11) <-> Antanambe/Mananara (10), Sahafafana/Soatanana (10/11), Antara/Sasomanga (12)
data$G_adj[11,c(9,10,12)] <- 1 # Antara/Sasomanga (12) <-> #Sahafafana/Soatanana (10/11), Sasomanga (11), Antara/Vohirandranina (13/14)
data$G_adj[12,c(11)] <- 1 # Antara/Vohirandranina (13/14) <-> Antara/Sasomanga (12)
## Run MCMC
time <- system.time(a <- run.mcmc(G_adj_known=TRUE, G_adj=data$G_adj, g_known=FALSE, const_coal=FALSE, H=data$H, h_se=data$h_se,seed=seed, preburn_iter=preburn_iter, burn_iter=burn_iter, iter=iter, noisy_H=FALSE, type="coal"))
## Estimate medians and save output
#### estimate medians and save output
g_med <- matrixStats::colMedians(a$ans$g)
gam_med <- matrixStats::colMedians(a$ans$gam)
write.table(g_med, file=paste0(wd,"/g_med.tsv"), col.names=FALSE, row.names=FALSE)
write.table(gam_med, file=paste0(wd,"/gam_med.tsv"), col.names=FALSE, row.names=FALSE)
save(list=ls(),file=paste0(wd, "/gene_flow_inference.RData"))


######################
## M. lehilahytsara ##
######################
## Set working directory
wd <- "/scratch/projects/nib00015/northeastProject/microcebus/lehilahytsara/geneFlow"
## Get data
data <- prepare_data(wd, "lehilahytsara_inds.txt", "lehilahytsara_coordinates.txt")
## Assign edges of adjacency matrix
data$G_adj[1,c(4)] <- 1 # Marojejy (N) <-> Anjanaharibe (4)
data$G_adj[2,c(3,4)] <- 1 # Fizono (2/3) <-> Ambodivoahangy/Andaparaty (4), Anjanaharibe (4)
data$G_adj[3,c(2,4,5,6)] <- 1 # Ambodivoahangy/Andaparaty (4) <-> Fizono (2/3), Anjanaharibe (4), Anjiahely/Marovovonana (5), Antsahabe (5)
data$G_adj[4,c(1,2,3,5,6)] <- 1 # Anjanaharibe (4) <-> Marojejy (N), Fizono (2/3), Ambodivoahangy/Andaparaty (4), Anjiahely/Marovovonana (5), Antsahabe (5)
data$G_adj[5,c(3,4,6,7)] <- 1 # Anjiahely/Marovovonana (5) <-> Ambodivoahangy/Andaparaty (4), Anjanaharibe (4), Antsahabe (5), Ambalajia (6)
data$G_adj[6,c(3,4,5,7)] <- 1 # Antsahabe (5) <-> Ambodivoahangy/Andaparaty (4), Anjanaharibe (4), Anjiahely/Marovovonana (5), Ambalajia (6)
data$G_adj[7,c(5,6,8,9)] <- 1 # Ambalajia (6) <-> Anjiahely/Marovovonana (5), Antsahabe (5), Ambavala/Madera (9), Ankoetrika (9)
data$G_adj[8,c(7,9,10)] <- 1 # Ambavala/Madera (9) <-> Ambalajia (6), Ankoetrika (9), Sahafafana/Soatanana (10/11)
data$G_adj[9,c(7,8,10,11)] <- 1 # Ankoetrika (9) <-> Ambalajia (6), Ambavala/Madera (9), Sahafafana/Soatanana (10/11), Riamalandy (11/12)
data$G_adj[10,c(8,9,11)] <- 1 # Sahafafana/Soatanana (10/11) <-> Ambavala/Madera (9), Ankoetrika (9), Riamalandy (11/12)
data$G_adj[11,c(7,9,10)] <- 1 # Riamalandy (11/12) <-> Ambalajia (6), Ankoetrika (9), Sahafafana/Soatanana (10/11)
## Run MCMC
time <- system.time(a <- run.mcmc(G_adj_known=TRUE, G_adj=data$G_adj, g_known=FALSE, const_coal=FALSE, H=data$H, h_se=data$h_se,seed=seed, preburn_iter=preburn_iter, burn_iter=burn_iter, iter=iter, noisy_H=FALSE, type="coal"))
## Estimate medians and save output
#### estimate medians and save output
g_med <- matrixStats::colMedians(a$ans$g)
gam_med <- matrixStats::colMedians(a$ans$gam)
write.table(g_med, file=paste0(wd,"/g_med.tsv"), col.names=FALSE, row.names=FALSE)
write.table(gam_med, file=paste0(wd,"/gam_med.tsv"), col.names=FALSE, row.names=FALSE)
save(list=ls(),file=paste0(wd, "/gene_flow_inference.RData"))


#################
## M. simmonsi ##
#################
## Set working directory
wd <- "/scratch/projects/nib00015/northeastProject/microcebus/simmonsi/geneFlow"
## Get data
data <- prepare_data(wd, "simmonsi_inds.txt", "simmonsi_coordinates.txt")
## Assign edges of adjacency matrix
data$G_adj[1,c(2,3,4)] <- 1 # Ambodiriana <-> Ste Marie, Befotaka, Tampolo
data$G_adj[2,c(1,3,4)] <- 1 # Ste Marie <-> Ambodiriana, Befotaka, Tampolo
data$G_adj[3,c(1,2,4,5)] <- 1 # Befotaka <-> Ambodiriana, Ste Marie, Tampolo, Zahamena
data$G_adj[4,c(1,2,3,5,6)] <- 1 #Tampolo <-> Ambodiriana, Ste Marie, Befotaka, Zahamena, Betampona
data$G_adj[5,c(3,4,6)] <- 1 # Zahamena <-> Befotaka, Tampolo, Betampona
data$G_adj[6,c(4,5)] <- 1 # Betampona <-> Tampolo, Zahamena
## Run MCMC
time <- system.time(a <- run.mcmc(G_adj_known=TRUE, G_adj=data$G_adj, g_known=FALSE, const_coal=FALSE, H=data$H, h_se=data$h_se,seed=seed, preburn_iter=preburn_iter, burn_iter=burn_iter, iter=iter, noisy_H=FALSE, type="coal"))
## Estimate medians and save output
#### estimate medians and save output
g_med <- matrixStats::colMedians(a$ans$g)
gam_med <- matrixStats::colMedians(a$ans$gam)
write.table(g_med, file=paste0(wd,"/g_med.tsv"), col.names=FALSE, row.names=FALSE)
write.table(gam_med, file=paste0(wd,"/gam_med.tsv"), col.names=FALSE, row.names=FALSE)
save(list=ls(),file=paste0(wd, "/gene_flow_inference.RData"))


#############################
## A. laniger/A. mooreorum ##
#############################
## Set working directory
wd <- "/scratch/projects/nib00015/northeastProject/avahi/lanigerMooreorum/geneFlow"
## Get data
data <- prepare_data(wd, "lanigerMooreorum_inds.txt", "lanigerMooreorum_coordinates.txt")
## Assign edges of adjacency matrix
data$G_adj[1,c(2,3,7)] <- 1 #Ambavala (9) <-> Antanambe (10), Antanetiambo (8), Sahafafana (10)
data$G_adj[2,c(1,7,8)] <- 1 #Antanambe (10) <-> Ambavala (9), Sahafafana (10), Sasomanga (12)
data$G_adj[3,c(1,5)] <- 1 #Antanetiambo (8) <-> Ambavala (9), Behovana (7)
data$G_adj[4,c(5,6)] <- 1 #Beanana (6) <-> Behovana (7), Marovovonana (5)
data$G_adj[5,c(3,4)] <- 1 #Behovana (7) <-> Antanetiambo (8), Beanana (6)
data$G_adj[6,c(4,9)] <- 1 #Marovovonana (5) <-> Beanana (6), Antsahabe (4)
data$G_adj[7,c(1,2,8)] <- 1 #Sahafafana (10) <-> Ambavala (9), Antanambe (10), Sasomanga (12)
data$G_adj[8,c(2,7)] <- 1 #Sasomanga (12) <-> Antanambe (10), Sahafafana (10)
data$G_adj[9,c(6,10)] <- 1 #Antsahabe (4) <-> Marovovonana (5), Fizono (2/3)
data$G_adj[10,c(9)] <- 1 #Fizono (2/3) <-> Antsahabe (4)
## Run MCMC
time <- system.time(a <- run.mcmc(G_adj_known=TRUE, G_adj=data$G_adj, g_known=FALSE, const_coal=FALSE, H=data$H, h_se=data$h_se,seed=seed, preburn_iter=preburn_iter, burn_iter=burn_iter, iter=iter, noisy_H=FALSE, type="coal"))
## Estimate medians and save output
#### estimate medians and save output
g_med <- matrixStats::colMedians(a$ans$g)
gam_med <- matrixStats::colMedians(a$ans$gam)
write.table(g_med, file=paste0(wd,"/g_med.tsv"), col.names=FALSE, row.names=FALSE)
write.table(gam_med, file=paste0(wd,"/gam_med.tsv"), col.names=FALSE, row.names=FALSE)
save(list=ls(),file=paste0(wd, "/gene_flow_inference.RData"))








