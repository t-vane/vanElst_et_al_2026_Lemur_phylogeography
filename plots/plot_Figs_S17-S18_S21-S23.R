### Script to generate Figures S17-S18 and S21-S23 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(ggplot2)
library(gridExtra)
library(scales)


## Set clustering plotting function
plot_clustering <- function(individuals, populations) {
  for(i in 2:populations) {
    samples <- read.table(individuals)
    admix <- t(as.matrix(read.table(paste("populations.snps.filtered06.mac3.admix.", i, ".Q",sep=""))))
    admix <- admix[,order(samples[,2])]
    samples <- samples[order(samples[,2]),]
    niveau <- levels(as.factor(samples[,2]))
    pop <- samples[order(match(samples[,2],niveau)),]
    tempo <- tapply(1:nrow(pop),pop[,2],mean)
    tempo <- tempo[order(tempo)]
    tempomax <- tapply(1:nrow(pop),pop[,2],max)
    tempomax <- tempomax[order(tempomax)]
    barplot(admix, col=rainbow(i), space=0, border=NA, ylab="admixture", main=paste("K =",i), xaxt='n', oma=c(12, 2, 2, 2))
    text(tempo-0.5,-0.15, labels=niveau, xpd=T, srt=90, cex=.7)
    abline(v=tempomax, lty=2, lwd=1.2, col="white")
    for(j in 1:length(samples[,1])) { 
      if(length(max(admix[,j]))<2) {
        text(j-0.5,1+(0.005*which(admix[,j]==max(admix[,j]))), labels=which(admix[,j]==max(admix[,j])), xpd=T, srt=90,c ex=.5)
      }
    }
  }
}


##############
## Fig. S17 ##
##############
## M. jonahi/M. macarthurii
crossval_jonahiMacarthurii <- data.frame(
c(2,3,4,5,6,7,8,9,10,11,12,13,14),
c(0.33424, 0.28118, 0.26225, 0.26735, 0.26199, 0.26650, 0.27000, 0.28367, 0.28245,0.29234,0.29610,0.30184,0.30819))
colnames(crossval_jonahiMacarthurii) <- c("K", "error")
p_jonahiMacarthurii <- ggplot(data=crossval_jonahiMacarthurii, aes(x=K, y=error)) +
    geom_point(size=3) + 
    geom_line() +
    xlab(expression(italic(K))) + 
    ylab("Cross-validation error") +
    theme_minimal() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_continuous(breaks = seq(2, 14, by = 2)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))

## M. lehilahytsara
crossval_lehilahytsara <- data.frame(
  c(2,3,4,5,6,7,8,9,10,11,12,13,14),
  c(0.41460, 0.42487, 0.44212, 0.47307, 0.50665, 0.52355, 0.54651, 0.58764, 0.64534,0.68343,0.69273,0.70654,0.75411))
colnames(crossval_lehilahytsara) <- c("K", "error")
p_lehilahytsara <- ggplot(data=crossval_lehilahytsara, aes(x=K, y=error)) +
    geom_point(size=3) + 
    geom_line() +
    xlab(expression(italic(K))) + 
    ylab("Cross-validation error") +
    theme_minimal() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_continuous(breaks = seq(2, 14, by = 2)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))

## M. simmonsi
crossval_simmonsi <- data.frame(
  c(2,3,4,5,6,7,8),
  c(0.47417, 0.54143, 0.60883, 0.66762, 0.75582,0.85250,0.97273))
colnames(crossval_simmonsi) <- c("K", "error")
p_simmonsi <- ggplot(data=crossval_simmonsi, aes(x=K, y=error)) +
    geom_point(size=3) + 
    geom_line() +
    xlab(expression(italic(K))) + 
    ylab("Cross-validation error") +
    theme_minimal() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_continuous(breaks = seq(2, 8, by = 2)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))

## A. laniger/A. mooreorum
crossval_lanigerMooreorum <- data.frame(
  c(2,3,4,5,6,7,8,9,10,11,12,13,14),
  c(0.51762,0.58211,0.67904,0.73372,0.94580,0.97190,1.06693,1.28748,1.40974,1.28033,1.36010,1.08078,1.63659))
colnames(crossval_lanigerMooreorum) <- c("K", "error")
p_lanigerMooreorum <- ggplot(data=crossval_lanigerMooreorum, aes(x=K, y=error)) +
    geom_point(size=3) + 
    geom_line() +
    xlab(expression(italic(K))) + 
    ylab("Cross-validation error") +
    theme_minimal() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_continuous(breaks = seq(2, 14, by = 2)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))

## Plot
svg("Fig_S17.svg", 8, 8)
grid.arrange(p_jonahiMacarthurii, p_lehilahytsara, p_simmonsi, p_lanigerMooreorum, nrow = 2)
dev.off()


#########################################
## Fig. S18 - M. jonahi/M. macarthurii ##
#########################################
## Plot (labels were added manually subsequently)
pdf("Fig_S18", 12, 24)
par(mfrow=c(13,1))
plot_clustering("individuals_jonahiMacarthurii.txt", 14)
dev.off()


#################################
## Fig. S21 - M. lehilahytsara ##
#################################
## Plot (labels were added manually subsequently)
pdf("Fig_S21", 12, 24)
par(mfrow=c(13,1))
plot_clustering("individuals_lehilahytsara.txt", 14)
dev.off()


############################
## Fig. S22 - M. simmonsi ##
############################
## Plot (labels were added manually subsequently)
pdf("Fig_S23", 12, 12)
par(mfrow=c(7,1))
plot_clustering("individuals_simmonsi.txt", 8)
dev.off()


########################################
## Fig. S23 - A. laniger/A. mooreorum ##
########################################
## Plot (labels were added manually subsequently)
pdf("Fig_S23", 12, 24)
par(mfrow=c(13,1))
plot_clustering("individuals_lanigerMooreorum.txt", 14)
dev.off()

