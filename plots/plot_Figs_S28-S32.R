### Script to generate Figures S28-S32 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)

library(car)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(readxl)
library(rstatix)

## Read heterozygosity data and process
data <- read_excel("het_master.xlsx")
data$O_HET <- data$N_SITES - data$O_HOM
data$E_HET <- data$N_SITES - data$E_HOM
data$F_O_HET <- data$O_HET / data$N_SITES
data$F_E_HET <- data$E_HET / data$N_SITES
data$Ele <- as.numeric(data$Ele)

####################
## Fig. S28 - all ##
####################

## Separate M. lehilahytsara populations
data[data$Species == "Mlehilahytsara" & data$IRS == "S",]$Species <- "MlehilahytsaraS"

## Statistics
# Check for normality
model  <- lm(F_O_HET ~ Species, data = data)
# Create QQ plot of residuals
shapiro.test(residuals(model))
ggqqplot(residuals(model))
# Check for homogeneity of variances
leveneTest(F_O_HET ~ Species, data = data)
# Kruskal-Wallis test
kruskal.test(F_O_HET ~ Species, data = data)
# Post-hoc Dunn's test with Bonferroni correction
dunn_test(F_O_HET ~ Species, data=data, p.adjust.method = "bonferroni")

## Plot
svg("Fig_S28.svg", 8, 8)
ggplot(data, aes(x=Species, y=F_O_HET, fill=Species)) + 
    geom_boxplot() +
    theme_minimal() + 
    xlab("Species") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_discrete(limits = c('Mjonahi', 'Mmacarthurii', 'Mlehilahytsara', 'MlehilahytsaraS','Msimmonsi','Mmurinus')) +
    scale_fill_manual(values=c("#F8766D","#00BA38", "black", "#fff00a", "black", "#619CFF"))
dev.off()
    

#########################################
## Fig. S29 - M. jonahi/M. macarthurii ##
#########################################

## Statistics
# Check for normality
vars <- c("F_O_HET", "Lat", "Lon", "Ele")
lapply(vars, function(v)
  shapiro.test(data[data$Species == "Mjonahi", v])
)
# Check for correlations
cor.test(data$Lat[data$Species == "Mjonahi"], data$F_O_HET[data$Species == "Mjonahi"], method="spearman")
cor.test(data$Lon[data$Species == "Mjonahi"], data$F_O_HET[data$Species == "Mjonahi"], method="spearman")
cor.test(data$Ele[data$Species == "Mjonahi"], data$F_O_HET[data$Species == "Mjonahi"], method="spearman")
cor.test(data$Ele[data$Species == "Mjonahi"], data$Lat[data$Species == "Mjonahi"],  method="spearman")
cor.test(data$Ele[data$Species == "Mjonahi"], data$Lon[data$Species == "Mjonahi"],  method="spearman")

## Plot
# Lump populations
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mmac_Marovovonana_5"] <- "Mmac_5"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mmac_Anjiahely_5"] <- "Mmac_5"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Antsiatsiaka_6"] <- "Mjon_AntsBean_6"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Beanana_6"] <- "Mjon_AntsBean_6"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Beanana_7"] <- "Mjon_AmboBean_7"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Ambodimanga_7"] <- "Mjon_AmboBean_7"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Antanetiambo_7"] <- "Mjon_AntaAmboBeho_7"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Ambolozatsy_7"] <- "Mjon_AntaAmboBeho_7"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Behovana_7"] <- "Mjon_AntaAmboBeho_7"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Madera_9"] <- "Mjon_Ambavala_9"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Antanambe_10"] <- "Mjon_AntaMana_10"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Mananara_10"] <- "Mjon_AntaMana_10"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Sahafafana_10"] <- "Mjon_SahaSoat_1011"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Sahafafana_11"] <- "Mjon_SahaSoat_1011"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Soatanana_11"] <- "Mjon_SahaSoat_1011"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Antara_12"] <- "Mjon_AntaSaso_12"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Sasomanga_12"] <- "Mjon_AntaSaso_12"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Vohirandranina_14"] <- "Mjon_AntaSaso_12"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Vohirandranina_13"] <- "Mjon_AntaSaso_12"
data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population[data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",]$Population == "Mjon_Antara_13"] <- "Mjon_AntaSaso_12"

# Per locality
p1 <- ggplot(data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",], aes(x=Population, y=F_O_HET, fill=Population)) + 
    geom_boxplot() +
    theme_minimal() + 
    xlab("Locality") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_x_discrete(limits = c('Mmac_5','Mjon_AntsBean_6','Mjon_AmboBean_7','Mjon_AntaAmboBeho_7','Mjon_Antanetiambo_8','Mjon_Ankoetrika_9','Mjon_Ambavala_9','Mjon_AntaMana_10','Mjon_SahaSoat_1011','Mjon_Sasomanga_11','Mjon_AntaSaso_12')) +
    scale_fill_manual(values=c("#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D", "#fff00a"))

# Per elevation
p2 <- ggplot(data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",], aes(x=Ele, y=F_O_HET, fill=Species)) + 
    geom_point(size=3, pch=21) +
    theme_minimal() + 
    xlab("Elevation (m)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14),panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_fill_manual(values=c("#F8766D", "#fff00a"))

# Per latitude
p3 <- ggplot(data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",], aes(x=Lat, y=F_O_HET, fill=Species)) + 
    geom_point(size=3, pch=21) +
    theme_minimal() + 
    xlab("Latitude (S <-> N)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_fill_manual(values=c("#F8766D", "#fff00a"))

# Per longitude
p4 <- ggplot(data[data$Species == "Mjonahi" | data$Species=="Mmacarthurii",], aes(x=Lon, y=F_O_HET, fill=Species)) + 
    geom_point(size=3, pch=21) +
    theme_minimal() + 
    xlab("Longitude (E <-> W)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_fill_manual(values=c("#F8766D", "#fff00a"))


svg(paste0(out_dir,"suppl25_jonahi.het.svg"), 8, 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()


############################
## Fig. S30 - M. simmonsi ##
############################

## Statistics
# Check for normality
vars <- c("F_O_HET", "Lat", "Lon", "Ele")
lapply(vars, function(v)
  shapiro.test(data[data$Species == "Msimmonsi", v])
)
# Check for correlations
cor.test(data$Lat[data$Species == "Msimmonsi"], data$F_O_HET[data$Species == "Msimmonsi"], method="spearman")
cor.test(data$Lon[data$Species == "Msimmonsi"], data$F_O_HET[data$Species == "Msimmonsi"], method="spearman")
cor.test(data$Ele[data$Species == "Msimmonsi"], data$F_O_HET[data$Species == "Msimmonsi"], method="spearman")
cor.test(data$Ele[data$Species == "Msimmonsi"], data$Lat[data$Species == "Msimmonsi"],  method="spearman")
cor.test(data$Ele[data$Species == "Msimmonsi"], data$Lon[data$Species == "Msimmonsi"],  method="spearman")

## Plot
# Per locality
p1 <- ggplot(data[data$Species=="Msimmonsi",], aes(x=Population, y=F_O_HET)) + 
    geom_boxplot(fill="#619CFF") +
    theme_minimal() + 
    xlab("Locality") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_discrete(limits = c('Msim_Ambodiriana_11', 'Msim_SteMarie_11a', 'Msim_Befotaka_15', 'Msim_Tampolo_16','Msim_Zahamena_16','Msim_Betampona_16'))

# Per elevation
p2 <- ggplot(data[data$Species=="Msimmonsi",], aes(x=Ele, y=F_O_HET)) + 
    geom_point(fill= "#619CFF", size=3, pch=21) +
    theme_minimal() + 
    xlab("Elevation (m)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA))

# Per latitude
p3 <- ggplot(data[data$Species=="Msimmonsi",], aes(x=Lat, y=F_O_HET)) + 
    geom_point(fill= "#619CFF", size=3, pch=21) +
    theme_minimal() + 
    xlab("Latitude (S <-> N)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA))

# Per longitude
p4 <- ggplot(data[data$Species=="Msimmonsi",], aes(x=Lon, y=F_O_HET)) + 
    geom_point(fill= "#619CFF", size=3, pch=21) +
    theme_minimal() + 
    xlab("Longitude (W <-> E)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA))

svg("Fig_S30.svg", 8, 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()


#################################
## Fig. S31 - M. lehilahytsara ##
#################################

## Statistics
# Check for normality
vars <- c("F_O_HET", "Lat", "Lon", "Ele")
lapply(vars, function(v)
  shapiro.test(data[data$Species == "Mlehilahytsara" & data$IRS != "S", v])
)
# Check for correlations
cor.test(data$Lat[data$Species == "Mlehilahytsara" & data$IRS!="S"], data$F_O_HET[data$Species == "Mlehilahytsara" & data$IRS!="S"], method="spearman")
cor.test(data$Lon[data$Species == "Mlehilahytsara" & data$IRS!="S"], data$F_O_HET[data$Species == "Mlehilahytsara" & data$IRS!="S"], method="spearman")
cor.test(data$Ele[data$Species == "Mlehilahytsara" & data$IRS!="S"], data$F_O_HET[data$Species == "Mlehilahytsara" & data$IRS!="S"], method="spearman")
cor.test(data$Ele[data$Species == "Mlehilahytsara" & data$IRS!="S"], data$Lat[data$Species == "Mlehilahytsara" & data$IRS!="S"],  method="spearman")
cor.test(data$Ele[data$Species == "Mlehilahytsara" & data$IRS!="S"], data$Lon[data$Species == "Mlehilahytsara" & data$IRS!="S"],  method="spearman")

## Plot
# Lump populations
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Fizono_2"] <- "Mleh_Fizono_23"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Fizono_3"] <- "Mleh_Fizono_23"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Ambodivoahangy_4"] <- "Mleh_AmboAnda_4"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Andaparaty_4"] <- "Mleh_AmboAnda_4"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Anjiahely_5"] <- "Mleh_AnjiMaro_5"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Marovovonana_5"] <- "Mleh_AnjiMaro_5"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Soatanana_11"] <- "Mleh_SoatSaha_1011"
data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population[data[data$Species=="Mlehilahytsara" & data$IRS!="S",]$Population == "Mleh_Sahafafana_10"] <- "Mleh_SoatSaha_1011"

# Per location
p1 <- ggplot(data[data$Species == "Mlehilahytsara" & data$IRS!="S",], aes(x=Population, y=F_O_HET)) + 
    geom_boxplot(fill="#00BA38") +
    theme_minimal() + 
    xlab("Locality") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA)) +
    scale_x_discrete(limits = c('Mleh_Marojejy_N','Mleh_Fizono_23', 'Mleh_Anjanaharibe_4', 'Mleh_AmboAnda_4', 'Mleh_Antsahabe_5','Mleh_AnjiMaro_5', 'Mleh_Ambalajia_6', 'Mleh_Ankoetrika_9', 'Mleh_AmbavalaMadera_9', 'Mleh_SoatSaha_1011', 'Mleh_Riamalandy_S'))

# Per elevation
p2 <- ggplot(data[data$Species == "Mlehilahytsara" & data$IRS!="S",], aes(x=Ele, y=F_O_HET)) + 
    geom_point(fill= "#00BA38", size=3, pch=21) +
    theme_minimal() + 
    xlab("Elevation (m)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA))

# Lat
p3 <- ggplot(data[data$Species == "Mlehilahytsara" & data$IRS!="S",], aes(x=Lat, y=F_O_HET)) + 
    geom_point(fill= "#00BA38", size=3, pch=21) +
    theme_minimal() + 
    xlab("Latitude (S <-> N)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA))

# Lon
p4 <- ggplot(data[data$Species == "Mlehilahytsara" & data$IRS!="S",], aes(x=Lon, y=F_O_HET)) + 
    geom_point(fill= "#00BA38", size=3, pch=21) +
    theme_minimal() + 
    xlab("Longitude (W <-> E)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA))

svg("Fig_S31.svg", 8, 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()


########################################
## Fig. S32 - A. laniger/A. mooreorum ##
#######################################

## Statistics
# Check for normality
vars <- c("F_O_HET", "Lat", "Lon", "Ele")
lapply(vars, function(v)
  shapiro.test(data[data$Species == "Alaniger", v])
)
# Check for correlations
cor.test(data$Lat[data$Species == "Alaniger"], data$F_O_HET[data$Species == "Alaniger"], method="spearman")
cor.test(data$Lon[data$Species == "Alaniger"], data$F_O_HET[data$Species == "Alaniger"], method="spearman")
cor.test(data$Ele[data$Species == "Alaniger"], data$F_O_HET[data$Species == "Alaniger"], method="spearman")
cor.test(data$Ele[data$Species == "Alaniger"], data$Lat[data$Species == "Alaniger"],  method="spearman")
cor.test(data$Ele[data$Species == "Alaniger"], data$Lon[data$Species == "Alaniger"],  method="spearman")

## Plot
# Per locality
p1 <- ggplot(data[data$Species == "Alaniger" | data$Species=="Amooreorum",], aes(x=Population, y=F_O_HET, fill=Population)) + 
    geom_boxplot() +
    theme_minimal() + 
    xlab("Locality") + ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_x_discrete(limits = c('Amoo_Fizono_2', 'Amoo_Antsahabe_4', 'Alan_Marovovonana_5', 'Alan_Beanana_6', 'Alan_Behovana_7', 'Alan_Antanetiambo_8', 'Alan_Ambavala_9', 'Alan_Sahafafana_10', 'Alan_Antanambe_10', 'Alan_Sasomanga_12')) +
    scale_fill_manual(values=c("#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#fff00a","#fff00a"))

# Per elevation
p2 <- ggplot(data[data$Species == "Alaniger" | data$Species=="Amooreorum",], aes(x=Ele, y=F_O_HET, fill=Species)) + 
    geom_point(size=3, pch=21) +
    theme_minimal() + 
    xlab("Elevation (m)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_fill_manual(values=c("#F8766D", "#fff00a"))

# Per latitude
p3 <- ggplot(data[data$Species == "Alaniger" | data$Species=="Amooreorum",], aes(x=Lat, y=F_O_HET, fill=Species)) + 
    geom_point(size=3, pch=21) +
    theme_minimal() + 
    xlab("Latitude (S <-> N)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_fill_manual(values=c("#F8766D", "#fff00a"))

# Per longitude
p4 <- ggplot(data[data$Species == "Alaniger" | data$Species=="Amooreorum",], aes(x=Lon, y=F_O_HET, fill=Species)) + 
    geom_point(size=3, pch=21) +
    theme_minimal() + 
    xlab("Longitude (E <-> W)") + 
    ylab("Observed heterozygosity") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA), legend.position="none") +
    scale_fill_manual(values=c("#F8766D", "#fff00a"))

svg("Fig_S32.svg", 8, 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
