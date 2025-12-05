### Script to generate tree structure for Figures 1 and S13-S16 of van Elst et al. (2026), Molecular Ecology (https://doi.org/10.1111/mec.70195)
### Labels were added manually

library(ggtree)
library(ape)
library(phytools)
library(dplyr)
library(ggnewscale)
library(ggplot2)

## Read full trees
tree_raw_microcebus <- read.tree("microcebus.populations.snps.filtered06.mac3.noinv.treefile")
tree_raw_avahi <- read.tree("avahi.populations.snps.filtered06.mac3.noinv.treefile")

## Set outgroups and ladderize
outgroup_microcebus=c("Mmur_RMR49", "Mmur_RMR44_S12", "Mmur_RMR45")
outgroup_avahi=c("Aocc_f01y14_jba", "Aocc_f05y14_jba", "Aocc_f09y08_jba", "Aocc_f10y08_jba")
tree_microcebus <- ladderize(ape::root(tree_raw_microcebus, outgroup_microcebus))
tree_avahi <- ladderize(ape::root(tree_raw_avahi, outgroup_avahi))


#################################################
## Figs. 1A and S13 - M. jonahi/M. macarthurii ##
#################################################
## Subset tree and plot
keep_inds <-  c("Mjon_06Msp005_iaka_S1A","Mjon_06Msp010_iaka_S1A","Mjon_06Msp020_iaka_S1A","Mjon_06Msp024_iaka_S1A","Mjon_06Msp025_iaka_S1A","Mjon_06Msp027_iaka_S1A",
  "Mjon_06Msp029_iaka_S1A","Mjon_06Msp031_iaka_S12","Mjon_06Msp032_iaka","Mjon_06Msp033_iaka_S1A","Mjon_07Msp001_nana_S12A","Mjon_07Msp002_nana_S12A","Mjon_07Msp005_bodiNE_S12A",
  "Mjon_07Msp007_bodiNE_S12A","Mjon_07Msp008_bodiNE_S12A","Mjon_07Msp009_bodiNE","Mjon_07Msp050_atsy","Mjon_07Msp053_atsy_S12A","Mjon_07Msp058_atsy_S12","Mjon_07Msp059_atsy_S12A",
  "Mjon_07Msp061_atsy_S1A","Mjon_07Msp064_atsy_S1A","Mjon_07Msp066_atsy_S1A","Mjon_07Msp072_atsy_S1A","Mjon_07Msp073_atsy_S12A","Mjon_07Msp080_atsy","Mjon_07Msp081_atsy",
  "Mjon_07Msp083_atsy","Mjon_07Msp201_bodiNE","Mjon_07Msp202_bodiNE","Mjon_07Msp206_bodiNE","Mjon_07Msp210_bodiNE","Mjon_08Msp003_tane","Mjon_08Msp004_tane","Mjon_08Msp005_tane",
  "Mjon_08Msp006_tane","Mjon_08Msp007_tane","Mjon_08Msp008_tane","Mjon_08Msp011_tane","Mjon_08Msp012_tane_S1A","Mjon_08Msp015_tane_S1A","Mjon_08Msp017_tane_S1A","Mjon_08Msp020_tane_S1A",
  "Mjon_09Msp019_made","Mjon_09Msp030_made_S12","Mjon_09Msp036_made_S12","Mjon_09Msp045_made","Mjon_09Msp110_koet_S12","Mjon_09Msp114_koet_S12","Mjon_09Msp115_koet_S12",
  "Mjon_09Msp121_koet_S12","Mjon_10Msp009_fana","Mjon_10Msp014_fana","Mjon_10Msp016_fana","Mjon_10Msp101_fana","Mjon_10Msp106_fana","Mjon_10Msp115_fana","Mjon_10Msp120_fana",
  "Mjon_11Msp001_saso_S12","Mjon_11Msp003_saso","Mjon_11Msp005_saso","Mjon_11Msp006_saso","Mjon_11Msp008_saso","Mjon_11Msp009_saso","Mjon_11Msp011_saso","Mjon_11Msp012_saso",
  "Mjon_11Msp020_saso_S12","Mjon_11Msp022_saso","Mjon_11Msp023_saso","Mjon_11Msp024_saso","Mjon_11Msp025_saso_S12","Mjon_11Msp027_saso","Mjon_11Msp041_soa_S12","Mjon_11Msp050_fana",
  "Mjon_11Msp051_fana","Mjon_11Msp052_fana","Mjon_11Msp053_fana","Mjon_11Msp054_fana","Mjon_11Msp056_fana","Mjon_11Msp102_saso","Mjon_11Msp202_fana","Mjon_11Msp204_fana",
  "Mjon_11Msp206_fana","Mjon_11Msp207_fana","Mjon_11Msp208_fana","Mjon_12Msp001_saso_S12","Mjon_12Msp004_saso","Mjon_12Msp005_saso","Mjon_12Msp015_tara_S12","Mjon_12Msp017_tara",
  "Mjon_12Msp019_tara","Mjon_12Msp024_tara","Mjon_12Msp025_tara","Mjon_12Msp026_tara","Mjon_12Msp111_tara_S12","Mjon_12Msp112_tara","Mjon_12Msp113_tara","Mjon_12Msp116_tara",
  "Mjon_12Msp117_tara","Mjon_12Msp118_tara","Mjon_12Msp119_tara","Mjon_12Msp122_tara_S12","Mjon_13Msp001_tara_S12","Mjon_13Msp101_tara","Mjon_13Msp102_tara","Mjon_13Msp103_rand",
  "Mjon_13Msp104_tara","Mjon_13Msp105_tara","Mjon_13Msp106_rand","Mjon_13Msp107_rand","Mjon_13Msp108_rand_S12","Mjon_13Msp112_rand","Mjon_14Msp006_rand","Mjon_14Msp007_rand",
  "Mjon_14Msp016_rand","Mjon_14Msp017_rand","Mjon_14Msp019_rand","Mjon_14Msp020_rand","Mjon_14Msp022_rand_S12","Mjon_14Msp023_rand","Mjon_14Msp025_rand_S12","Mjon_14Msp100_rand_S12",
  "Mjon_14Msp103_rand","Mjon_14Msp108_rand","Mjon_14Msp111_rand","Mjon_14Msp112_rand","Mjon_14Msp113_rand","Mjon_A12_S12","Mjon_A13_S12","Mjon_A23_S12","Mjon_A24_S12","Mjon_A34_S12",
  "Mjon_AB1_S12","Mjon_B13","Mjon_B24","Mjon_B34","Mjon_BC1_S12","Mjon_BD1","Mjon_f43y21_nana","Mjon_f45y21_nana","Mjon_f46y21_nana","Mjon_f47y21_nana","Mjon_f48y21_nana_S12",
  "Mjon_f50y21_beho","Mjon_f51y21_beho","Mjon_f52y21_beho","Mjon_f53y21_beho","Mjon_f56y21_beho","Mjon_f57y21_beho","Mjon_f59y21_beho","Mjon_f61y21_beho","Mjon_f62y21_beho",
  "Mjon_f64y21_beho","Mjon_f65y21_beho","Mjon_f66y21_tane","Mjon_f68y21_tane","Mjon_f70y21_tane","Mjon_f71y21_tane","Mjon_f72y21_tane","Mjon_f75y21_tane","Mjon_f77y21_tane",
  "Mjon_f79y21_tane","Mjon_f81y21_tane","Mjon_m44y21_nana_S12","Mjon_m49y21_nana","Mjon_m63y21_beho","Mjon_m74y21_tane","Mjon_m76y21_tane","Mjon_m78y21_tane","Mjon_m80y21_tane",
  "Mjon_m82y21_tane","Mjon_MBB019","Mjon_MBB020","Mjon_MBB021","Mjon_MBB022","Mjon_MBB025","Mjon_MBB027","Mjon_MBB029","Mmac_01y06_hely_S12","Mmac_01y07_hely_S12","Mmac_04y06_hely_S12",
  "Mmac_04y13_hely_S12","Mmac_05y08_hely_S12","Mmac_06y08_hely_S12","Mmac_07y08_hely_S12","Mmac_08y08_hely_S12","Mmac_f32y21_vovo","Mmac_f35y21_vovo","Mmac_m01y06_hely",
  "Mmac_m31y21_vovo_S12","Mmac_m34y21_vovo","Mmac_m38y21_vovo","Mmac_m40y21_vovo")

svg("Fig_S13.svg", 6, 8)
p <- ggtree(keep.tip(tree_microcebus, tip = keep_inds)) + 
  geom_tree(size = 0.75) + geom_treescale() + 
  xlim_tree(0.1) + 
  geom_text2(aes(subset=(label!="100/100" & !isTip), label = label, hjust = -0.08), size=3.5)
p %>% collapse(201, mode="max") %>% collapse(221, mode="max") %>% collapse(233, mode="max") %>% collapse(238, mode="max") %>%
  collapse(256, mode="max") %>% collapse(327, mode="max") %>% collapse(322, mode="max") %>% collapse(324, mode="max") %>%
  collapse(304, mode="max") %>% collapse(318, mode="max") %>% collapse(280, mode="max") %>% collapse(332, mode="max") %>%
  collapse(356, mode="max") %>% collapse(372, mode="max")
dev.off()


#########################################
## Figs. 1B and S14 - M. lehilahytsara ##
#########################################
## Subset tree and plot
keep_inds <- c("Mmit_02y07_hely_S12","Mmit_04y07_hely_S12","Mmit_07y06_habe_S12","Mmit_08y07_hely_S12","Mmit_10y07_hely_S12","Mmit_MBB012","Mmit_MBB013","Mmit_MBB016","Mmit_PBZT115",
  "Mmit_RMR186","Mmit_RMR187","Mleh_02y00_man_S12","Mleh_04Msp001_angy_S12","Mleh_04Msp002_angy_S1A","Mleh_06Msp202_bala","Mleh_06Msp203_bala","Mleh_09Msp015_made_S12",
  "Mleh_09Msp074_made_S12","Mleh_09Msp075_made","Mleh_09Msp118_koet_S12","Mleh_09Msp123_koet","Mleh_09Msp127_koet","Mleh_10Msp018_fana_S12","Mleh_10Msp022_fana_S12","Mleh_10Msp107_fana",
  "Mleh_10Msp108_fana","Mleh_11Msp040_soa","Mleh_ANJZ11","Mleh_B12_S12","Mleh_B14_S12","Mleh_B23_S12","Mleh_BC2_S12","Mleh_BC3","Mleh_C12","Mleh_C23_S12","Mleh_C24","Mleh_DWW3235",
  "Mleh_DWW3236","Mleh_DWW3243","Mleh_DWW3244","Mleh_DWW3249","Mleh_f03y19_mana","Mleh_f04y21_anda","Mleh_f05y21_anda","Mleh_f06y06_hely_S1A","Mleh_f07y21_anda","Mleh_f08y06_habe_S12A",
  "Mleh_f08y21_anda","Mleh_f09y21_anda","Mleh_f10y21_anda","Mleh_f13y21_anda","Mleh_f15y21_fizo","Mleh_f16y21_fizo","Mleh_f17y07_hely_S12AB","Mleh_f18y21_fizo","Mleh_f19y21_fizo",
  "Mleh_f20y21_fizo","Mleh_f22y21_fizo","Mleh_f23y21_fizo","Mleh_f24y21_fizo","Mleh_f25y21_fizo","Mleh_f26y21_vovo","Mleh_f27y21_vovo","Mleh_f33y21_vovo","Mleh_f36y21_vovo",
  "Mleh_f37y21_vovo","Mleh_f39y21_vovo","Mleh_f41y21_vovo","Mleh_f42y21_vovo","Mleh_JMR001_S12","Mleh_JMR002_S12","Mleh_JMR092","Mleh_m01y08_hely","Mleh_m01y19_mana","Mleh_m01y21_anda",
  "Mleh_m02y19_mana","Mleh_m02y21_anda","Mleh_m03y21_anda","Mleh_m04y19_bohi","Mleh_m05y06_hely","Mleh_m05y19_bohi","Mleh_m06y19_bohi","Mleh_m06y21_anda","Mleh_m07y07_hely_S12A",
  "Mleh_m11y07_hely_S12A","Mleh_m11y21_anda","Mleh_m12y21_anda","Mleh_m14y07_hely_S12A","Mleh_m14y21_anda","Mleh_m15y07_hely_S12A","Mleh_m17y21_fizo","Mleh_m21y21_fizo","Mleh_m28y21_vovo",
  "Mleh_m29y21_vovo","Mleh_m30y21_vovo","Mleh_MBB001_S123","Mleh_MBB002","Mleh_MBB003_S123","Mleh_MBB037","Mleh_MBB038","Mleh_MBB039","Mleh_MBB040","Mleh_MBB041","Mleh_MBB042",
  "Mleh_MBB043","Mleh_MBB044","Mleh_MBB045","Mleh_MicroFizono09_fizo_S12","Mleh_ODY6.9","Mleh_RMR95","Mleh_RMR97","Mleh_SIB7.1","Mleh_TAD4.32")

svg("Fig_S14.svg", 6, 8)
p <- ggtree(keep.tip(tree_microcebus, tip = keep_inds)) + 
  geom_tree(size = 0.75) + 
  geom_treescale() + 
  xlim_tree(0.1) + 
  geom_text2(aes(subset=(label!="100/100" & !isTip), label = label, hjust = -0.08), size=3.5)
p1 %>% collapse(123, mode="max") %>% collapse(138, mode="max") %>% collapse(149, mode="max") %>% collapse(165, mode="max") %>%
  collapse(155, mode="max") %>% collapse(177, mode="max") %>% collapse(179, mode="max") %>% collapse(189, mode="max") %>%
  collapse(192, mode="max") %>% collapse(203, mode="max") %>% collapse(199, mode="max") %>% collapse(217, mode="max") %>%
  collapse(208, mode="max") %>% collapse(225, mode="max")
dev.off()


####################################
## Figs. 1C and S15 - M. simmonsi ##
####################################
## Subset tree and plot
keep_inds <-  c("Mbor_RMR115","Mbor_RMR116","Mbor_RMR124","Mbor_RMR129","Msim_15Msp001_taka","Msim_15Msp004_taka","Msim_15Msp005_taka","Msim_15Msp008_taka","Msim_15Msp009_taka",
  "Msim_15Msp013_taka","Msim_15Msp014_taka","Msim_15Msp016_taka","Msim_15Msp018_taka","Msim_15Msp020_taka","Msim_15Msp021_taka","Msim_15Msp101_taka","Msim_15Msp102_taka",
  "Msim_15Msp103_taka","Msim_15Msp105_taka","Msim_15Msp106_taka","Msim_15Msp107_taka","Msim_15Msp112_taka","Msim_15Msp114_taka","Msim_A01_2014","Msim_A02_2014","Msim_A03_2014",
  "Msim_A08_2014","Msim_BET87","Msim_PBZT117","Msim_POLO522","Msim_SteMsp001_stm","Msim_SteMsp002_stm_S12","Msim_SteMsp003_stm","Msim_SteMsp004_stm","Msim_ZAH2","Msim_ZAH5")

svg("Fig_S15.svg", 6, 8)
p <- ggtree(keep.tip(tree_microcebus, tip = keep_inds)) + 
geom_tree(size = 0.75) + 
geom_treescale() + 
  xlim_tree(0.1) + 
  geom_text2(aes(subset=(label!="100/100" & !isTip), label = label, hjust = -0.08), size=3.5)
p1 %>% collapse(69, mode="max") %>% collapse(62, mode="max") %>% collapse(40, mode="max") %>% collapse(58, mode="max") 
dev.off()


################################################
## Figs. 1D and S16 - A. laniger/A. mooreorum ##
################################################
## Subset tree and plot
tree_avahi_subset <- drop.tip(tree_avahi, tip = outgroup_avahi)

svg("Fig_S16.svg", 6, 8)
p <- ggtree(tree_avahi_subset) + 
  geom_tree(size = 0.75) + 
  geom_treescale() + 
  xlim_tree(0.1) + 
  geom_text2(aes(subset=(label!="100/100" & !isTip), label = label, hjust = -0.08), size=3.5)
flip(p,34,60) %>% collapse(39, mode="max") %>% collapse(45, mode="max") %>% collapse(49, mode="max") %>% collapse(53, mode="max") %>% collapse(56, mode="max") %>% 
  collapse(57, mode="max") %>% collapse(58, mode="max") %>% collapse(59, mode="max") %>% collapse(61, mode="max") %>% collapse(63, mode="max")
dev.off()
