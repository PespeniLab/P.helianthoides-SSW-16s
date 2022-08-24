# ANCOMBC 
# 4/19/22

library(microbiome)
library(phyloseq)
library(qiime2R)
library(ANCOMBC)
library(knitr)
library(kableExtra)

#Working Directory
getwd()
setwd("C:/Users/andre/Desktop/Pespeni Lab/SSWD Pycnopodia 16s Networks Manuscript")

#Use table.qza for differential ASVs and table-lev7.qza for ASVs collapse to species
physeq<-qza_to_phyloseq(
  features="table-lev7.qza",
  tree=,
  taxonomy=,
  metadata = "pyc_manifest.txt")

# Read in manifest
Manifest <- read.table("pyc_manifest.txt", header = TRUE, row.names = 1)


### Subsisting the files for pairwise analysis of differential abundance ###

#Subset samples, just healthy animals (inclusive of Naive and Exposed samples)
pseq_healthy <- subset_samples(physeq, animal.health == "healthy")

#Subset samples, just impacted sites (inclusive of Exposed and Wasting samples)
pseq_impacted <- subset_samples(physeq, site.status == "impacted")



#ANCOMBC for site-animal-health Naive vs Exposed (Taxa-lev7 : Species)
#Compared by Site-animal-health status as (HH=Naive, SH=Exposed)
outHHtoSH <- ancombc(
  phyloseq = pseq_healthy, 
  group = "site.animal.health",
  formula = "site.animal.health", #dots not hyphens for file name?
  p_adj_method = "fdr", 
  zero_cut = 0.9, # by default prevalence filter of 10% is applied, this causes us to lose too many taxa that aren't in many of the samples
  lib_cut = 0,
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
resHHtoSH <- as.data.frame(outHHtoSH$res)
colnames(resHHtoSH) <- c("Exp_beta", "Exp_se", "Exp_W", "Exp_p_val", "Exp_q_val", "Exp_DA")


write.csv(as.data.frame(resHHtoSH),file="resHHtoSH_lev7.csv")


####################################################################

#ANCOMBC for site-animal-health Exposed vs Wasting (Taxa-lev7 : Species)
#Compared by Site-animal-health status as (SH=Exposed, SH=Wasting)
outSHtoSS <- ancombc(
  phyloseq = pseq_impacted, 
  group = "site.animal.health",
  formula = "site.animal.health", #dots not hyphens for file name?
  p_adj_method = "fdr", 
  zero_cut = 0.9, # by default prevalence filter of 10% is applied, this causes us to lose too many taxa that aren't in many of the samples
  lib_cut = 0,
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
resSHtoSS <- as.data.frame(outSHtoSS$res)
colnames(resSHtoSS) <- c("Waste_beta", "Waste_se", "Waste_W", "Waste_p_val", "Waste_q_val", "Waste_DA")

write.csv(as.data.frame(resSHtoSS),file="resSHtoSS_lev7.csv")


###Merge ANCOM-BC lists and write to CSV
Both <- merge(resHHtoSH, resSHtoSS, by = 0, all.x = TRUE, all.y = TRUE) 
write.csv(Both, "ANCOMBC_Lev7_All.csv")



#================= Ven-Diagrams =======================

#Read in Differential abundant (species) lists

Up_SH <- resHHtoSH %>% subset(resHHtoSH$Exp_q_val< 0.05 & resHHtoSH$Exp_beta > 0)

Down_SH <- resHHtoSH %>% subset(resHHtoSH$Exp_q_val < 0.05 & resHHtoSH$Exp_beta < 0)

Up_SS <- resSHtoSS %>% subset(resSHtoSS$Waste_q_val < 0.05 & resSHtoSS$Waste_beta > 0)
Down_SS <- resSHtoSS %>% subset(resSHtoSS$Waste_q_val < 0.05 & resSHtoSS$Waste_beta < 0)

List <- list("Up: Naive vs Exposed"=rownames(Up_SH), "Down: Naive vs Exposure"=rownames(Down_SH),
             "Up: Exposed vs Wasting"=rownames(Up_SS), "Down: Exposed vs Wasting"=rownames(Down_SS))

### Make Euler diagram
library(eulerr)
Eu <- euler(List, shape = "ellipse")

plot(Eu, fills = list(fill = c("Orange", "Orange3", "magenta", "magenta"), alpha = 1),
     shape = "ellipse",)
