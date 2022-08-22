# Taxa Bar plots form Plyloseq Data
# Andrew McCracken
# 4/13/22

setwd("C:/Users/andre/Desktop/Pespeni Lab/SSWD Pycnopodia 16s Networks Manuscript/Data/Taxa-bar-plots")

# # Load Libraries (idk if you need them all but here we go) 
# library(microbiome)
# library(phyloseq)
# library(qiime2R)
# library(knitr)
# library(kableExtra)
# #install.packages("remotes")
# #remotes::install_github("gmteunisse/Fantaxtic")
# library(remotes)
# library(fantaxtic)
# library(reshape2)
library(tidyverse)
library(sjmisc)
# 
# 
# import phyloseq files and create phyloseq object
physeq_all<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree.qza",
  taxonomy="taxonomy.qza",
  metadata = "pyc_manifest")

#Filter to Family Level (or Genus or Species)
filter_fam <- tax_glom(physeq_all, "Family")

#normalize per individual sample for percent abundance.
norm_abund_fam <- transform_sample_counts(filter_fam, function(x) x / sum(x))

#Filter top taxa - top 15 move abundant taxa
#https://rdrr.io/github/gmteunisse/Fantaxtic/man/get_top_taxa.html
top_tax_merge <- get_top_taxa(norm_abund_fam, 15, relative = TRUE, discard_other = FALSE, other_label = "Other")


# Phyloseq Plot
plot_bar(top_tax_merge, "animalID", "Abundance", fill="Family")


### Bar-Plots with ggplot2

########################################################
#convert to data table
source('phyloseq_to_df.R')
phylo_table_fam_all <- phyloseq_to_df(norm_abund_fam)
phylo_table_fam_15 <- phyloseq_to_df(top_tax_merge)

# write.csv(phylo_table_fam_15,"C:/Users/andre/Desktop/Pespeni Lab/SSWD Pycnopodia 16s Networks Manuscript/Taxa-bar-plots/Physeq_DF_Fam_top15", row.names = FALSE)

#phylo_table_fam_15 <- read.csv('Physeq_DF_Fam_top15')

# Randomly sample 5 from HH, 5 SH, 5 SS
# Col 9-55 -> HH
# Col 56-75 -> SH
# Col 76-93 -> SS

subset_fam_15 <- cbind(phylo_table_fam_15[,6], 
                       sample(x = phylo_table_fam_15[,9:55], size=5),
                       sample(x = phylo_table_fam_15[,56:75], size=5),
                       sample(x = phylo_table_fam_15[,76:93], size=5))


fam_15_long <- pivot_longer(subset_fam_15, 
                            col=2:16,  
                            names_to = 'Sample', 
                            values_to = 'Abundance')

names(fam_15_long)[1] <- 'taxa'

#Relative Abundance of 5 randomly sampled from each site-helath group
group_by(fam_15_long, taxa)
p <- ggplot(data=fam_15_long) + geom_bar(aes(x=Sample, y=Abundance, fill=taxa), stat="identity")

p <- p + scale_fill_manual(values=c("darkslateblue", "coral", "chocolate", "chartreuse","blue", "darkgreen", "aquamarine", "goldenrod2", "springgreen3", "firebrick3", "violet", "purple", "grey4", "gold","lightskyblue3", "tomato")) + theme(axis.text.x=element_text(angle = 90, hjust = 0)) 

p +  theme(axis.text.x=element_text(size=12))
p + theme_classic() + theme_minimal()







