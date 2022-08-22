# Pathway enrichment Bar and KEGG Analysis from Picrust2 Output
# 5/11/22

library(tidyverse)

setwd("C:/Users/andre/Desktop/Pespeni Lab/SSWD Pycnopodia 16s Networks Manuscript/picrust stuff")

#using picrust output files "path_abun_unstrat.tsv"
pi2_path <- read.table(file = 'path_abun_unstrat.tsv', sep = '\t', header = TRUE)

# Collapse KEGG pathways SAH Groups and grab mean enrichment of Pathway

pi2_avg <- transmute(pi2_path, 
               pi2_path[1],
               naive = rowMeans(pi2_path[2:48]), 
               exposed = rowMeans(pi2_path[49:68]),
               wasting = rowMeans(pi2_path[69:86]))

pi2_avg <- data.frame(pi2_avg)

#write.csv(pi2_avg, file = 'avg_path_SAH.csv', row.names=F)


#############
# collapsing KO enrichment to SAH group 
# using picrust output file "KOpred_metagenome_unstrat.tsv"
KO <- read.table(file = 'KOpred_metagenome_unstrat.tsv', sep = '\t', header = TRUE)

KO_avg <- transmute(KO, 
                     KO[1],
                     naive = rowMeans(KO[2:48]), 
                     exposed = rowMeans(KO[49:68]),
                     wasting = rowMeans(KO[69:86]))

write.csv(KO_avg, file = 'avg_KO_SAH.csv', row.names=F)

