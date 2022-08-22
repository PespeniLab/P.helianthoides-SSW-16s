# 8/5/22
# Andrew McCracken
# sequence from Feautre ID and Taxon ID - Qiime

library(tidyverse)
library(phylotools)


setwd("C:/Users/andre/Desktop/Pespeni Lab/SSWD Pycnopodia 16s Networks Manuscript/data")

taxa <- read.table(file = "id_taxa_table.tsv", sep = '\t', header = TRUE)
seqs<- read.fasta(file = "repseqs.fasta", clean_name = FALSE)


# merge tables based on sequence IDs

x <-  merge(taxa, seqs, by.x = "Feature.ID", by.y = "seq.name")

#write.csv(as.data.frame(x),file="taxa_ID_Seqs.csv")

# k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio

