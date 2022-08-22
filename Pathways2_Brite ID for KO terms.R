# Getting Brite catagories for KO terms
# 5/17/22
# Andrew McCracken


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# Need to delete the LOCK00 Folder in R for proper download -- Folow error instructions
# BiocManager::install("KEGGREST")

library(KEGGREST)
library(tidyverse)

KO <- read.csv('avg_KO_SAH.csv') #Picrust output
KO <- mutate(KO, Brite3 = NA)

ko_list <- KO$function.
#ko_list <- head(KO$function.)

#### get URL error on K00050
# for (i in 149:length(ko_list)){
#   x <- keggGet(KO$function.[i])
#   x <- unlist(x)
#   x <- as.list(x)
#   KO$Brite3[i] <- x$BRITE3
# }


############ FINISHED  ~ 2h runtime ~ ########################
# BRITE3 : avoiding url error 
 for (i in 1:length(ko_list)){
   x <- tryCatch(keggGet(KO$function.[i]), error=function(e) NA)
   if(is.na(x) == TRUE){
     next
   }
   else{
     x <- unlist(x)
     x <- as.list(x)
     KO$Brite3[i] <- x$BRITE3
   }
 }

write.csv(KO, file = 'KO_Brite_enrichment.csv', row.names = F)

########################################################################
# # Brite 4: 

KO <- mutate(KO, Brite4 = NA)
for (i in 1:length(ko_list)){
  x <- tryCatch(keggGet(KO$function.[i]), error=function(e) NA)
  if(is.na(x) == TRUE){
    next
  }
  else{
    x <- unlist(x)
    x <- as.list(x)
    KO$Brite4[i] <- x$BRITE4
  }
}
write.csv(KO, file = 'KO_Brite_enrichment_B4.csv', row.names = F)

##############################################################

### Check individual incidences
x <- tryCatch(keggGet(KO$function.[41]), error=function(e) NA)
x <- keggGet(KO$function.[41])
x <- unlist(x)
x <- as.list(x)
KO$Brite3[41] <- x$BRITE3

### Check Individual KO terms
x <- keggGet('K00356')
x <- unlist(x)
x <- as.list(x)
x$BRITE3

