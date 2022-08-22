# Plot KO_Brite Enrichment values

#load Libraries
library(tidyverse)

#load data

############ Raw KO Abundance BRITE3 ###################
KO_brite <- read.csv('KO_Brite_enrichment_B3.csv')

# significant Brite pathways differential between groups calculated with t-test with BH correction using Morpheus tool (https://software.broadinstitute.org/morpheus) based off of average brite enrichment accross groups after brite IDs were collected  (p<0.05)

BRITE3_sig <- read.csv('KO_BRITE3_sig.csv', header=F)

KO_brite_long <- pivot_longer(KO_brite, 
                              cols = 3:5,
                              names_to = "Health",
                              values_to = "B3_Abundance",
                              values_drop_na = TRUE)

KO_brite_long$Health <- factor(KO_brite_long$Health, levels = c("naive", "exposed", "wasting"))


ggplot(data=KO_brite_long, aes(fill=Health, x=Brite3, y=B3_Abundance)) +
  geom_bar(position="dodge", stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=10, angle=90,hjust=0.95,vjust=0.2)) +
  xlab('Brite Hierarchy') + ylab('Abundance')



# Taking subset of Brites that are differential between two groups
KO_sub <- subset(KO_brite_long, Brite3 %in% BRITE3_sig$V1)

ggplot(data=KO_sub, aes(fill=Health, x=Brite3, y=B3_Abundance)) +
  geom_bar(position="dodge", stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=12, angle=90,hjust=0.95,vjust=0.2))+
  coord_flip() + ylab("Abundance") + xlab("Level-3 Brite Catagory") +
  theme_minimal()



# Log normalize values
KO_sub_log <- mutate(KO_sub, B3_log_abun = log10(1+B3_Abundance))

ggplot(data=KO_sub_log, aes(fill=Health, x=Brite3, y=B3_log_abun)) +
  geom_bar(position="dodge", stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=10, angle=90,hjust=0.95,vjust=0.2))





##################### Average KO abundance BRITE3 ########################
brite_avg <- KO_brite %>%
  group_by(Brite3) %>%
  summarise_all(funs(mean))

avg_brite_long <- pivot_longer(brite_avg, 
                               cols = 3:5,
                               names_to = "Health",
                               values_to = "B3_Abundance",
                               values_drop_na = TRUE)

avg_brite_long$Health <- factor(avg_brite_long$Health, levels = c("naive", "exposed", "wasting"))

KO_avg_sub <- subset(avg_brite_long, Brite3 %in% BRITE3_sig$V1)


ggplot(data=KO_avg_sub, aes(fill=Health, x=Brite3, y=B3_Abundance)) +
  geom_bar(position="dodge", stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=10, angle=90,hjust=0.95,vjust=0.2)) +
  xlab('Brite Hierarchy') + ylab('Average Abundance')

#log-normalize 
avg_brite_long <- mutate(avg_brite_long, B3_log_abun = log10(1+B3_Abundance))

ggplot(data=avg_brite_long, aes(fill=Health, x=Brite3, y=B3_log_abun)) +
  geom_bar(position="dodge", stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=10, angle=90,hjust=0.95,vjust=0.2))





# ======================= Brite 4 ===========================================

KO_brite_4 <- read.csv('KO_Brite_enrichment_B4.csv')

KO_brite_long <- pivot_longer(KO_brite_4, 
                              cols = 2:4,
                              names_to = "Health",
                              values_to = "B3_Abundance",
                              values_drop_na = TRUE)

KO_brite_long$Health <- factor(KO_brite_long$Health, levels = c("naive", "exposed", "wasting"))


ggplot(data=KO_brite_long, aes(fill=Health, x=Brite4, y=B3_Abundance)) +
  geom_bar(position="dodge", stat = 'identity') +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=2, angle=90,hjust=0.95,vjust=0.2))


# Log normalize values
KO_brite_long <- mutate(KO_brite_long, B3_log_abun = log10(1+B3_Abundance))

ggplot(data=KO_brite_long, aes(fill=Health, x=Brite3, y=B3_log_abun)) +
  geom_bar(position="dodge", stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("blue", "orange", "Magenta")) +
  theme(axis.text.x=element_text(size=10, angle=90,hjust=0.95,vjust=0.2))


########## sum of each brite #################
brite4_avg <- KO_brite_4 %>%
  group_by(Brite4) %>%
  summarise_all(funs(sum))

