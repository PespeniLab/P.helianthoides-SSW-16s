library(tidyverse)
library(reticulate)
library(ggthemes)
library(ggprism)
library(gridExtra)

source("hsbm-analysis-helpers.R")

# Load the results dictionary and save it as a dataframe
fit <- py_load_object("fit-dict.pickle", pickle="pickle")

hmem <- imap_dfr(fit$membership, tidy_level)
abun <- read_csv("level-7.csv") |>
    select(index, contains("k__"))

write_csv(hmem, "hmembers-tidy.csv")

# Level 1 seastar memberships: which types of animals were grouped together?------
hh_subblocks <- c(25, 48, 30, 117, 22, 120, 139, 94) |> as.character()
se_subblocks <- c(131, 97, 31, 33, 8, 122, 135, 99) |> as.character()

taxa_block_order <- c(
    23, 110, 63, 51, 20, 62, 12, 29, 43, 77, 82, 130, 145, 40, 27, 121, 17, 57, 
    146, 50, 128, 49, 78, 66, 5, 96
) |>
    as.character()

hmem |>
    filter(type == "⭐", level == 1) |>
    mutate(treatment = factor(str_extract(label, "[A-Z]+"), labels=c("Healthy", "Exposed", "Sick"))) |>
    # group_by(level, block) |>
    count(block, treatment) |>
    mutate(block=factor(block, levels=c(rev(se_subblocks), rev(hh_subblocks)))) |>
    ggplot(aes(block, n, fill=fct_rev(treatment))) +
    geom_col(position=position_stack(), width=0.5, col="gray40") +
    geom_vline(xintercept=8.5, col="gray40", linetype="dashed") +
    geom_vline(xintercept=13.5, col="gray40", linetype="dashed") +
    # facet_wrap(~block, scales="fixed", ncol=1) +
    scale_fill_manual(values=c(Healthy="#6f9dcf", Exposed="#eb9e33", Sick="#bf3f48e6")) +
    scale_y_continuous(guide = guide_prism_minor(), breaks=c(0, 3, 6, 9), minor_breaks=1:9) +
    coord_flip() +
    labs(x="Sample block", y="Num. samples", fill=NULL) +
    theme_bw() +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

ggsave("plots/level1-samples.pdf", w=3.7, h=6)

# Level 1 community contributions-------------------------------------------------
# which seastar types did microbial communities appear on?------------------------
sub_comm31 <- block_pair_communities(hmem, 3, 1)

sub_comm31 |>
    bind_cols(map_dfr(sub_comm31$sub_comm, sub_taxa_counts)) |>
    pivot_longer(c(HH, SH, SS), "treatment", values_to="counts") |>
    ggplot(aes(fct_relevel(block_otu, rev(taxa_block_order)), counts, fill=fct_rev(treatment))) +
    geom_col(position="fill", width=.5, col="gray60") +
    scale_x_discrete(expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    scale_fill_manual(values=c(SS="#bf3f48e6", SH="#eb9e33", HH="#6f9dcf")) +
    coord_flip() +
    labs(x=NULL, y=NULL, fill=NULL) +
    theme_bw() +
    theme(axis.ticks=element_blank(), panel.border=element_blank())

ggsave("plots/contributions.pdf", width=4, height=7)

# Entropy analysis: identify subcommunities driving alpha and beta diversity------
sub_comm13 <- block_pair_communities(hmem, 1, 3)

alpha_divs <- sub_comm13 |>
    mutate(alpha_div=map(sub_comm, alpha_diversity, method="shannon")) |>
    unnest(alpha_div) |>
    mutate(
        block_star=fct_relevel(block_star, rev(c(hh_subblocks, se_subblocks)))
    )

samp_inter <- cumsum(c(8, 5)) + 0.5

ggplot(alpha_divs, aes(factor(block_star, labels=LETTERS[16:1]), alpha_div)) +
    geom_boxplot(col="#6f9dcf", alpha=0.7, outlier.shape=NA) +
    geom_point(size=0.76) +
    geom_vline(aes(xintercept=x), tibble(x=samp_inter), col="gray40", linetype="dashed") +
    coord_flip() +
    labs(x="Sample Block", y="Entropy") +
    theme_bw()

beta_divs <- sub_comm31 |>
    mutate(beta_div=map(sub_comm, beta_diversity, method="shannon")) |>
    unnest(beta_div) |>
    mutate(block_otu=fct_relevel(block_otu, rev(taxa_block_order)))

otu_inter <- cumsum(c(2, 9, 7, 3)) + 0.5

ggplot(beta_divs, aes(factor(block_otu, labels=length(taxa_block_order):1), beta_div)) +
    geom_boxplot(col="#eb9e33", alpha=0.7, outlier.shape=NA) +
    geom_point(size=0.76) +
    geom_vline(aes(xintercept=x), tibble(x=otu_inter), col="gray40", linetype="dashed") +
    coord_flip() +
    labs(x="OTU Block", y="Entropy") +
    theme_bw()

# significance tests
oneway.test(alpha_div ~ block_star, alpha_divs, var.equal=FALSE)
oneway.test(beta_div ~ block_otu, beta_divs, var.equal=FALSE)
