library(tidyverse)
library(reticulate)
library(ggthemes)
library(ggprism)
library(gridExtra)

source("hsbm-analysis-helpers.R")

# mdict <- py_load_object("module-members.pickle", pickle="pickle")
fit <- py_load_object("fit-dict2.pickle", pickle="pickle")

get_block_types <- function(block) ifelse(str_detect(block, "k__"), "🦠", "⭐")

get_level_df <- function(lblocks) {
    print(lblocks)
    imap_dfr(lblocks, ~tibble(block=.y, type=get_block_types(.x), label=.x))
}

tidy_level <- function(lblocks, lvl) tibble(level=as.integer(lvl), get_level_df(lblocks))

hmem <- imap_dfr(fit$membership, tidy_level)
abun <- read_csv("level-7.csv") |>
    select(index, contains("k__"))

filter(hmem, level == 0, label == "k__Bacteria;p__Spirochaetes;c__Spirochaetes;o__Spirochaetales;f__Spirochaetaceae;g__;s__")
filter(hmem, level == 1, block == 5)$label

filter(abun, `k__Bacteria;p__Spirochaetes;c__Spirochaetes;o__Spirochaetales;f__Spirochaetaceae;g__;s__` < 18000)$index

write_csv(hmem, "hmembers-tidy.csv")

hh_subblocks <- c(25, 48, 30, 117, 22, 120, 139, 94) |> as.character()
se_subblocks <- c(131, 97, 31, 33, 8, 122, 135, 99) |> as.character()

taxa_block_order <- c(
    23, 110, 63, 51, 20, 62, 12, 29, 43, 77, 82, 130, 145, 40, 27, 121, 17, 57, 
    146, 50, 128, 49, 78, 66, 5, 96
) |>
    as.character()

## Based on this, can classify blocks 1-3 as HH1, HH2, and S
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
    # theme(
    #     strip.background=element_blank(), strip.text=element_blank(),
    #     axis.text.y = element_blank(), panel.spacing=unit(1, "mm"), axis.ticks.y=element_blank()
    # )

ggsave("level1-samples.pdf", w=3.7, h=6)
# tikz_plot(gg, "level1-samples", w=3.7, h=6)

block_pair_communities <- function(hdf, slevel=3, olevel=3, sblocks=NULL, oblocks=NULL) {
    if (is.null(sblocks))
        spred <- expression(level == slevel & type == "⭐")
    else 
        spred <- expression(level == slevel & type == "⭐" & block %in% sblocks)
    if (is.null(oblocks))
        opred <- expression(level == olevel & type == "🦠")
    else 
        opred <- expression(level == olevel & type == "🦠" & block %in% oblocks)
    
    hdf_nest <- hdf |>
        filter(eval(spred) | eval(opred)) |>
        nest(labels=c(label))

    left_join(hdf_nest, hdf_nest, by=character()) |>
        filter(block.x != block.y, type.x != type.y, type.x == "⭐") |>
        transmute(block_star=block.x, block_otu=block.y, sub_comm=map2(labels.x, labels.y, get_sub_comm))
}

get_sub_comm <- function(c1, c2) {
    filter(abun, index %in% c1$label) |>
        select(index, any_of(c2$label))
}

sub_comm_otus <- function(sub_comms) {
    labs <- map(sub_comms$sub_comm, colnames)
    map2_dfr(sub_comms$block_otu, labs, ~tibble(otu_block=.x, label=.y[2:length(.y)]))
}

sub_comm31 <- block_pair_communities(hmem, 3, 1)
    # mutate(block_otu=factor(block_otu, levels=taxa_block_order)) |>
    # arrange(block_otu)

sub_comm11 <- block_pair_communities(hmem, 1, 1) |> 
    mutate(
        block_star=fct_relevel(block_star, c(hh_subblocks, se_subblocks)),
        block_otu=factor(block_otu, levels=taxa_block_order)
    )

# l1_otus <- sub_comm_otus(sub_comm31) |>
    mutate(otu_block=factor(otu_block, levels=taxa_block_order)) |>
    arrange(otu_block)

write_csv(l1_otus, "level-1-memberships.csv")

sub_comm11 |>
    mutate(
        counts=map_dbl(sub_comm11$sub_comm, ~sum(.x[,2:ncol(.x)])),
        par_block=ifelse(block_star %in% se_subblocks, "sick exposed", "healthy")
    ) |>
    filter(counts > 0) |>
    arrange(par_block) |>
    ggplot(aes(fct_relevel(block_otu, rev(taxa_block_order)), counts, fill=fct_inorder(par_block))) +
    geom_col(position="fill", col="gray40", width=.5) +
    scale_x_discrete(expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    scale_fill_manual(values=c(`sick exposed`="#eb9e33", healthy="#6f9dcf")) +
    coord_flip() +
    labs(x=NULL, y=NULL, fill=NULL) +
    theme_bw() +
    theme(axis.ticks=element_blank(), panel.border=element_blank())

ggsave("plots/prop-participation.pdf", width=4, height=7)

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

## Matrix visualization highlighting L1 groups of interest (feel pretty meh about this)

l1samps <- hmem |> 
    filter(level == 1, type == "⭐") |> 
    mutate(block=fct_relevel(block, c(hh_subblocks, se_subblocks))) |> 
    arrange(block)

l1otus <- hmem |> 
    filter(level == 1, type == "🦠") |> 
    mutate(block=fct_relevel(block, taxa_block_order)) |> 
    arrange(block)

absort <- abun |>
    mutate(index=fct_relevel(index, rev(l1samps$label))) |>
    select(index, l1otus$label) |> 
    arrange(index) |> 
    pivot_longer(contains("k__"), "otu")
gg <- ggplot(absort, aes(otu, index, fill=log(value+1))) +
    geom_tile() +
    scale_fill_viridis_c(option="inferno") +
    theme(axis.text.x=element_blank(), axis.ticks=element_blank())

xrange <- which(l1otus$block %in% c(23, 110, 63, 51, 20))
yrange <- 86 - which(l1samps$block %in% c())

## Entropy, again (don't like this)

sub_comm21 <- block_pair_communities(hmem, 2, 1) |> 
    mutate(
        # block_star=fct_relevel(block_star, c(hh_subblocks, se_subblocks)),
        block_otu=factor(block_otu, levels=taxa_block_order)
    )

alpha_divs <- sub_comm21 |>
    mutate(
        alpha_divs=map(sub_comm, alpha_diversity, method="shannon"),
        alpha_div=map_dbl(alpha_divs, mean)
    )

samp_inter <- cumsum(c(3, 5)) + 0.5

ggplot(alpha_divs, aes(block_star, fct_rev(block_otu), fill=alpha_div)) +
    geom_tile() +
    # geom_vline(aes(xintercept=x), tibble(x=samp_inter), col="gray40", linetype="dashed") +
    scale_fill_viridis_c(option="inferno") +
    labs(x="Sample Block", y="Evenness") +
    theme_bw()

beta_divs <- sub_comm21 |>
    mutate(
        beta_divs=map(sub_comm, beta_diversity, method="shannon"),
        beta_div=map_dbl(beta_divs, mean)
    )

otu_inter <- cumsum(c(2, 9, 7, 3)) + 0.5

ggplot(beta_divs, aes(block_star, fct_rev(block_otu), fill=ifelse(beta_div == 0, NA, beta_div))) +
    geom_tile() +
    # geom_vline(aes(xintercept=x), tibble(x=samp_inter), col="white", linetype="longdash") +
    geom_hline(aes(yintercept=x), tibble(x=otu_inter), col="white", linetype="longdash") +
    scale_fill_viridis_c(option="inferno") +
    labs(x="Sample Block", y="OTU Block", fill="Evenness") +
    theme_bw()

## Entropy for each individual, labelled by group (TODO: just remove 0 pts)

sub_comm13 <- block_pair_communities(hmem, 1, 3)

alpha_divs <- sub_comm13 |>
    mutate(alpha_div=map(sub_comm, alpha_diversity, method="evenness")) |>
    unnest(alpha_div) |>
    mutate(
        block_star=fct_relevel(block_star, rev(c(hh_subblocks, se_subblocks)))
    )

samp_inter <- cumsum(c(8, 5)) + 0.5

gg1 <- ggplot(alpha_divs, aes(factor(block_star, labels=LETTERS[16:1]), alpha_div)) +
    geom_boxplot(col="#6f9dcf", alpha=0.7, outlier.shape=NA) +
    geom_point(size=0.76) +
    geom_vline(aes(xintercept=x), tibble(x=samp_inter), col="gray40", linetype="dashed") +
    coord_flip() +
    labs(x="Sample Block", y="Evenness") +
    theme_bw()

beta_divs <- sub_comm31 |>
    mutate(beta_div=map(sub_comm, beta_diversity, method="evenness")) |>
    unnest(beta_div) |>
    mutate(block_otu=fct_relevel(block_otu, rev(taxa_block_order)))

otu_inter <- cumsum(c(2, 9, 7, 3)) + 0.5

gg2 <- ggplot(beta_divs, aes(factor(block_otu, labels=length(taxa_block_order):1), beta_div)) +
    geom_boxplot(col="#eb9e33", alpha=0.7, outlier.shape=NA) +
    geom_point(size=0.76) +
    geom_vline(aes(xintercept=x), tibble(x=otu_inter), col="gray40", linetype="dashed") +
    # ylim(NA, 32) +
    coord_flip() +
    labs(x="OTU Block", y="Evenness") +
    theme_bw()

gg <- grid.arrange(gg1, gg2, nrow=1)
ggsave("plots/l1-evenness.pdf", gg, width=7, height=6)

# significance tests

oneway.test(alpha_div ~ block_star, alpha_divs, var.equal=FALSE)
oneway.test(beta_div ~ block_otu, beta_divs, var.equal=FALSE)

## Treatment contributions x entropy

prop_sick_even <- function(sc1x, sc3x) {
    evens <- sc1x |>
        mutate(
            counts=map_dbl(sc1x$sub_comm, ~sum(.x[,2:ncol(.x)]))
        ) |>
        filter(counts > 0) |>
        group_by(block_otu) |>
        mutate(prop_group=counts/sum(counts)) |>
        summarise(evenness=-sum(log(prop_group) * prop_group) / log(n()))

    sc3x |>
        bind_cols(map_dfr(sc3x$sub_comm, sub_taxa_counts)) |>
        rowwise() |>
        mutate(tot=sum(c_across(c(HH, SH, SS))), prop_healthy=HH/tot, prop_exposed=SH/tot, prop_sick=SS/tot) |>
        ungroup() |>
        left_join(evens)
}

bind_rows(
    mutate(prop_sick_even(sub_comm11, sub_comm31), level="Level 1"),
    mutate(prop_sick_even(sub_comm00, sub_comm30), level="Level 0")
) |>
    ggplot(aes(prop_sick, prop_healthy, col=evenness, label=block_otu)) +
    geom_point(size=1.5) +
    # geom_text(size=2) +
    scale_color_viridis_c(option="inferno") +
    facet_wrap(~fct_inorder(level), scales="free_y") +
    # scale_color_gradient(low="#6f9dcf", high="#eb9e33") +
    labs(x="Proportion sick", y="Evenness", col="Prop. exposed") +
    theme_bw()

ggsave("plots/prop-sick-x-evenness.pdf", width=7.2, height=4)

sub_comm11 |>
    mutate(counts=map_dbl(sub_comm11$sub_comm, ~sum(.x[,2:ncol(.x)]))) %>%
    # group_by(block_otu) %>%
    split(.$block_otu) |>
    map(counts_to_matrix)

# TODO girl.....its only one otu per group, modularity makes no sense. Thoughts:
# Max modularity in 1 vs all rest sample groups. Clearly not very useful if v present in 2 but absent elsewhere
# Max modularity over subset on sample groups. Reflects how well abunadnace can be isolated amoung the samples
# Only use evenness: one idea might be to have x axis be l 1 groups, big dot eveness, small dots evenness in each l0 child group


## Connectance

sub_comm11 |>
    mutate(
        connect=map_dbl(sub_comm, connectance)
        # block_star=factor(block_star, levels=c(hh_subblocks, se_subblocks))
    ) |>
    group_by(block_otu) |>
    summarize(mconn=mean(connect))


    ggplot(aes(block_star, connect)) +
    geom_col()

## Community diversity (attempt at entropy)

sub_comms |>
    mutate(div_healthy=map(sub_comm, ~community_diversity(filter(.x, str_detect(index, "HH"))))) |>
    mutate(div_exposed=map(sub_comm, ~community_diversity(filter(.x, str_detect(index, "SH"))))) |>
    mutate(div_sick=map(sub_comm, ~community_diversity(filter(.x, str_detect(index, "SS"))))) |>
    unnest(c("div_healthy", "div_exposed", "div_sick")) |>
    pivot_longer(contains("div"), "divs") |>
    ggplot(aes(value, col=divs)) +
    geom_density() +
    facet_wrap(~block_otu, scales="free", ncol=1)

## old family stuff

get_taxonomy <- function(labels, tax_level) {
    # if (level == "order")
    #     c <- "o__"
    c <- str_sub(tax_level, 1, 1)
    taxs <- str_extract(labels, paste0("(?<=", c, "__)\\w+(?=;)"))
    # fam <- str_extract(labels, "(?<=f__)\\w+(?=;)")
    # fam <- ifelse(
    #     is.na(fam), 
    #     str_c("Order: ", str_extract(labels, "(?<=o__)\\w+(?=;)")), 
    #     fam
    # )
    map_chr(taxs, ~if(is.na(.x)) paste0("Unknown ", tax_level) else .x)
}

top_fam_list <- function(df, lvl, by_hand) {
    top_fam <- df |>
        filter(type == "🦠", level == lvl) |>
        mutate(family=get_family(label)) |>
        count(block, family) |>
        group_by(block) |>
        mutate(pct=round(100 * n/sum(n), digits=2)) |>
        filter(pct >= 50)
        # slice_max(n, n=1) |>
        # add_count(block) |>
        # filter(nn == 1)

    top_fam
    # ret <- factor(top_fam$block, labels=top_fam$family)
    # fct_c(by_hand, ret)
}

hmem |>
    filter(type == "🦠", level == 1) |>
    mutate(order=get_taxonomy(label)) |>
    count(block, order) |>
    filter(n > 1) |>
    ggplot(aes(str_trunc(order, 10, ellipsis="."), n)) +
    geom_col() +
    facet_wrap(~fct_relevel(block, taxa_block_order), scales="free") +
    scale_x_discrete(guide=guide_axis(angle=45)) +
    theme(legend.position="none")
