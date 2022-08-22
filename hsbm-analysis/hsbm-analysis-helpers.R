library(tidyverse)
library(tikzDevice)

hsbm_summary <- function(comm, lvl=0) {
    ret <- comm |>
        filter(level == lvl) |>
        group_by(type)

    ret |>
        count(block) |>
        summarise(mean_block_size=mean(n), max_block_size=max(n), min_block_size=min(n)) |>
        right_join(summarise(ret, samp_total=n(), nblocks=length(unique(block))))
}

block_members <- function(tidy_member_df, block, level=1, tax_level=NULL) {
    ret <- filter(tidy_member_df, level == !!level, block == !!block) |>
        pull(label)
    if (!is.null(tax_level))
        ret <- get_taxonomy(ret, tax_level)
    ret
}

counts_to_matrix <- function(comm_counts) {
    samp <- unique(comm_counts$block_star)
    otu <- unique(comm_counts$block_otu)
    ret <- matrix(comm_counts$counts, nrow=length(otu), ncol=length(samp), dimnames=list(otu, samp))
    t(ret)
}

abun_per_taxa <- function(comm) map_dbl(comm[,2:ncol(comm)], sum)

mean_per_taxa <- function(comm) map_dbl(comm[,2:ncol(comm)], mean)

abun_per_sample <- function(abun_df) {
    ret <- abun_df |>
        select(2:last_col()) |>
        rowwise() |>
        summarize(sum=sum(c_across())) |>
        pull(sum)

    names(ret) <- abun_df$index
    return(ret)
}

num_samples_with_otu <- function(abun_df, otu) {
    sum(pull(abun_df, otu) > 0)
}

prob_taxa_in_sample <- function(abun_df) {
    at <- abun_per_taxa(abun_df)
    as <- abun_per_sample(abun_df)
    total <- sum(as)
    expand_grid(
        enframe(at/total, "otu", "prop_otu"),
        enframe(as/total, "sample", "prop_sample")
    ) |>
        mutate(prob=1 - exp(-prop_otu*prop_sample))
}

dissemination <- function(abun_df, comm_df=NULL) {
    taxa <- colnames(select(comm_df, 2:last_col()))
    at <- abun_per_taxa(comm_df)
    st <- map_dbl(taxa, ~num_samples_with_otu(abun_df, .x))

    prob_taxa_in_sample(abun_df) |>
        group_by(otu) |>
        summarize(prob_occur=sum(prob)) |>
        filter(otu %in% taxa) |>
        mutate(diss_score=st / (prob_occur * at))
}

col_diversity <- function(comm, method="shannon") {
    norm <- function(col) {
        # lc <- ifelse(col > 0, log(col) + 1, 0)
        if (sum(col) == 0)
            return(0)
        p <- col / sum(col)
        p[p > 0]
    }
    if (method == "shannon")
        fun <- function(p) if (length(p) > 1) -sum(p * log(p)) else 0
    else if (method == "evenness")
        fun <- function(p) if (length(p) > 1) -sum(p * log(p)) / log(length(p)) else 0
    else
        fun <- function(p) 1 / sum(p^2)
    
    comm |>
        map(norm) |>
        map_dbl(fun)
}

beta_diversity <- function(comm, method="shannon") {
    comm |>
        select(2:last_col()) |>
        col_diversity(method)
}

alpha_diversity <- function(comm, method="shannon") {
    cmat <- select(comm, 2:last_col()) |> as.matrix()

    as_tibble(t(cmat), .name_repair="minimal") |>
        col_diversity(method)
    # comm |>
    #     rowwise() |>
    #     mutate(tot=sum(c_across(2:last_col()))) |>
    #     ungroup() |>
    #     transmute(across(2:last_col(1), ~.x / tot)) |>
    #     rowwise() |>
    #     summarize()
}

connectance <- function(comm) {
    comm_bin <- ifelse(comm[,2:ncol(comm)] > 0, 1, 0)
    links <- sum(comm_bin)
    links / (nrow(comm) * ncol(comm_bin))
}

sub_taxa_counts <- function(comm) {
    comm |>
        mutate(treat=str_extract(index, "[A-Z]+")) %>%
        split(.$treat) |>
        map_dbl(~sum(.x[2:(ncol(.x)-1)]))
}
