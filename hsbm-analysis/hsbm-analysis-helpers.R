library(tidyverse)

# Methods to convert Python results dictionary to a dataframe---------------------
tidy_level <- function(lblocks, lvl) tibble(level=as.integer(lvl), get_level_df(lblocks))

get_block_types <- function(block) ifelse(str_detect(block, "k__"), "🦠", "⭐")

get_level_df <- function(lblocks) {
    imap_dfr(lblocks, ~tibble(block=.y, type=get_block_types(.x), label=.x))
}

# Methods to extract the sub-abundance matrix for pairs of blocks-----------------
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

sub_taxa_counts <- function(comm) {
    comm |>
        mutate(treat=str_extract(index, "[A-Z]+")) %>%
        split(.$treat) |>
        map_dbl(~sum(.x[2:(ncol(.x)-1)]))
}

# Methods for summarizing and inspecting membership blocks------------------------
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

# Helpers for entropy analysis----------------------------------------------------
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
}
