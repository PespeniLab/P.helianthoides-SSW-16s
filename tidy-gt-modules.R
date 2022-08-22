library(tidyverse)
library(reticulate)

mdict <- py_load_object("module-members.pickle", pickle="pickle")

get_block_types <- function(block) ifelse(str_detect(block, "k__"), "🦠", "⭐")

get_level_df <- function(lblocks) {
    imap_dfr(lblocks, ~tibble(block=.y, type=get_block_types(.x), label=.x))
}

tidy_level <- function(lblocks, lvl) tibble(level=lvl, get_level_df(lblocks))

imap_dfr(mdict, tidy_level) |> write_csv("hmembers-tidy.csv")
