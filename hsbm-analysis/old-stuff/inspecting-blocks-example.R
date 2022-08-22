library(tidyverse)

block_df <- read_csv("hmembers-tidy.csv")

# replace '23' with block of interest
block_members(block_df, 23)
block_members(block_df, 25)

# extract only the family, then count number of occurences
block_members(block_df, 96, tax_level="family") |>
    table()

