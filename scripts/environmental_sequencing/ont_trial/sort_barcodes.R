#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)

# Function to add sequencing barcode ids to a table with dna ids
# dna_id_df         A dafa frame with the DNA ids. 
# sample_indices    Integer vector with the sample indices to which the barcodes will be added
# barcode_indices   Integer vector with the barcode indices used. Order must match sample_indices
add_barcode_ids <- function(dna_id_df, sample_indices, barcode_indices) {
  # Generate table realting barcodes with sample indices
  plate_col <- LETTERS[seq(1,8)]
  plate_row <- seq(1,12)
  barcode_df <- crossing(plate_col, plate_row) %>%
    arrange(plate_row) %>%
    mutate(plate_well = paste(plate_col, plate_row, sep = "")) %>%
    add_column( barcode_id = paste("bc", seq(1:96), sep = "")) %>%
    group_by(plate_row) %>%
    arrange(desc(plate_col), .by_group = T) %>%
    add_column(barcode_index = seq(1:96)) %>%
    ungroup() %>%
    select(barcode_id, barcode_index)
  #
  index_key <- tibble(sample_index = sample_indices, barcode_index = barcode_indices)
  dna_id_df <- dna_id_df %>%
    add_column(sample_index = seq(1, nrow(.))) %>%
    left_join(index_key, by = "sample_index") %>%
    left_join(barcode_df, by = "barcode_index")
  dna_id_df
}


# Load table with dna codes for ont env trial samples
env_dna_codes_trial <- read_csv(file = "documents/tables/env_dna_codes_trial.csv") %>%
  add_barcode_ids(sample_indices = c(seq(1, 48), seq(97, 104)),
                  barcode_indices = c(seq(1, 48), seq(89, 96)))

