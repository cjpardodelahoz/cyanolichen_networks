#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(spgs)

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

# Gerate a table with the barcode sequences in both orientations
barcoded_16s_forward <- read_csv(file = "documents/tables/barcoded_16s_forward.csv") %>%
  mutate(barcode_id = str_replace(Name, "bc0*", "bc")) %>%
  mutate(barcode_id = str_remove(barcode_id, "F")) %>%
  mutate(f_primer_f_sequence = str_to_upper(Sequence)) %>%
  mutate(f_barcode_f_sequence = str_remove(f_primer_f_sequence, "AGAGTTTGATCCTGGCTCAG")) %>%
  mutate(f_barcode_r_sequence = reverseComplement(f_barcode_f_sequence)) %>%
  mutate(f_barcode_r_sequence = str_to_upper(f_barcode_r_sequence)) %>%
  select(barcode_id, f_barcode_f_sequence, f_barcode_r_sequence, f_primer_f_sequence)
barcoded_16s_reverse <- read_csv(file = "documents/tables/barcoded_16s_reverse.csv") %>%
  mutate(barcode_id = str_replace(Name, "bc0*", "bc")) %>%
  mutate(barcode_id = str_remove(barcode_id, "R")) %>%
  mutate(r_primer_f_sequence = str_to_upper(Sequence)) %>%
  mutate(r_primer_r_sequence = reverseComplement(r_primer_f_sequence)) %>%
  mutate(r_primer_r_sequence = str_to_upper(r_primer_r_sequence)) %>%
  mutate(r_barcode_f_sequence = str_remove(r_primer_f_sequence, "GGTTACCTTGTTACGACTT")) %>%
  mutate(r_barcode_r_sequence = reverseComplement(r_barcode_f_sequence)) %>%
  mutate(r_barcode_r_sequence = str_to_upper(r_barcode_r_sequence)) %>%
  select(barcode_id, r_barcode_f_sequence, r_barcode_r_sequence, r_primer_r_sequence)
barcoded_16s_full <- barcoded_16s_forward %>%
  left_join(barcoded_16s_reverse, by = "barcode_id")
# Load table with dna codes for ont env trial samples and add barcode info
env_dna_codes_trial_barcoded <- read_csv(file = "documents/tables/env_dna_codes_trial.csv") %>%
  add_barcode_ids(sample_indices = c(seq(1, 48), seq(97, 104), seq(49, 96), seq(105, 117)),
                  barcode_indices = c(seq(1, 48), seq(89, 96), seq(49, 96), seq(1, 13))) %>%
  left_join(barcoded_16s_full, by = "barcode_id") %>%
  mutate(pool = case_when(
    between(sample_index, 1, 48) ~ "pool1",
    between(sample_index, 49, 96) ~ "pool2",
    between(sample_index, 97, 104) ~ "pool1",
    between(sample_index, 105, 117) ~ "pool2"
  ))
# Export barcode files for demultiplexing
if (!dir.exists("analyses/environmental_sequencing/ont_trial/demultiplex/")) {
  dir.create("analyses/environmental_sequencing/ont_trial/demultiplex/", recursive = T)
}
## Pool 1
env_dna_codes_trial_barcoded %>%
  filter(pool == "pool1") %>%
  select(dna_code, f_barcode_f_sequence) %>%
  mutate(dna_code = paste(dna_code, "_f", sep = "")) %>%
  write_delim(file = "analyses/environmental_sequencing/ont_trial/demultiplex/barcodes_pool1_f_f.tsv", 
              delim = "\t", col_names = F)
env_dna_codes_trial_barcoded %>%
  filter(pool == "pool1") %>%
  select(dna_code, r_barcode_r_sequence) %>%
  mutate(dna_code = paste(dna_code, "_r", sep = "")) %>%
  write_delim(file = "analyses/environmental_sequencing/ont_trial/demultiplex/barcodes_pool1_r_r.tsv", 
              delim = "\t", col_names = F)
## Pool 2
env_dna_codes_trial_barcoded %>%
  filter(pool == "pool2") %>%
  select(dna_code, f_barcode_f_sequence) %>%
  mutate(dna_code = paste(dna_code, "_f", sep = "")) %>%
  write_delim(file = "analyses/environmental_sequencing/ont_trial/demultiplex/barcodes_pool2_f_f.tsv", 
              delim = "\t", col_names = F)
env_dna_codes_trial_barcoded %>%
  filter(pool == "pool2") %>%
  select(dna_code, r_barcode_r_sequence) %>%
  mutate(dna_code = paste(dna_code, "_r", sep = "")) %>%
  write_delim(file = "analyses/environmental_sequencing/ont_trial/demultiplex/barcodes_pool2_r_r.tsv", 
              delim = "\t", col_names = F)
# Print files with list of samples in each pool
env_dna_codes_trial_barcoded %>%
  filter(pool == "pool1") %>%
  pull(dna_code) %>%
  write(file = "misc_files/environmental_sequencing/ont_trial/pool1_samples.txt")
env_dna_codes_trial_barcoded %>%
  filter(pool == "pool2") %>%
  pull(dna_code) %>%
  write(file = "misc_files/environmental_sequencing/ont_trial/pool2_samples.txt")
# Print file with sample ids and barcode-16s primer combinations for trimming
env_dna_codes_trial_barcoded %>%
  select(dna_code, f_primer_f_sequence, r_primer_r_sequence) %>%
  write_delim(file = "analyses/environmental_sequencing/ont_trial/demultiplex/barcoded_primers.tsv", 
              delim = "\t", col_names = F)
