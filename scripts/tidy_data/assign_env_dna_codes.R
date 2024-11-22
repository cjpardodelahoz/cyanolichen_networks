#!/usr/bin/env Rscript

# Required packages and custom functions

# Tidy
library(tidyverse)

# Function to insert row to a dataframe
insert_row <- function(data, new_row, r_index) {
  data_new <- rbind(data[1:r_index-1, ],            
                    new_row,                
                    data[- (1:r_index-1), ])        
  data_new
}

# Function to insert row with negative controls
insert_negs <- function(data, neg_rows, neg_indices) {
  n_rows <- length(neg_rows)
  na_string <- rep(NA, ncol(data)-1)
  col_labels <- colnames(data)
  new_df <- data 
  for (i in 1:n_rows) {
    row_index <- neg_rows[i]
    neg_index <- neg_indices[i]
    neg_label <- paste("neg", neg_index, sep = "")
    neg_row <- c(na_string, neg_label)
    names(neg_row) <- col_labels
    new_df <- insert_row(new_df, neg_row, row_index)
  }
  new_df
}

# Assign DNA codes to environmental samples

# Load gap and nxtpelt metadata
gap_metadata <- read_csv("data/tables/gap_metadata.csv")
nxtpelt_metadata <- read_csv("data/tables/nxtpelt_metadata.csv") %>%
  mutate(moss = 
           as.character(moss))
# Make a vector with integers every 24 starting at 24 until 768
neg_rows_top <- c(seq(24, 768, 24), 781)
neg_indices_top <- seq(6, 38)
# Add unique code identifiers for each DNA sample
env_dna_codes_top <- gap_metadata %>%
  select(!c("gap_distance")) %>% 
  bind_rows(nxtpelt_metadata) %>%
  group_by(site, plot, quadrant) %>%
  arrange(.by_group = T) %>%
  ungroup() %>%
  filter(!is.na(moss)) %>%
  filter(layer == "top") %>%
  group_by(site) %>%
  mutate(site_n_samples = n()) %>%
  mutate(dna_code = 
           paste(site, "_", seq(1:site_n_samples[1]), "t", sep = "")) %>%
  select(colnames(nxtpelt_metadata), dna_code) %>%
  insert_negs(neg_rows = neg_rows_top, neg_indices = neg_indices_top) %>% # Add negative controls
  mutate(extraction_batch = ifelse(grepl("^neg", dna_code), dna_code, NA)) %>%
  ungroup() %>%
  fill(extraction_batch, .direction = "up")
env_dna_codes_bott <- gap_metadata %>%
  select(!c("gap_distance")) %>% 
  bind_rows(nxtpelt_metadata) %>%
  group_by(site, plot, quadrant) %>%
  arrange(.by_group = T) %>%
  ungroup() %>%
  filter(!is.na(moss)) %>%
  filter(layer == "bott") %>%
  group_by(site) %>%
  mutate(site_n_samples = n()) %>%
  mutate(dna_code = 
           paste(site, "_", seq(1:site_n_samples[1]), "b", sep = ""))


# Get table with sample and dna codes for env samples included in trial
env_dna_codes_all <- bind_rows(select(env_dna_codes_top, -extraction_batch), env_dna_codes_bott) %>%
  group_by(site, plot, quadrant) %>%
  arrange(.by_group = T)
trial_sites <- c("s8", "s10", "s12", "s13", "s15")
env_bench_codes <- paste(rep(c("c", "e", "m", "n"), each = 3), seq(1:3), "_top",
                         sep = "")
env_dna_codes_trial <- env_dna_codes_all %>%
  filter(site %in% trial_sites & plot == "high") %>%
  filter(ifelse(site == "s12", quadrant == "q1", quadrant == "q2")) %>%
  mutate(site_number = 
           parse_number(site)) %>%
  arrange(site_number) %>%
  select(colnames(nxtpelt_metadata), dna_code) %>%
  insert_negs(neg_rows = c(24, 48, 72, 96, 105), neg_indices = c(1, 2, 3, 4, 5)) %>%
  ungroup() %>%
  add_row(dna_code = env_bench_codes)
# Save trial env dna codes
dir.create("documents/tables", recursive = T)
write_csv(env_dna_codes_trial, "documents/tables/env_dna_codes_trial.csv")
write_csv(env_dna_codes_top, "documents/tables/env_dna_codes_top.csv")

# Barcode assignments for PacBio sequencing

# Generate table with barcode combinations
plate_col <- LETTERS[seq(1,8)]
plate_row <- seq(1,12)
plate_layout <- crossing(plate_col, plate_row) %>%
    arrange(plate_row) %>%
    mutate(plate_well = paste(plate_col, plate_row, sep = ""))
barcode_df <- bind_rows(plate_layout, plate_layout, plate_layout, plate_layout) %>%
    add_column(plate_id = rep(c("p1", "p2", "p3", "p4"), each = 96)) %>%
    arrange(plate_row) %>%
    add_column(rev_barcode_id = rep(paste("Kinnex16S_Rev_", seq(13,44,1), sep = ""), 12)) %>%
    arrange(plate_col, plate_id) %>%
    add_column(fwd_barcode_id = rep(paste("Kinnex16S_Fwd_", sprintf("%02d", seq(1:12)), sep = ""), 32)) %>%
    group_by(plate_id, plate_row) %>%
    arrange(desc(plate_col), .by_group = T) %>%
    ungroup()

# Vectors with dna codes from first batch of 384 envs samples
env_dna_codes_first_batch <- c(env_dna_codes_top$dna_code[1:96],      # env1
                                env_dna_codes_top$dna_code[97:144],   # This breakdown accounts for the env2 mistake I made
                                env_dna_codes_top$dna_code[169:192],
                                env_dna_codes_top$dna_code[145:168],
                                env_dna_codes_top$dna_code[193:288],  # env3
                                env_dna_codes_top$dna_code[289:384])  # env4
# Vector with dna codes from second batch of 384 envs samples
env_dna_codes_second_batch <- c(env_dna_codes_top$dna_code[385:480],  # env5
                                env_dna_codes_top$dna_code[481:576],   # env6
                                env_dna_codes_top$dna_code[577:672],  # env7
                                env_dna_codes_top$dna_code[673:768])  # env8

# Add dna codes to barcode tables
first_batch_barcoded <- barcode_df %>%
    add_column(dna_code = env_dna_codes_first_batch)
second_batch_barcoded <- barcode_df %>%
    add_column(dna_code = env_dna_codes_second_batch)
# Save barcode assignments for demultiplexing
first_batch_barcoded %>%
  mutate(Barcode = paste(fwd_barcode_id, rev_barcode_id, sep = "--"),
          Sample = dna_code) %>%
  select(Barcode, Sample) %>%
  write_csv("documents/tables/revio_order_10070_barcode.csv")
second_batch_barcoded %>%
  mutate(Barcode = paste(fwd_barcode_id, rev_barcode_id, sep = "--"),
          Sample = dna_code) %>%
  select(Barcode, Sample) %>%
  write_csv("documents/tables/revio_order_10413_barcode.csv")

# Join all metadata for top env samples

# Compile barcode data from all batches
all_batches_barcoded <- bind_rows(first_batch_barcoded, .id = "sequencing_batch") %>%
  select(-c(plate_col, plate_row))
# Join sequencing and env metadata
env_top_metadata <- env_dna_codes_top %>%
  left_join(all_batches_barcoded, by = "dna_code")
# Save metadata for top env samples
write_csv(env_top_metadata, "documents/tables/env_top_metadata.csv")
