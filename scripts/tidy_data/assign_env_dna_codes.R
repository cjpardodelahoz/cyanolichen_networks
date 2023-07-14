#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
#
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

# Load gap and nxtpelt metadata
gap_metadata <- read_csv("data/tables/gap_metadata.csv")
nxtpelt_metadata <- read_csv("data/tables/nxtpelt_metadata.csv") %>%
  mutate(moss = 
           as.character(moss))
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
           paste(site, "_", seq(1:site_n_samples[1]), "t", sep = ""))
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
env_dna_codes_all <- bind_rows(env_dna_codes_top, env_dna_codes_bott) %>%
  group_by(site, plot, quadrant) %>%
  arrange(.by_group = T)
# Get table with sample and dna codes for env samples included in trial
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
  