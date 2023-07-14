#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)

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
env_dna_codes_trial <- env_dna_codes_all %>%
  filter(site %in% trial_sites & plot == "high") %>%
  filter(ifelse(site == "s12", quadrant == "q1", quadrant == "q2")) %>%
  mutate(site_number = 
           parse_number(site)) %>%
  arrange(site_number) %>%
  select(colnames(nxtpelt_metadata), dna_code)

# Save trial env dna codes
dir.create("documents/tables", recursive = T)
write_csv(env_dna_codes_trial, "documents/tables/env_dna_codes_trial.csv")
write_csv(env_dna_codes_top, "documents/tables/env_dna_codes_top.csv")
  