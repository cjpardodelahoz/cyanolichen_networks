#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)

# Load table with ID ranges
pa_ranges <- read_csv("data/tables/pa_range_tidy.csv") %>%
  select(1:5)
# Load site voucher data
site_voucher <- read_csv("data/tables/site_voucher.csv")
# Expand ranges into specimen id table
lichen_dna_codes <- pa_ranges %>%
  rowwise() %>%
  mutate(dna_code = list(seq(first, last))) %>%
  unnest(dna_code) %>%
  select(-first, -last) %>%
  mutate(dna_code = paste("PA", dna_code, sep = "")) %>%
  left_join(site_voucher, by = "site")
# Save table with lichen DNA codes
write_csv(lichen_dna_codes, "documents/tables/lichen_dna_codes.csv")
