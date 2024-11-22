#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(spgs)
# Function to add DNA codes to barcode tables
barcode_dna_codes <- function(lichen_dna_codes, pap_set_values, barcode_type = c("ITS", "rbclx", "both"), its_barcode_df = its_barcode_df, rbclx_barcode_df = rbclx_barcode_df) {
  # Filter DNA codes based on pap_set values
  pool_dna_codes <- lichen_dna_codes %>%
    filter(pap_set %in% pap_set_values) %>%
    pull(dna_code)
  
  # Initialize result list
  result <- list()
  
  # Process ITS barcode table if requested
  if (barcode_type %in% c("ITS", "both") && !is.null(its_barcode_df)) {
    result$its_barcoded_df <- its_barcode_df[1:length(pool_dna_codes), ] %>%
      add_column(dna_code = pool_dna_codes) %>%
      left_join(primer_barcode_df, by = c("fwd_barcode_id" = "full_primer_id")) %>%
      rename(fwd_primer_barcode_plus = full_primer_sequence) %>%
      left_join(primer_barcode_df, by = c("rev_barcode_id" = "full_primer_id")) %>%
      rename(rev_primer_barcode_minus = full_primer_sequence) %>%
      mutate(rev_primer_barcode_plus = reverseComplement(rev_primer_barcode_minus) %>% str_to_upper(),
             fwd_primer_barcode_minus = reverseComplement(fwd_primer_barcode_plus) %>% str_to_upper())
  }
  
  # Process rbclx barcode table if requested
  if (barcode_type %in% c("rbclx", "both") && !is.null(rbclx_barcode_df)) {
    result$rbclx_barcoded_df <- rbclx_barcode_df[1:length(pool_dna_codes), ] %>%
      add_column(dna_code = pool_dna_codes) %>%
      left_join(primer_barcode_df, by = c("fwd_barcode_id" = "full_primer_id")) %>%
      rename(fwd_primer_barcode_plus = full_primer_sequence) %>%
      left_join(primer_barcode_df, by = c("rev_barcode_id" = "full_primer_id")) %>%
      rename(rev_primer_barcode_minus = full_primer_sequence) %>%
      mutate(rev_primer_barcode_plus = reverseComplement(rev_primer_barcode_minus) %>% str_to_upper(),
             fwd_primer_barcode_minus = reverseComplement(fwd_primer_barcode_plus) %>% str_to_upper())
  }
  
  return(result)
}

# Barcode assignments for ONT sequencing

# Generate table with barcode combinations
plate_col <- LETTERS[seq(1,8)]
plate_row <- seq(1,12)
plate_layout <- crossing(plate_col, plate_row) %>%
    arrange(plate_row) %>%
    mutate(plate_well = paste(plate_col, plate_row, sep = ""))
its_barcode_df <- bind_rows(plate_layout, plate_layout, plate_layout, plate_layout) %>%
    add_column(plate_id = rep(c("p1", "p2", "p3", "p4"), each = 96)) %>%
    arrange(plate_row) %>%
    add_column(rev_barcode_id = rep(paste("ITS4_R", seq(1:32), sep = ""), 12)) %>%
    arrange(plate_col, plate_id) %>%
    add_column(fwd_barcode_id = rep(paste("ITS1F_F", seq(1:12), sep = ""), 32)) %>%
    group_by(plate_id, plate_row) %>%
    arrange(desc(plate_col), .by_group = T) %>%
    ungroup()
rbclx_barcode_df <- bind_rows(plate_layout, plate_layout, plate_layout, plate_layout) %>%
    add_column(plate_id = rep(c("p1", "p2", "p3", "p4"), each = 96)) %>%
    arrange(plate_row) %>%
    add_column(rev_barcode_id = rep(paste("CX_R", seq(1:32), sep = ""), 12)) %>%
    arrange(plate_col, plate_id) %>%
    add_column(fwd_barcode_id = rep(paste("CW_F", seq(1:12), sep = ""), 32)) %>%
    group_by(plate_id, plate_row) %>%
    arrange(desc(plate_col), .by_group = T) %>%
    ungroup()

# Indices of empty wells
empty_wells <- c(
  "PA192", "PA193", "PA1052", "PA1455", "PA1937",
  paste("PA", seq(565, 566, 1), sep = ""),
  paste("PA", seq(2277, 2280, 1), sep = "")
)
# DNA codes from Jola arctic specimens (rbcLX)
jola1_dna_codes <- c(
  paste("PL", seq(647, 674, 1), sep = ""),
  "PL676",
  paste("PL", seq(682, 691, 1), sep = ""),
  paste("PL", seq(693, 701, 1), sep = ""),
  paste("PL", seq(703, 750, 1), sep = "")
)
# DNA codes from arctic specimens (ITS and rbcLX)
arctic_dna_codes <- c(paste("PL", seq(751, 796, 1), sep = ""))
# DNA codes from Diane's lichen specimens
diane_dna_codes <- c(paste("PA", seq(4101, 4117, 1), sep = ""))

# Sort DNA codes into the PAP sets
lichen_dna_codes <- read_csv("documents/tables/lichen_voucher_jm.csv") %>% # Note that this table does not include PA1937 even thoug it is included in the range for site 14
    select(dna_code) %>%
    add_row(dna_code = empty_wells) %>%
    add_row(dna_code = diane_dna_codes) %>%
    mutate(dna_index = str_remove(dna_code, "PA") %>% as.integer(),
           pap_set = case_when(
               dna_index <= 96 ~ "PAP1",
               dna_index <= 192 ~ "PAP2",
               dna_index <= 288 ~ "PAP3",
               dna_index <= 384 ~ "PAP4",
               dna_index <= 480 ~ "PAP5",
               dna_index <= 566 ~ "PAP6", # 564 is the last sample with DNA. Added 2 extra to account for empty wells
               dna_index <= 696 ~ "PAP7",
               dna_index <= 792 ~ "PAP8",
               dna_index <= 888 ~ "PAP9",
               dna_index <= 984 ~ "PAP10",
               dna_index <= 1080 ~ "PAP11",
               dna_index <= 1176 ~ "PAP12",
               dna_index <= 1272 ~ "PAP13",
               dna_index <= 1368 ~ "PAP14",
               dna_index <= 1464 ~ "PAP15",
               dna_index <= 1560 ~ "PAP16",
               dna_index <= 1656 ~ "PAP17",
               dna_index <= 1752 ~ "PAP18",
               dna_index <= 1848 ~ "PAP19",
               dna_index <= 1944 ~ "PAP20",
               dna_index <= 2040 ~ "PAP21",
               dna_index <= 2136 ~ "PAP22",
               dna_index <= 2232 ~ "PAP23",
               dna_index <= 3047 ~ "PAP24",
               dna_index <= 3143 ~ "PAP25",
               dna_index <= 4107 ~ "PAP26",
               dna_index <= 4117 ~ "PAP6",
           )) %>%
    add_row(dna_code = arctic_dna_codes, 
            pap_set = "arctic") %>%
    add_row(dna_code = jola1_dna_codes, pap_set = "jola1") %>%
    group_by(pap_set) %>%
    arrange(dna_index)

# Diretory to store demultiplex barcode keys
if (!dir.exists("analyses/lichen_sequencing/ont_lichen/demultiplex/")) {
  dir.create("analyses/lichen_sequencing/ont_lichen/demultiplex/", recursive = T)
}

# Define pool sets
pool1_sets <- c("PAP1", "PAP2", "PAP3")
pool2_sets <- c("PAP4", "PAP5", "PAP6", "arctic")
pool3_sets <- c("PAP7", "PAP8", "PAP9", "jola1")
pool4_sets <- c("PAP10", "PAP11", "PAP12")
pool5_sets <- c("PAP13", "PAP14", "PAP15")
pool6_sets <- c("PAP16", "PAP17", "PAP18", "PAP26")
pool7_sets <- c("PAP19", "PAP20", "PAP21")
pool8_sets <- c("PAP22", "PAP23", "PAP24", "PAP25")

# Load table with primer and barcode sequences and remove the padding sequence from the barcodes
primer_barcode_df <- read_csv("documents/tables/its_rbclx_barcoded_primers.csv") %>%
    mutate(full_primer_sequence = str_remove(full_primer_sequence, "GGTAG"))
    
# Assign sequence barcodes to the DNA codes and save the demultiplex keys (pools 1-8)
for (i in 1:8) {
  pool_name <- paste0("pool", i)
  pool_sets <- get(paste0("pool", i, "_sets"))

  # Add barcodes to DNA codes
  pool_barcoded <- barcode_dna_codes(lichen_dna_codes, pool_sets, barcode_type = "both", 
                                     its_barcode_df = its_barcode_df, rbclx_barcode_df = rbclx_barcode_df)
  
  # Define the file paths
  its_plus_path <- paste0("analyses/lichen_sequencing/ont_lichen/demultiplex/", pool_name, "_its_barcoded_plus.tsv")
  #its_minus_path <- paste0("analyses/lichen_sequencing/ont_lichen/demultiplex/", pool_name, "_its_barcoded_minus.tsv")
  rbclx_plus_path <- paste0("analyses/lichen_sequencing/ont_lichen/demultiplex/", pool_name, "_rbclx_barcoded_plus.tsv")
  #rbclx_minus_path <- paste0("analyses/lichen_sequencing/ont_lichen/demultiplex/", pool_name, "_rbclx_barcoded_minus.tsv")
  
  # Write ITS barcoded data
  pool_barcoded$its_barcoded_df %>% 
    select(dna_code, fwd_primer_barcode_plus, rev_primer_barcode_plus) %>%
    mutate(dna_code = paste(dna_code, "_its_plus", sep = "")) %>%
    write_delim(file = its_plus_path, delim = "\t", col_names = FALSE)
  
  #pool_barcoded$its_barcoded_df %>% 
    #select(dna_code, fwd_primer_barcode_minus, rev_primer_barcode_minus) %>%
    #mutate(dna_code = paste(dna_code, "_its_minus", sep = "")) %>%
    #write_delim(file = its_minus_path, delim = "\t", col_names = FALSE)
  
  # Write rbclx barcoded data
  pool_barcoded$rbclx_barcoded_df %>% 
    select(dna_code, fwd_primer_barcode_plus, rev_primer_barcode_plus) %>%
    mutate(dna_code = paste(dna_code, "_rbclx_plus", sep = "")) %>%
    write_delim(file = rbclx_plus_path, delim = "\t", col_names = FALSE)
  
  #pool_barcoded$rbclx_barcoded_df %>% 
    #select(dna_code, fwd_primer_barcode_minus, rev_primer_barcode_minus) %>%
    #mutate(dna_code = paste(dna_code, "_rbclx_minus", sep = "")) %>%
    #write_delim(file = rbclx_minus_path, delim = "\t", col_names = FALSE)
}
