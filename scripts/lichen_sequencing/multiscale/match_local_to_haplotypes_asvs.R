#!/usr/bin/env Rscript

########## REQUIRED PACKAGES AND CUSTOM FUNCTIONS ##########

# Load required packages
library(tidyverse)

########## MATCH LOCAL VOUCHER DATA WITH ITS, RBCLX HAPLOTYPES AND 16S ASVS ##########

# Site voucher
site_voucher <- read_csv("data/tables/site_voucher.csv")

# Apothecia data
apothecia_data <- read_csv("analyses/lichen_sequencing/ont_lichen/id/lichen_voucher_jm.csv") %>%
    select(dna_code, apothecia)

# Load local voucher V6
local_voucher_v6 <- read_csv("analyses/lichen_sequencing/ont_lichen/id/lichen_voucher_jm_curated_v6.csv")

# Load the table with IDs for the ASVs from lichenized Nostoc
lichenized_asv_table <- read_csv("analyses/lichen_sequencing/pacbio_lichen/selected_specimens_id_table_curated.csv") %>%
    select(c("asv_16s_first", "dna_code")) %>%
    rename("asv_16s" = "asv_16s_first")

# Load the haplotype clusters
its_haplotype_clusters <- read_delim("analyses/lichen_sequencing/multiscale/haplotypes/its/its_multiscale_haplotypes.tsv", delim = "\t", col_names = F) %>%
    filter(X1 != "C") %>%
    select(c(X2, X9)) %>%
    mutate(X2 = paste0("its_", X2))
rbclx_haplotype_clusters <- read_delim("analyses/lichen_sequencing/multiscale/haplotypes/rbclx/rbclx_multiscale_haplotypes.tsv", delim = "\t", col_names = F) %>%
    filter(X1 != "C") %>%
    select(c(X2, X9)) %>%
    mutate(X2 = paste0("rbclx_", X2))
colnames(its_haplotype_clusters) <- c("its_haplotype", "dna_code")
colnames(rbclx_haplotype_clusters) <- c("rbclx_haplotype", "dna_code")

# Join voucher with haplotype clusters and 16S ASVs for publication
local_voucher_v7 <- local_voucher_v6 %>%
    left_join(site_voucher, by = "site") %>%
    left_join(apothecia_data, by = "dna_code") %>%
    mutate(`Nostoc species complex` = ifelse(!is.na(`Nostoc species complex`) & !str_detect(`Nostoc species complex`, "low_depth|contaminant"), 
                                                    paste0("spp. complex ", `Nostoc species complex`), 
                                                    `Nostoc species complex`),
            `Nostoc section` = ifelse(!is.na(`Nostoc section`) & !str_detect(`Nostoc section`, "low_depth|contaminant"), 
                                                    paste0("section ", `Nostoc section`), 
                                                    `Nostoc section`),
            `Nostoc OTU` = ifelse(is.na(`Nostoc phylogroup`), 
                                coalesce(`Nostoc species complex`,`Nostoc phylogroup`),
                                `Nostoc phylogroup`)) %>%
    left_join(its_haplotype_clusters, by = "dna_code") %>%
    left_join(rbclx_haplotype_clusters, by = "dna_code") %>%
    mutate(its_haplotype = ifelse(str_detect(`Mycobiont molecular ID`, "low_depth|contaminant"), 
                                NA, 
                                its_haplotype), 
            rbclx_haplotype = ifelse(str_detect(`Nostoc OTU`, "low_depth|contaminant"),
                                NA,
                                rbclx_haplotype)) %>% # to prevent extrapolation to low depth samples
    left_join(lichenized_asv_table, by = "dna_code") %>%
    mutate(asv_16s_source = case_when(
        !is.na(asv_16s) ~ "sequenced",
        is.na(asv_16s) ~ "extrapolated",
    )) %>%
    group_by(rbclx_haplotype) %>%
    mutate(asv_16s = ifelse(is.na(asv_16s), 
                                  first(na.omit(asv_16s)), 
                                  asv_16s)) %>%
    ungroup() %>%
    mutate(asv_16s = case_when(`Nostoc phylogroup` == "low_depth" ~ NA,
        `Nostoc phylogroup` == "contaminant" ~ NA,
        .default = asv_16s)) # to prevent extrapolation to low depth samples

# Create version of local dataset with haplotype data for analyses only
local_data <- local_voucher_v6 %>%
    left_join(apothecia_data, by = "dna_code") %>%
    mutate(`Nostoc species complex` = ifelse(!is.na(`Nostoc species complex`) & !str_detect(`Nostoc species complex`, "low_depth|contaminant"), 
                                                    paste0("sppcomplex", str_remove(`Nostoc species complex`, "\\.")), 
                                                    `Nostoc species complex`),
            nostoc_otu = ifelse(is.na(`Nostoc phylogroup`), 
                                coalesce(`Nostoc species complex`,`Nostoc phylogroup`),
                                `Nostoc phylogroup`),
            `Mycobiont molecular ID` = str_remove_all(`Mycobiont molecular ID`, "_")) %>%
    left_join(its_haplotype_clusters, by = "dna_code") %>%
    left_join(rbclx_haplotype_clusters, by = "dna_code") %>%
    mutate(its_haplotype = ifelse(str_detect(`Mycobiont molecular ID`, "low_depth|contaminant"), 
                                NA, 
                                its_haplotype), 
            rbclx_haplotype = ifelse(str_detect(nostoc_otu, "low_depth|contaminant"),
                                NA,
                                rbclx_haplotype)) %>% # to prevent extrapolation to low depth samples
    left_join(lichenized_asv_table, by = "dna_code") %>%
    mutate(asv_16s_source = case_when(
        !is.na(asv_16s) ~ "sequenced",
        is.na(asv_16s) ~ "extrapolated",
    )) %>%
    group_by(rbclx_haplotype) %>%
    mutate(asv_16s = ifelse(is.na(asv_16s), 
                                  first(na.omit(asv_16s)), 
                                  asv_16s)) %>%
    ungroup() %>%
    mutate(asv_16s = case_when(`Nostoc phylogroup` == "low_depth" ~ NA,
        `Nostoc phylogroup` == "contaminant" ~ NA,
        .default = asv_16s)) # to prevent extrapolation to low depth samples

# Which dna_codes have no haplotype cluster?
no_its_haplotype <- local_voucher_v7 %>%
    filter(is.na(its_haplotype) & `Mycobiont molecular ID` != "low_depth") %>%
    select(dna_code)
no_rbclx_haplotype <- local_voucher_v7 %>%
    filter(is.na(rbclx_haplotype) & `Nostoc phylogroup` != "low_depth") %>%
    select(dna_code)

# Which dna_codes have no 16S ASV?
no_16s_asv <- local_voucher_v7 %>%
    filter(is.na(asv_16s) & `Nostoc phylogroup` != "low_depth") %>%
    select(dna_code)

# Write voucher table to file to curate IDs
write_csv(local_voucher_v7, "documents/tables/lichen_voucher_jm_curated_v7.csv")
write_csv(local_data, "data/tables/local_data.csv")

########## DETECTION TABLE OF ASVS IN PELTIGERA THALLI ##########

# Filter out NA and "contamination"
filtered_voucher <- local_data %>%
  filter(!is.na(asv_16s) & asv_16s != "contamination")

# Generate the binary scoring table
pelt_detection <- filtered_voucher %>%
  select(asv_16s, site) %>%
  distinct() %>%
  mutate(presence = "present") %>%
  pivot_wider(names_from = site, values_from = presence, values_fill = "absent")

# Save the detection table to RData
save(pelt_detection, file = "analyses/lichen_sequencing/multiscale/pelt_detection.RData")
