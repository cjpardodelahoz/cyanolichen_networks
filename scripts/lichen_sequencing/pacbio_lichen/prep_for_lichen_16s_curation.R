#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)

# Load list of specimens selected for 16S sequencing
selected_specimens <- scan("analyses/lichen_sequencing/ont_lichen/consensus/unique_nostoc_ids.txt", what = "character") 

# Load curated lichen voucher and subset to the selected specimens
selected_specimens_id_table <- read_csv("documents/tables/lichen_voucher_jm_curated_v4.csv") %>%
   select(dna_code, `Mycobiont molecular ID`, `Nostoc section`, `Nostoc species complex`, `Nostoc phylogroup`) %>%
   filter(dna_code %in% selected_specimens)

# Write selected specimen ID table
write_csv(selected_specimens_id_table, "analyses/lichen_sequencing/pacbio_lichen/selected_specimens_id_table.csv")

# Read the OTU table
otu_table_raw <- read_delim(file = "analyses/environmental_sequencing/pacbio_env/unoise/counts/batch123_otutab.txt", delim = "\t") %>%
   column_to_rownames("#OTU ID") %>%
   t() %>% 
   as.data.frame()

# Load the list of Nostocale ASVs
nostocales_asvs <- scan("analyses/environmental_sequencing/pacbio_env/unoise/taxonomy/kraken/batch123_asvs_nostocales.txt", what = "character") %>%
   str_remove(";.*")

# Subset OTU table to PA samples (i.e., Peltigera specimens)
pa_otu_table <- otu_table_raw[grepl("^PA", rownames(otu_table_raw)), ]

# Subset the OTU table to include only ASVs identified as Nostocales by Kraken
pa_otu_table_nostocales <- pa_otu_table %>% select(all_of(which(colnames(.) %in% nostocales_asvs)))

# Remove rows and columns with 0 counts
pa_otu_table_nostocales <- pa_otu_table_nostocales[rowSums(pa_otu_table_nostocales) > 0, colSums(pa_otu_table_nostocales) > 10]

# Get the top 10 ASVs for each lichen sample
top10_nostocales_asvs_per_lichen <- pa_otu_table_nostocales %>%
   rownames_to_column("sample") %>%
   gather(asv, count, -sample) %>%
   arrange(desc(count)) %>%
   group_by(sample) %>%
   slice(1:10) %>%
   ungroup()

# Print the top 10 ASVs for each lichen sample
write_csv(top10_nostocales_asvs_per_lichen, "analyses/lichen_sequencing/pacbio_lichen/top10_nostocales_asvs_per_lichen.csv")
