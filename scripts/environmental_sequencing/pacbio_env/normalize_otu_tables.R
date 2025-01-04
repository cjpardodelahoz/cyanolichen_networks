#!/usr/bin/env Rscript

########## REQUIRED PACKAGES AND CUSTOM FUNCTIONS ##########

# Load required packages
library(tidyverse)
library(phyloseq)
library(ranacapa)

# Custom function to pool samples
pool_samples <- function(physeq, group_vars, original_taxa_names) {
  # Pool samples by specified variables
  physeq_pooled <- merge_samples(physeq, group = interaction(sample_data(physeq)[, group_vars], drop = TRUE))
  
  # Update sample data to reflect the new pooled samples
  pooled_sample_data <- sample_data(physeq_pooled)
  for (var in group_vars) {
    pooled_sample_data[[var]] <- sapply(strsplit(sample_names(physeq_pooled), "\\."), `[`, which(group_vars == var))
  }
  sample_data(physeq_pooled) <- pooled_sample_data
  
  # Set original taxa names
  taxa_names(physeq_pooled) <- original_taxa_names
  
  return(physeq_pooled)
}

# Custom function to create rarefaction plots
create_rarefaction_plot <- function(physeq, plot_title) {
  # Create rarefaction curves faceted by site
  rarefaction_plot <- ggrare(physeq, step = 1000, color = NULL, label = "Sample", se = FALSE) +
    facet_wrap(~ site) +
    theme_minimal() +
    scale_color_manual(values = "black") +
    labs(
      title = plot_title,
      x = "Number of Sequences",
      y = "Number of OTUs"
    )
  
  # Print the rarefaction plot
  print(rarefaction_plot)
}

# Custom function to normalize OTU table with total sum scaling
normalize_otu_table <- function(physeq, original_taxa_names) {
  otu_table_normalized <- transform_sample_counts(physeq, function(x) x / sum(x))
  taxa_names(otu_table_normalized) <- original_taxa_names
  return(otu_table_normalized)
}

# Custom function to convert OTU table to binary with a read threshold
binary_otu_table <- function(physeq, threshold = 10, original_taxa_names) {
  otu_table_binary <- transform_sample_counts(physeq, function(x) as.numeric(x >= threshold))
  taxa_names(otu_table_binary) <- original_taxa_names
  return(otu_table_binary)
}

# Function to process OTU tables
process_otu_table <- function(otu_table, label) {
  # Filter the metadata to only include samples present in the OTU table
  sample_ids <- rownames(otu_table)
  metadata_filtered <- metadata %>%
    filter(dna_code %in% sample_ids)
  
  # Remove samples with plate_id p4_3
  metadata_filtered <- metadata_filtered %>%
    filter(plate_id != "p4_3")
  
  # Create the phyloseq object
  otu_table <- otu_table(otu_table[metadata_filtered$dna_code, ], taxa_are_rows = FALSE)
  sample_data <- sample_data(metadata_filtered)
  sample_names(sample_data) <- metadata_filtered$dna_code
  physeq <- phyloseq(otu_table, sample_data)
  
  # Set original taxa names
  original_taxa_names <- colnames(otu_table)
  taxa_names(physeq) <- original_taxa_names
  
  # Print the phyloseq object to verify
  print(physeq)
  
  # Assign the original phyloseq object to the environment
  assign(paste0("physeq_", label), physeq, envir = .GlobalEnv)
  
  # Create normalized version of the non-pooled phyloseq object
  physeq_normalized <- normalize_otu_table(physeq, original_taxa_names)
  assign(paste0("physeq_", label, "_normalized"), physeq_normalized, envir = .GlobalEnv)
  
  # Pool samples and create rarefaction plots for different pooling versions
  physeq_quadrant_type <- pool_samples(physeq, c("site", "plot", "quadrant", "sample_type"), original_taxa_names)
  assign(paste0("physeq_", label, "_quadrant_type"), physeq_quadrant_type, envir = .GlobalEnv)
  assign(paste0("rarefaction_", label, "_quadrant_type"), create_rarefaction_plot(physeq_quadrant_type, paste0("Rarefaction Curves by Site, Plot, Quadrant, Sample Type (", label, ")")), envir = .GlobalEnv)
  
  physeq_plot_type <- pool_samples(physeq, c("site", "plot", "sample_type"), original_taxa_names)
  assign(paste0("physeq_", label, "_plot_type"), physeq_plot_type, envir = .GlobalEnv)
  assign(paste0("rarefaction_", label, "_plot_type"), create_rarefaction_plot(physeq_plot_type, paste0("Rarefaction Curves by Site, Plot, Sample Type (", label, ")")), envir = .GlobalEnv)
  
  physeq_site_type <- pool_samples(physeq, c("site", "sample_type"), original_taxa_names)
  assign(paste0("physeq_", label, "_site_type"), physeq_site_type, envir = .GlobalEnv)
  assign(paste0("rarefaction_", label, "_site_type"), create_rarefaction_plot(physeq_site_type, paste0("Rarefaction Curves by Site, Sample Type (", label, ")")), envir = .GlobalEnv)
  
  physeq_site <- pool_samples(physeq, c("site"), original_taxa_names)
  assign(paste0("physeq_", label, "_site"), physeq_site, envir = .GlobalEnv)
  assign(paste0("rarefaction_", label, "_site"), create_rarefaction_plot(physeq_site, paste0("Rarefaction Curves by Site (", label, ")")), envir = .GlobalEnv)
  
  # Normalize OTU table with total sum scaling for pooled versions
  assign(paste0("physeq_", label, "_quadrant_type_normalized"), normalize_otu_table(physeq_quadrant_type, original_taxa_names), envir = .GlobalEnv)
  assign(paste0("physeq_", label, "_plot_type_normalized"), normalize_otu_table(physeq_plot_type, original_taxa_names), envir = .GlobalEnv)
  assign(paste0("physeq_", label, "_site_type_normalized"), normalize_otu_table(physeq_site_type, original_taxa_names), envir = .GlobalEnv)
  assign(paste0("physeq_", label, "_site_normalized"), normalize_otu_table(physeq_site, original_taxa_names), envir = .GlobalEnv)
  
  # Convert OTU table to binary with a read threshold
  physeq_binary <- binary_otu_table(physeq, 10, original_taxa_names)
  assign(paste0("physeq_", label, "_binary"), physeq_binary, envir = .GlobalEnv)
  
  # Pool the binary phyloseq object for different pooling versions
  assign(paste0("physeq_", label, "_binary_quadrant_type"), pool_samples(physeq_binary, c("site", "plot", "quadrant", "sample_type"), original_taxa_names), envir = .GlobalEnv)
  assign(paste0("physeq_", label, "_binary_plot_type"), pool_samples(physeq_binary, c("site", "plot", "sample_type"), original_taxa_names), envir = .GlobalEnv)
  assign(paste0("physeq_", label, "_binary_site_type"), pool_samples(physeq_binary, c("site", "sample_type"), original_taxa_names), envir = .GlobalEnv)
  assign(paste0("physeq_", label, "_binary_site"), pool_samples(physeq_binary, c("site"), original_taxa_names), envir = .GlobalEnv)
}

# Function to check normalization
check_normalization <- function(physeq_normalized) {
  otu_table_normalized <- otu_table(physeq_normalized)
  sample_sums <- rowSums(otu_table_normalized)
  print("Sample sums (should be 1 for each sample):")
  print(sample_sums)
}

# Function to check binary conversion
check_binary_conversion <- function(physeq_binary) {
  otu_table_binary <- otu_table(physeq_binary)
  print("OTU table (binary format, should contain only 0s and 1s):")
  print(otu_table_binary)
}

########## POOLING AND PRAREFACTION CURVES ##########

# Load the OTU tables
load("analyses/environmental_sequencing/pacbio_env/unoise/counts/otu_tables.RData")

# Load the environmental sequencing metadata
metadata <- read_csv("documents/tables/env_top_metadata.csv")

# Process otu_table_strict
process_otu_table(otu_table_strict, "strict")

# Process otu_table_scrub
process_otu_table(otu_table_scrub, "scrub")

# Verify normalization, and binary conversion for strict
print("Checking normalization for strict:")
check_normalization(physeq_strict_normalized)

print("Checking binary conversion for strict:")
check_binary_conversion(physeq_strict_binary)

# Verify normalization, and binary conversion for scrub
print("Checking normalization for scrub:")
check_normalization(physeq_scrub_normalized)

print("Checking binary conversion for scrub:")
check_binary_conversion(physeq_scrub_binary)

# Save all the rarefaction plots to "documents/plots/" as pdfs
ggsave("documents/plots/rarefaction_strict_quadrant_type.pdf", plot = rarefaction_strict_quadrant_type, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_strict_plot_type.pdf", plot = rarefaction_strict_plot_type, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_strict_site_type.pdf", plot = rarefaction_strict_site_type, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_strict_site.pdf", plot = rarefaction_strict_site, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_scrub_quadrant_type.pdf", plot = rarefaction_scrub_quadrant_type, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_scrub_plot_type.pdf", plot = rarefaction_scrub_plot_type, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_scrub_site_type.pdf", plot = rarefaction_scrub_site_type, device = "pdf", width = 30, height = 20, units = "cm")
ggsave("documents/plots/rarefaction_scrub_site.pdf", plot = rarefaction_scrub_site, device = "pdf", width = 30, height = 20, units = "cm")

# Save the processed phyloseq objects to "analyses/environmental_sequencing/pacbio_env/unoise/counts/phyloseq_objects.RData"
save(
  physeq_strict, physeq_strict_normalized, physeq_strict_quadrant_type, physeq_strict_plot_type, physeq_strict_site_type, physeq_strict_site,
  physeq_strict_quadrant_type_normalized, physeq_strict_plot_type_normalized, physeq_strict_site_type_normalized, physeq_strict_site_normalized,
  physeq_strict_binary, physeq_strict_binary_quadrant_type, physeq_strict_binary_plot_type, physeq_strict_binary_site_type, physeq_strict_binary_site,
  file = "analyses/environmental_sequencing/pacbio_env/unoise/counts/phyloseq_objects_strict.RData"
)
save(
  physeq_scrub, physeq_scrub_normalized, physeq_scrub_quadrant_type, physeq_scrub_plot_type, physeq_scrub_site_type, physeq_scrub_site,
  physeq_scrub_quadrant_type_normalized, physeq_scrub_plot_type_normalized, physeq_scrub_site_type_normalized, physeq_scrub_site_normalized,
  physeq_scrub_binary, physeq_scrub_binary_quadrant_type, physeq_scrub_binary_plot_type, physeq_scrub_binary_site_type, physeq_scrub_binary_site,
  file = "analyses/environmental_sequencing/pacbio_env/unoise/counts/phyloseq_objects_scrub.RData"
)