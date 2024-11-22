#!/usr/bin/env Rscript


########## REQUIRED PACKAGES AND CUSTOM FUNCTIONS ##########

# Load required packages
library(tidyverse)
library(RColorBrewer)
library(devtools)
library(SCRuB)
library(vegan)

# Function to run SCRuB in a single batch of samples
scrub_batch <- function(otu_table, metadata, batch_column, batch) {
    # Filter the metadata to only include samples from the current batch
    metadata_batch <- metadata %>%
        filter(!!sym(batch_column) == batch & dna_code %in% rownames(otu_table)) %>%
        mutate(is_control = case_when(
                    grepl("neg", dna_code) ~ TRUE,
                    .default = FALSE),
                sample_type = case_when(
                    grepl("neg", dna_code) ~ "extraction_control",
                    .default = "env")) %>%
        rename(sample_well = plate_well) %>%
        select(dna_code, is_control, sample_type, sample_well) %>%
        column_to_rownames("dna_code")
    # Filter the OTU table to only include samples from the current batch
    otu_table_batch <- otu_table[rownames(metadata_batch), ]
    # Run SCruB
    scrub_out <- SCRuB(otu_table_batch, metadata_batch)
    # Decontaminated OTU table with samples as columns
    decontaminated_samples <- scrub_out$decontaminated_samples %>% as.data.frame()
    # Return the decontaminated samples
    return(decontaminated_samples)
}
# Function to run scrub for all batches
scrub_wrapper <- function(otu_table, metadata, batch_type) {
    # If batch_type is "extraction", filter metadata using the extraction_batch column, if it is "sequencing", use the plate_id column
    if (batch_type == "extraction") {
        batch_column <- "extraction_batch"
    } else if (batch_type == "sequencing") {
        batch_column <- "plate_id"
    } else {
        stop("batch_type must be either 'extraction' or 'sequencing'")
    }
    # Get unique batch values
    batches <- metadata %>% filter(!is.na(plate_well)) %>% pull(batch_column) %>% unique()
    # Initialize a list to store the decontaminated samples
    decontaminated_samples_list <- list()
    # Map over the batches and run scrub_batch
    decontaminated_samples_list <- map(batches, ~scrub_batch(otu_table, metadata, batch_column, .x))
    # Bind the rows of all the decontaminated samples
    decontaminated_samples <- bind_rows(decontaminated_samples_list, otu_table[grepl("neg", rownames(otu_table)), ])
    # Return the decontaminated samples
    return(decontaminated_samples)
}
# Function to get NMDS (based on binary bray-curtis) from OTU table
plot_nmds <- function(otu_table, labels = TRUE) {
    # Remove empty rows and columns
    otu_table <- otu_table[rowSums(otu_table) > 0, colSums(otu_table) > 10]
    # Calculate Bray-Curtis dissimilarity
    bray_curtis <- vegdist(otu_table, method = "bray", binary = TRUE)
    # Perform NMDS analysis
    nmds <- metaMDS(bray_curtis)
    # Extract NMDS scores (site scores)
    nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
    # Add sample names to the NMDS scores
    nmds_scores$sample <- rownames(otu_table)
    # Calculate the log of the total number of reads (row sums)
    nmds_scores$total_reads_log <- log(rowSums(otu_table))
    # Add sample type information
    nmds_scores <- nmds_scores %>%
        mutate(sample_type = case_when(str_detect(sample, "neg") ~ "neg",
            .default = str_extract(sample, "^[^_]+"))
        )
    # Get custom color scale from the Paired palette
    n_colors <- nmds_scores$sample_type %>% unique() %>% length()
    palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_colors)
    # Plot NMDS
    plot  <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = sample_type, size = total_reads_log)) +
        geom_point() +
        scale_color_manual(values = palette) +
        labs(x = "NMDS1", y = "NMDS2", size = "Log(Total reads)", color = "DNA source") +
        theme(legend.position = "right") +
        theme(legend.text = element_text(size = 24),
            legend.title = element_text(size = 24),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 24)) +
        theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
            panel.grid.minor = element_line(color = "lightgray")) +
        guides(color = guide_legend(override.aes = list(size = 6)))
    # Add sample names to the plot
    if (labels) {
        plot <- plot + geom_text(aes(label = sample), hjust = 1, vjust = 1)
    }
    # Return the plot
    plot
}

########## REMOVE CONTAMINATION ##########

# Load data

# Read the environmental seuqencing layout metadata
metadata <- read_csv("documents/tables/env_top_metadata.csv")
# Read the OTU table read_delim with first column as row names
otu_table_raw <- read_delim(file = "analyses/environmental_sequencing/pacbio_env/unoise/counts/batch123_otutab.txt", delim = "\t") %>%
    column_to_rownames("#OTU ID") %>%
    t() %>% 
    as.data.frame()
# Remove samples with less than 7500 reads, unless they are negative controls
otu_table <- otu_table_raw[rowSums(otu_table_raw) > 7500 | grepl("neg", rownames(otu_table_raw)), ]
otu_table <- otu_table[!grepl("^PA", rownames(otu_table)), ]
#rownames(otu_table)

# Decontaminate

# Strategy 1: decontaminate with scrub grouping negative controls by sequencing batch
otu_table_scrub <- scrub_wrapper(otu_table, metadata, "sequencing")
# Strategy 2: remove all ASVs present in the negative controls
neg_controls <- otu_table[grep("^neg", rownames(otu_table)), ]
asvs_to_remove <- colnames(neg_controls)[colSums(neg_controls > 0) > 0]
otu_table_strict <- otu_table[, !colnames(otu_table) %in% asvs_to_remove]

########## VISUALIZE COMPOSITIONAL DIFFERENCES ##########

# Subset OTU tables to include only Nostocales ASVs

# Load the list of Nostocale ASVs
nostocales_asvs <- scan("analyses/environmental_sequencing/pacbio_env/unoise/taxonomy/kraken/batch123_asvs_nostocales.txt", what = "character") %>%
    str_remove(";.*")
# Subset the OTU tables to only include Nostocales ASVs
otu_table_scrub_nostocales <- otu_table_scrub %>% select(all_of(nostocales_asvs))
otu_table_strict_nostocales <- otu_table_strict %>% select(all_of(which(colnames(.) %in% nostocales_asvs)))

# Remove spurious counts and samples
otu_table_scrub <- otu_table_scrub[rowSums(otu_table_scrub) > 0, colSums(otu_table_scrub) > 10]
otu_table_strict <- otu_table_strict[rowSums(otu_table_strict) > 0, colSums(otu_table_strict) > 10]
otu_table_scrub_nostocales <- otu_table_scrub_nostocales[rowSums(otu_table_scrub_nostocales) > 500 | grepl("neg", rownames(otu_table_scrub_nostocales)), colSums(otu_table_scrub_nostocales) > 0]
otu_table_strict_nostocales <- otu_table_strict_nostocales[rowSums(otu_table_strict_nostocales) > 500, colSums(otu_table_strict_nostocales) > 0]

# NMDS for full community composition and Nostocales ASVs
plot_nmds(otu_table_raw, labels = TRUE)
ggsave("documents/plots/env_nmds_raw_plot.pdf")
plot_nmds(otu_table, labels = TRUE)
ggsave("documents/plots/env_nmds_filtered_plot1.pdf")
plot_nmds(otu_table_scrub, labels = TRUE)
ggsave("documents/plots/env_nmds_scrub_plot.pdf")
plot_nmds(otu_table_strict, labels = TRUE)
ggsave("documents/plots/env_nmds_strict_plot.pdf")
plot_nmds(otu_table_scrub_nostocales, labels = TRUE)
ggsave("documents/plots/env_nmds_scrub_nostocales_plot.pdf")
plot_nmds(otu_table_strict_nostocales, labels = TRUE)
ggsave("documents/plots/env_nmds_strict_nostocales_plot.pdf")