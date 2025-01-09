#!/usr/bin/env Rscript

###### LOAD REQUIRED PACKAGES AND CUSTOM FUNCTIONS ######

# Load required packages and functions
library(tidyverse)

# Function to calculate haplotype sharing index (HSI) for a pseudoreplicate
calculate_hsi <- function(data) {
  n_records <- nrow(data)
  max_hap_freq <- max(table(data$haplotype_pair, useNA = "no"))
  hsi <- max_hap_freq / n_records
  return(data.frame(n_records = n_records, max_hap_freq = max_hap_freq, hsi = hsi))
}

###### PREPARE DATA FOR CLONALITY ANALYSIS ######

# Load regional ID data (from ABMIs sites)
abmi_data <- read_csv("data/tables/abmi_id_data.csv")

# Load curated lichen voucher
local_data <- read_csv("data/tables/local_data.csv")

# Load haplotype cluster data and add column names
its_haplotype_clusters <- read_delim("analyses/lichen_sequencing/multiscale/haplotypes/its/its_multiscale_haplotypes.tsv", delim = "\t", col_names = F) %>%
    filter(X1 != "C") %>%
    select(c(X2, X9)) %>%
    mutate(X2 = paste("its", X2, sep = "_")) %>%
    distinct(X2, X9)
rbclx_haplotype_clusters <- read_delim("analyses/lichen_sequencing/multiscale/haplotypes/rbclx/rbclx_multiscale_haplotypes.tsv", delim = "\t", col_names = F) %>%
    filter(X1 != "C") %>%
    select(c(X2, X9)) %>%
    mutate(X2 = paste("rbclx", X2, sep = "_")) %>%
    distinct(X2, X9)
colnames(its_haplotype_clusters) <- c("its_haplotype", "dna_code")
colnames(rbclx_haplotype_clusters) <- c("rbclx_haplotype", "dna_code")

# Add haplotype data to ABMI data
abmi_data <- abmi_data %>%
    left_join(its_haplotype_clusters, by = c("dna_id" = "dna_code")) %>%
    left_join(rbclx_haplotype_clusters, by = c("dna_id" = "dna_code")) %>%
    mutate(haplotype_pair = case_when(!is.na(its_haplotype) & !is.na(rbclx_haplotype) ~ paste(its_haplotype, rbclx_haplotype, sep = "_"),
                                      !is.na(its_haplotype) | !is.na(rbclx_haplotype) ~ NA))

# Pair haplotypes in local data
local_data <- local_data %>%
    mutate(haplotype_pair = case_when(!is.na(its_haplotype) & !is.na(rbclx_haplotype) ~ paste(its_haplotype, rbclx_haplotype, sep = "_"),
                                      !is.na(its_haplotype) | !is.na(rbclx_haplotype) ~ NA))

###### CALCULATE CLONALITY (HSI) ######

# Calculate clonality (HSI) at the local scale
local_clonality <- local_data %>%
    filter(!is.na(haplotype_pair) & str_detect(`Mycobiont molecular ID`, "Peltigera.*")) %>%
    group_by(`Mycobiont molecular ID`, site) %>%
    summarise(n_haplotype_pairs = n_distinct(haplotype_pair, na.rm = T),
              n_records = n(),
              max_hap_freq = max(table(haplotype_pair, useNA = "no"))) %>%
    ungroup() %>%
    mutate(hsi = max_hap_freq / n_records)

# Filter local clonality data to include only Peltigera species with at least 10 records per site
filtered_local_clonality <- local_clonality %>%
    filter(n_records >= 10)

# Filter the data to include only Peltigera species that occur in at least 10 sites
filtered_abmi_data <- abmi_data %>%
  filter(str_detect(`mycobiont_molecular_id`, "Peltigera.*") & !is.na(haplotype_pair)) %>%
  #filter(mycobiont_molecular_id == "Peltigeraaphthosa5") %>%
  group_by(`mycobiont_molecular_id`) %>%
  filter(n_distinct(site) >= 10) %>%
  ungroup()

# Generate pseudoreplicates and calculate HSI
set.seed(8260997)  # For reproducibility
regional_clonality <- filtered_abmi_data %>%
  group_by(`mycobiont_molecular_id`) %>%
  nest() %>%
  mutate(pseudo_replicates = map(data, ~ {
    sites <- unique(.x$site)
    replicate(1000, {
      sampled_data <- map_dfr(sites, function(site) {
        .x %>% filter(site == site) %>% sample_n(1)
      })
      calculate_hsi(sampled_data)
    }, simplify = FALSE)
  })) %>%
  unnest(pseudo_replicates) %>%
  unnest(pseudo_replicates) %>%
  mutate(pseudo_replicate_number = rep(1:1000, times = n() / 1000)) %>%
  select(`mycobiont_molecular_id`, pseudo_replicate_number, n_records, max_hap_freq, hsi)

# Save regional clonality as RData
save(regional_clonality, file = "analyses/ecology/regional_clonality.RData")

###### PREPARE CLONAILITY DATA FOR PLOTTING ######

# Load regional clonality data
load("analyses/ecology/regional_clonality.RData")

# Add module assignments to the data
load("analyses/ecology/peltigera_module_assignments.RData")
regional_clonality <- regional_clonality %>%
  left_join(peltigera_module_assignments, by = "mycobiont_molecular_id") %>%
    mutate(module = as.factor(module))

# Summarize regional clonality by calculating the median HSI and percentiles
summarized_regional_clonality <- regional_clonality %>%
    group_by(mycobiont_molecular_id, module) %>%
    summarise(median_hsi = median(hsi),
              lower_percentile = quantile(hsi, 0.0275),
              upper_percentile = quantile(hsi, 0.975)) %>%
    ungroup()

# Join the summarized regional clonality with the filtered local clonality
combined_clonality <- summarized_regional_clonality %>%
    left_join(filtered_local_clonality, by = c("mycobiont_molecular_id" = "Mycobiont molecular ID"))

# Reorder the x axis by module and median HSI
combined_clonality <- combined_clonality %>%
    mutate(module = factor(module, levels = c("1", "2", "4", "3", "6"))) %>%
    group_by(module) %>%
    mutate(mycobiont_molecular_id = fct_reorder(mycobiont_molecular_id, median_hsi, .desc = FALSE)) %>%
    ungroup()


###### PLOT REGIONAL AND LOCAL SCALE HSI ######

# Define custom colors for the modules
module_colors <- c("1" = "#8a708a", "2" = "#7e937b", "3" = "#3277b0", "4" = "#be8551", "5" = "gray70", "6" = "gray20")

# Define a custom theme for plots
custom_theme <- theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), 
    axis.text = element_text(color = "black", size = 7)
)

# Update the plot with the reordered x axis and error bars
clonality_pairs <- ggplot(combined_clonality, aes(x = mycobiont_molecular_id, y = median_hsi, color = module)) +
    geom_pointrange(aes(ymin = lower_percentile, ymax = upper_percentile)) +
    geom_point(aes(y = hsi), position = position_jitter(width = 0.05), size = 2, shape = 8) +
    scale_color_manual(values = module_colors) +
    labs(x = "Peltigera species", y = "Max. haplotype proportion") +
    geom_errorbar(aes(ymin = lower_percentile, ymax = upper_percentile), width = 0.25) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
    #theme_minimal() +
    custom_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the plot as a PDF
ggsave("documents/plots/clonality_pairs.pdf", units = "cm", width = 20, height = 10)
