#!/usr/bin/env Rscript

########## REQUIRED PACKAGES AND CUSTOM FUNCTIONS ##########

# Load required packages
library(tidyverse)
library(phyloseq)
library(ape)
library(treeio)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(aplot)

########## SUBSET PLACEMENT TREE ##########

# Load data

# Load tree with TBAS placement of Nostocales ASVs
tbas_tree <- read.tree("analyses/environmental_sequencing/pacbio_env/unoise/taxonomy/tbas/tbas_nostocales_asvs.tree")
tbas_tree$tip.label <- str_remove(tbas_tree$tip.label, "size.*")

# Load metadata for TBAS reference tree and join genome ID key
genome_id_key <- read_csv("data/tables/genome_id_key.csv")
tbas_metadata <- read_csv("data/tables/tbas_metadata.csv") %>%
    left_join(genome_id_key, by = c("taxon_name" = "Taxon name"))

# Load the table with IDs for the ASVs from lichenized Nostoc
lichenized_asv_table <- read_csv("analyses/lichen_sequencing/pacbio_lichen/selected_specimens_id_table_curated.csv")

# Load the results of the blastn search of Nostocales ASVs against the reference Nostoc database
blast_to_ref_nostoc <- read_delim("analyses/environmental_sequencing/pacbio_env/unoise/blast/out/batch123_asvs_nostocales_to_nostoc_ref.txt",
                                col_names = F) %>%
    mutate(X1 = str_remove(X1, ";.*")) %>%
    left_join(genome_id_key, by = c("X2" = "genome_id")) %>%
    mutate(X2 = `Taxon name`) %>%
    select(-c(`Taxon name`))
colnames(blast_to_ref_nostoc) <- c("query", "subject", "alignment_length", "p_identity", "n_identity", "mismatches", "gaps")

# Load the phyloseq objects and remove the unnecessary pools
load("analyses/environmental_sequencing/pacbio_env/unoise/counts/phyloseq_objects_strict.RData")
load("analyses/environmental_sequencing/pacbio_env/unoise/counts/phyloseq_objects_scrub.RData")
rm(physeq_strict, physeq_strict_normalized, physeq_strict_quadrant_type, physeq_strict_plot_type, physeq_strict_site_type, physeq_strict_site,
  physeq_strict_quadrant_type_normalized, physeq_strict_plot_type_normalized, physeq_strict_site_type_normalized, physeq_strict_site_normalized,
  physeq_strict_binary, physeq_strict_binary_quadrant_type, physeq_strict_binary_plot_type, physeq_strict_binary_site_type,
  physeq_scrub, physeq_scrub_normalized, physeq_scrub_quadrant_type, physeq_scrub_plot_type, physeq_scrub_site_type, physeq_scrub_site,
  physeq_scrub_quadrant_type_normalized, physeq_scrub_plot_type_normalized, physeq_scrub_site_type_normalized, physeq_scrub_site_normalized,
  physeq_scrub_binary, physeq_scrub_binary_quadrant_type, physeq_scrub_binary_plot_type, physeq_scrub_binary_site_type
  )

# Load module assignments
load("analyses/ecology/peltigera_module_assignments.RData")

# Sort ASVs for trimming

# Get the lichenized ASVs
lichenized_asvs <- unique(lichenized_asv_table$asv_16s_first)
lichenized_asvs <- lichenized_asvs[lichenized_asvs != "contamination" & lichenized_asvs != "asv17464"]

# Get all the ASVs placed
all_asvs <- tbas_tree$tip.label[str_detect(tbas_tree$tip.label, "asv")] %>%
    str_remove("size.*")

# Non-lichenized ASVs
misc_asvs <- c("asv4028", "asv1844") # from section 3.7 and 3.11, not identical to any genome
non_lichenized_asvs <- all_asvs[!all_asvs %in% lichenized_asvs]
non_lichenized_asvs <- non_lichenized_asvs[!non_lichenized_asvs %in% misc_asvs]

# Keys for ASVs that are identical to a sequence from the Nostoc reference
excluded_asvs <- c("asv9394", "asv16729") # From section 2.1, which is not included in the regional network
genome_asv_key <- blast_to_ref_nostoc %>%
    filter(p_identity == 100) %>%
    select(query, subject) %>% 
    filter(!query %in% excluded_asvs)
genome_asv_key_unique <- genome_asv_key %>%
    distinct(query, .keep_all = T)

# Get ASVs that are identical to a sequence from the Nostoc reference
genome_asvs <- genome_asv_key_unique$query

# Sort reference taxa for trimming

# Get the sections that include lichenized Nostoc
#lichenized_sections <- lichenized_asv_table$`Nostoc section`[]
#lichenized_sections <- lichenized_sections[lichenized_sections != "contamination"] %>%
#    unique()
#lichenized_sections <- c(lichenized_sections, "section_3.7", "section_3.11") # Present in the regional network but not identical to any genome

# Which asvs are in genome_asvs and non_lichenized_asvs
#genome_asvs_in_non_lichenized <- genome_asvs[genome_asvs %in% non_lichenized_asvs]

# Get the taxa from the tbas metadata which are part of sections that are not in lichenized_sections
#excluded_ref_taxa <- tbas_metadata %>%
#    filter(!section %in% lichenized_sections | is.na(section)) %>%
#    pull(taxon_name)

# Get the genomes that match an ASV
#genomes_matching_asvs <- genome_asv_key_unique$subject %>%
#    unique()

# Prep module data for plotting

# Get module data for the ASVs
module_key <- lichenized_asv_table %>%
    left_join(nostoc_module_assignments, by = "nostoc_otu") %>%
    filter(!is.na(module)) %>%
    mutate(module = factor(as.character(module))) %>%
    select(asv_16s_first, module) %>%
    distinct(asv_16s_first, .keep_all = T)

# Prep OTU table for trimming an plotting 

# Keep only lichenized ASVs in the OTU tables
physeq_strict_binary_site <- prune_taxa(lichenized_asvs, physeq_strict_binary_site)
physeq_scrub_binary_site <- prune_taxa(lichenized_asvs, physeq_scrub_binary_site)

# asvs are in lichenized_asvs but not in the colnames of the OTU tables
asvs_missing_strict <- lichenized_asvs[!lichenized_asvs %in% colnames(otu_table(physeq_strict_binary_site))]
asvs_missing_scrub <- lichenized_asvs[!lichenized_asvs %in% colnames(otu_table(physeq_scrub_binary_site))]

# Convert the scrub OTU table to dataframe and transpose it
otu_table_scrub <- otu_table(physeq_scrub_binary_site)
env_detection <- t(otu_table_scrub) %>%
    apply(., 2, function(x) ifelse(x > 0, "present", "absent")) %>%
    as.data.frame() %>%
    rownames_to_column(var = "tip.label") %>%
    left_join(module_key, by = c("tip.label" = "asv_16s_first")) #%>%
    # In the asv column, replace the asvs that are in the genome_asv_key with the genome name (subject column)
    #mutate(tip.label = ifelse(tip.label %in% genome_asv_key_unique$query,
    #    genome_asv_key_unique$subject[match(tip.label, genome_asv_key_unique$query)],
    #    tip.label))

# Trim the tree

# Exclude non-lichenized ASVs, ASVs identical to a genome, and reference taxa that are not part of lichenized sections
#tbas_tree_trimmed <- drop.tip(tbas_tree, c(non_lichenized_asvs, excluded_ref_taxa, genome_asvs))
tbas_tree_trimmed <- keep.tip(tbas_tree, c(lichenized_asvs[!lichenized_asvs %in% asvs_missing_scrub]))

# Trimmed tree with full reference for inset
tbas_tree_inset <- drop.tip(tbas_tree, c(non_lichenized_asvs, genome_asvs))

# Plot and save the inset tree plot
ggtree(tbas_tree_inset)
ggsave("documents/plots/tree_inset.pdf", unit = "cm", width = 10, height = 20)

########## GVP AND SCY GENE DATA FOR HEATMAPS ##########

# Mock gvp and scy data
gvp_data <- env_detection %>%
    select(tip.label) %>%
    mutate(gvpA = sample(0:3, nrow(.), replace = T),
           gvpC = sample(0:3, nrow(.), replace = T))
scy_data <- env_detection %>%
    select(tip.label) %>%
    mutate(scyA = sample(0:3, nrow(.), replace = T),
           scyC = sample(0:3, nrow(.), replace = T))

# Reshape gvp and scy data to long format for plotting
gvp_data_long <- gvp_data %>%
    pivot_longer(cols = -tip.label, names_to = "gene", values_to = "copies")
scy_data_long <- scy_data %>%
    pivot_longer(cols = -tip.label, names_to = "gene", values_to = "copies")

########## PLOT TREE WITH ENV DETECTION ##########

# Define custom colors for the modules
module_colors <- c("1" = "#8a708a", "2" = "#7e937b", "3" = "#3277b0", "4" = "#be8551", "5" = "gray70", "6" = "gray20")

# Define the order of the sites
sites <- c("s5", "s4", "s8", "s3", "s6", "s13", "s14", "s2", "s15", "s7", "s1", "s12", "s11", "s10", "s9")

# Plot the tbas trimmed tree with the environmental detection
# Sites are ordered by NR and elevation, except in the grassland
# RM, FH, BO, PK, and GR
# s5, s4, s8, s3, s6, s13, s14, s2, s15, s7, s1, s12, s11, s10, s9
base_tree_plot <- ggtree(tbas_tree_trimmed) %<+% env_detection
tree_env_detection_plot <- base_tree_plot #+
    #geom_tiplab(aes(label = label), align = T, offset = 0.03)
for (site in sites) {
    tree_env_detection_plot <- tree_env_detection_plot +
        new_scale(new_aes = "shape") +
        new_scale_color() +
        geom_fruit(geom = geom_point, 
                             mapping = aes_string(y = "label", shape = site, color = "module"),
                             offset = 0.03, size = 3) +
        scale_shape_manual(values = c(46, 19)) +
        scale_color_manual(values = module_colors)
}
tree_env_detection_plot <- tree_env_detection_plot + theme(legend.position = "none")

# Add a heatmap of the GVP and SCY genes to the right of the tree
gene_copy_colors <- c("0" = "#fbf3e5", "1" = "#ead096", "2" = "#ad9a70", "3" = "#74694d")
tree_env_detection_gvp_plot <- tree_env_detection_plot +
    geom_fruit(data = gvp_data_long, geom = geom_tile, 
        mapping = aes(x = gene, y = tip.label, fill = as.character(copies)), width = 2.1, offset = 0.08) +
    geom_fruit(data = scy_data_long, geom = geom_tile, 
        mapping = aes(x = gene, y = tip.label, fill = as.character(copies)), width = 2.1, offset = 0.11) +
    scale_fill_manual(values = gene_copy_colors)

ggsave("tree.pdf", width = 10, height = 10)
