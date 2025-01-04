#!/usr/bin/env Rscript

########## REQUIRED PACKAGES AND CUSTOM FUNCTIONS ##########

# Load required packages
library(tidyverse)
library(phyloseq)
library(ape)
library(treeio)
library(ggtree)

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

# Sort ASVs for trimming

# Get the lichenized ASVs
lichenized_asvs <- unique(lichenized_asv_table$asv_16s_first)
lichenized_asvs <- lichenized_asvs[lichenized_asvs != "contamination"]

# Get all the ASVs placed
all_asvs <- tbas_tree$tip.label[str_detect(tbas_tree$tip.label, "asv")] %>%
    str_remove("size.*")

# Non-lichenized ASVs
#misc_asvs <- c("asv4028", "asv1844") # from section 3.7 and 3.11, not identical to any genome
non_lichenized_asvs <- all_asvs[!all_asvs %in% lichenized_asvs]
#non_lichenized_asvs <- non_lichenized_asvs[!non_lichenized_asvs %in% misc_asvs]

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
lichenized_sections <- lichenized_asv_table$`Nostoc section`[]
lichenized_sections <- lichenized_sections[lichenized_sections != "contamination"] %>%
    unique()
#lichenized_sections <- c(lichenized_sections, "section_3.7", "section_3.11") # Present in the regional network but not identical to any genome

# Which asvs are in genome_asvs and non_lichenized_asvs
genome_asvs_in_non_lichenized <- genome_asvs[genome_asvs %in% non_lichenized_asvs]

# Get the taxa from the tbas metadata which are part of sections that are not in lichenized_sections
excluded_ref_taxa <- tbas_metadata %>%
    filter(!section %in% lichenized_sections | is.na(section)) %>%
    pull(taxon_name)

# Trim the tree

# Exclude non-lichenized ASVs, ASVs identical to a genome, and reference taxa that are not part of lichenized sections
tbas_tree_trimmed <- drop.tip(tbas_tree, c(non_lichenized_asvs, excluded_ref_taxa, genome_asvs))



ggtree(tbas_tree_trimmed) +
    geom_tiplab(aes(label = label), size = 2) +
    theme_tree2()
ggsave("tree.pdf", width = 10, height = 10)
