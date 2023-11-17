#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(treeio)
library(ape)

# Parse placement tree

# Load trees with EPA placements
nostoc_corrected_tree <- read.tree("analyses/environmental_sequencing/ont_trial/placement/nostoc/trees/RAxML_labelledTree.corrected_nostoc_epa_result")
# Save placement tree in newick to open with FigTree
write.tree(nostoc_corrected_tree, "analyses/environmental_sequencing/ont_trial/placement/nostoc/trees/nostoc_corrected_placement.tree")

# BLAST results

# Load Nostoc ref clade assignments
nostoc_ref_clade_assignments <- read_csv("analyses/placement_refs/nosto_ref_clade_assignments.csv") %>%
  filter(str_detect(dna_id, ".*fa"))
# Blast output colnames
blast_colnames <- c("qseqid", "sseqid", "length", "nident", "gaps", "pident")
# Load blast output
nostoc_corrected_blast <- read_delim("analyses/environmental_sequencing/ont_trial/placement/nostoc/blast/all_nostoc_consensus_blast.txt",
                                     col_names = blast_colnames)
# Get IDS for high identity matches
nostoc_corrected_blast_99 <- nostoc_corrected_blast %>%
  filter(pident >= 99) %>%
  left_join(nostoc_ref_clade_assignments, by = c("sseqid" = "dna_id")) %>%
  mutate(read_pool = str_remove(qseqid, "_nostocales.*")) %>%
  mutate(read_pool = factor(read_pool,levels = c("all", "s10", "s10_t", "s10_b", "s12", "s12_t", "s13", "s13_t", "s15", "s15_t")))
# Plot results
section_diversity_plot <- nostoc_corrected_blast_99 %>%
  ggplot(aes(read_pool)) +
  geom_bar(aes(fill = factor(section))) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Read pool", y = "No. of consensus seqs", fill = "Nostoc section") +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.background = element_blank())
ggsave("docs/ont_trial_sections.pdf")
