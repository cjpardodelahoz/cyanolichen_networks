#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(treeio)
library(ape)

# Parse placement tree

# Nostocales outgroup
outgroup <- c("Scytonema_millei_VB511283", "Chroococcidiopsis_thermalis_PCC_7203", "Gloeocapsa_sp_PCC_7428", "Synechocystis_sp_PCC_7509", "QUERY___s10_t_nostocales_cluster55")
# Load trees with EPA placements
nostocales_corrected_tree <- read.tree("analyses/environmental_sequencing/ont_trial/placement/nostocales/trees/RAxML_labelledTree.corrected_nostocales_epa_result") %>%
  root.phylo(outgroup, resolve.root = T)
# Save placement tree in newick to open with FigTree
write.tree(nostocales_corrected_tree, "analyses/environmental_sequencing/ont_trial/placement/nostocales/trees/nostocales_corrected_placement.tree")

# Get IDs of sequences that fell within Nostoc

# Subset tree to Nostoc node
nostoc_node <- MRCA(nostocales_corrected_tree, c("Nostoc_sp_JC1668", "Nostoc_sp_KVJ20"))
nostoc_16s_queries <- tree_subset(nostocales_corrected_tree, node = nostoc_node, 
                                  levels_back = 0)$tip.label %>%
  as_tibble() %>%
  filter(str_detect(value, "QUERY")) %>%
  pull(value) %>%
  str_remove("QUERY___") %>%
  write("analyses/environmental_sequencing/ont_trial/placement/nostocales/nostoc_corrected_labels.txt")
