#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)

# Load Reference IDs from ABMI sequences
abmi_ref_ids <- read_csv("data/tables/nostoc_datas1c_includes_ined.csv") %>%
    select(c("DNA ID", "Section", "Species complex", "Phylogroup", "Mycobiont molecular ID"))
# Load blast results
its_blast_results <- read_delim("analyses/lichen_sequencing/ont_lichen/id/its/its_haplotypes_blast.tsv", delim = "\t", col_names = F) %>%
    group_by(X1) %>%
    filter(X4 == max(X4)) %>%
    ungroup()
rbclx_blast_results <- read_delim("analyses/lichen_sequencing/ont_lichen/id/rbclx/rbclx_haplotypes_blast.tsv", delim = "\t", col_names = F) %>%
    group_by(X1) %>%
    filter(X4 == max(X4)) %>%
    ungroup()
colnames(its_blast_results) <- c("query", "its_subject", "its_identity", "its_alignment_length", "its_mismatches")
colnames(rbclx_blast_results) <- c("query", "rbclx_subject", "rbclx_identity", "rbclx_alignment_length", "rbclx_mismatches")
# Load read counts
its_read_counts <- read_csv("analyses/lichen_sequencing/ont_lichen/consensus/its_read_counts.csv")
rbclx_read_counts <- read_csv("analyses/lichen_sequencing/ont_lichen/consensus/rbclx_read_counts.csv")
# Load haplotype cluster information
its_haplotype_clusters <- read_delim("analyses/lichen_sequencing/ont_lichen/id/its/its_haplotypes.tsv", delim = "\t", col_names = F) %>%
    select(c(X1, X2, X3, X9)) %>%
    filter(X1 != "C")
colnames(its_haplotype_clusters) <- c("its_record_type", "its_cluster", "its_length", "dna_code")
rbclx_haplotype_clusters <- read_delim("analyses/lichen_sequencing/ont_lichen/id/rbclx/rbclx_haplotypes.tsv", delim = "\t", col_names = F) %>%
    select(c(X1, X2, X3, X9)) %>%
    filter(X1 != "C")
colnames(rbclx_haplotype_clusters) <- c("rbclx_record_type", "rbclx_cluster", "rbclx_length", "dna_code")

# Join voucher with haplotype clusters, read counts, and blast results
voucher <- read_csv("documents/tables/lichen_voucher_jm.csv") %>%
    left_join(its_haplotype_clusters, by = "dna_code") %>%
    left_join(rbclx_haplotype_clusters, by = "dna_code") %>%
    left_join(its_read_counts, by = c("dna_code" = "sample")) %>%
    left_join(rbclx_read_counts, by = c("dna_code" = "sample")) %>%
    left_join(its_blast_results, by = c("dna_code" = "query")) %>%
    left_join(rbclx_blast_results, by = c("dna_code" = "query")) %>%
    left_join(abmi_ref_ids, by = c("its_subject" = "DNA ID")) %>%
    select(-c("Section", "Species complex", "Phylogroup")) %>%
    left_join(abmi_ref_ids, by = c("rbclx_subject" = "DNA ID"))
# Write voucher table to file to curate IDs
write_csv(voucher, "documents/tables/lichen_voucher_jm_curated_v0.csv")
