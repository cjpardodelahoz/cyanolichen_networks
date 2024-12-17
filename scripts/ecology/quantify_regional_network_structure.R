#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load required libraries
library(tidyverse)
library(bipartite)
library(data.table)
library(ggnewscale)

# Load custom functions

# Function to compute the metrics for the full network starting from the interaction matrix
quantify_network_topology <- function(interaction_matrix, network_name) {
    # Simulate null matrices for modularity test
    free_null_matrices <- bipartite::nullmodel(interaction_matrix, method = "vaznull", N = 999)
    # Compute modularity for the free null matrices
    free_null_modularity <- sapply(free_null_matrices, function(x) bipartite::metaComputeModules(x, method = "Beckett")@likelihood)
    # Compute modularity for the network and get module partitions
    observed_modularity <- bipartite::metaComputeModules(interaction_matrix, N = 5, method = "Beckett")
    # Comute expected modularity with the free null model
    expected_modularity <- mean(free_null_modularity)
    # Compute modularity z-score
    z_score_modularity <- (observed_modularity@likelihood - expected_modularity) / sd(free_null_modularity)
    # Compute modularity p-value
    p_value_modularity <- 1 / rank(c(observed_modularity@likelihood, free_null_modularity))[1]
    # If the interaction matrix is modular, continue with low-level nestedness analyses
    if (p_value_modularity < 0.5 & z_score_modularity > 2.0) {
        # Get module partitions
        module_partitions <- bipartite::module2constraints(observed_modularity)
        module_row_partitions <- module_partitions[1:nrow(interaction_matrix)]
        module_col_partitions <- module_partitions[(nrow(interaction_matrix)+1):(nrow(interaction_matrix)+ncol(interaction_matrix))]
        # Compute nestedness for the network with the WNODA method
        observed_nestedness <- bipartite::nest.smdm(interaction_matrix, 
            constraints = module_partitions, 
            weighted = TRUE,
            decreasing = "abund")
        # Simulate null matrics under the equiprobable and proportional null models for nestedness test
        equiprobable_null_matrices <- bipartite::restrictednull(interaction_matrix,
            Prior.Pij = "equiprobable",
            conditional.level = "modules",
            R.partitions = module_row_partitions,
            C.partitions = module_col_partitions,
            N = 999)
        proportional_null_matrices <- bipartite::restrictednull(interaction_matrix,
            Prior.Pij = "degreeprob",
            conditional.level = "modules",
            R.partitions = module_row_partitions,
            C.partitions = module_col_partitions,
            N = 999)
        # Compute modularity for restricted null matrices and get module partitions
        equiprobable_null_modularity <- sapply(equiprobable_null_matrices, function(x) bipartite::metaComputeModules(x, method = "Beckett"))
        proportional_null_modularity <- sapply(proportional_null_matrices, function(x) bipartite::metaComputeModules(x, method = "Beckett"))
        equiprobable_null_module_partitions <- lapply(equiprobable_null_modularity, function(x) bipartite::module2constraints(x))
        proportional_null_module_partitions <- lapply(proportional_null_modularity, function(x) bipartite::module2constraints(x))
        # Compute nestedness for the equiprobable and proportional null matrices
        equiprobable_null_nestedness <- mapply(function(x, y) bipartite::nest.smdm(x, constraints = y, weighted = TRUE, decreasing = "abund"), 
            equiprobable_null_matrices, equiprobable_null_module_partitions)
        proportional_null_nestedness <- mapply(function(x, y) bipartite::nest.smdm(x, constraints = y, weighted = TRUE, decreasing = "abund"),
            proportional_null_matrices, proportional_null_module_partitions)
        # Compute metrics for the network and store in a table
        network_metrics <- tibble::tibble(
            # Add the network name
            network = network_name,
            #richness = sum(rowSums(interaction_matrix) > 0),
            #connectance = sum(interaction_matrix > 0) / (sum(rowSums(interaction_matrix) > 0) * sum(colSums(interaction_matrix) > 0)),
            # Get modularity using the Beckett method
            observed_modularity = observed_modularity@likelihood,
            # Compute expected modularity with the free null model
            expected_modularity = mean(free_null_modularity),
            # Compute modularity z-score
            z_score_modularity = (observed_modularity - expected_modularity) / sd(free_null_modularity),
            # Compute modularity p-value
            p_value_modularity = 1 / rank(c(observed_modularity, free_null_modularity))[1],
            # Get full matrix nestedness (WNODA)
            observed_wnoda_full = observed_nestedness$WNODAmatrix,
            # Get same-module (SM) and different-module (DM) nestedness
            observed_wnoda_sm = observed_nestedness$WNODA_SM_matrix,
            observed_wnoda_dm = observed_nestedness$WNODA_DM_matrix,
            # Get expected nestedness for  SM, and DM using the equiprobable and proportional null models
            expected_wnoda_sm_equiprobable = mean(equiprobable_null_nestedness["WNODA_SM_matrix", ] %>% unlist()),
            expected_wnoda_dm_equiprobable = mean(equiprobable_null_nestedness["WNODA_DM_matrix", ] %>% unlist()),
            expected_wnoda_sm_proportional = mean(proportional_null_nestedness["WNODA_SM_matrix", ] %>% unlist()),
            expected_wnoda_dm_proportional = mean(proportional_null_nestedness["WNODA_DM_matrix", ] %>% unlist()),
            # Compute SM and DM z-scores
            z_score_wnoda_sm_equiprobable = (observed_wnoda_sm - expected_wnoda_sm_equiprobable) / sd(equiprobable_null_nestedness["WNODA_SM_matrix", ] %>% unlist()),
            z_score_wnoda_dm_equiprobable = (observed_wnoda_dm - expected_wnoda_dm_equiprobable) / sd(equiprobable_null_nestedness["WNODA_DM_matrix", ] %>% unlist()),
            z_score_wnoda_sm_proportional = (observed_wnoda_sm - expected_wnoda_sm_proportional) / sd(proportional_null_nestedness["WNODA_SM_matrix", ] %>% unlist()),
            z_score_wnoda_dm_proportional = (observed_wnoda_dm - expected_wnoda_dm_proportional) / sd(proportional_null_nestedness["WNODA_DM_matrix", ] %>% unlist()),
            # Calculate p- values for SM nestedness with the equiprobable and proportional null models
            p_value_wnoda_sm_equiprobable = 1 / rank(c(observed_wnoda_sm, equiprobable_null_nestedness["WNODA_SM_matrix", ] %>% unlist()))[1],
            p_value_wnoda_dm_equiprobable = 1 / rank(c(observed_wnoda_dm, equiprobable_null_nestedness["WNODA_DM_matrix", ] %>% unlist()))[1],
            p_value_wnoda_sm_proportional = 1 / rank(c(observed_wnoda_sm, proportional_null_nestedness["WNODA_SM_matrix", ] %>% unlist()))[1],
            p_value_wnoda_dm_proportional = 1 / rank(c(observed_wnoda_dm, proportional_null_nestedness["WNODA_DM_matrix", ] %>% unlist()))[1]
        )
    # If the interaction matrix is not modular, calculate only observed nestedness and return a table with observed_wnoda_full and the rest columns with NA values
    } else {
        # Compute nestedness for the network with the WNODA method
        observed_nestedness <- bipartite::nest.smdm(interaction_matrix, 
            weighted = TRUE,
            decreasing = "abund")
        # Compute metrics for the network and store in a table
        network_metrics <- tibble::tibble(
            # Add the network name
            network = network_name,
            #richness = sum(rowSums(interaction_matrix) > 0),
            #connectance = sum(interaction_matrix > 0) / (sum(rowSums(interaction_matrix) > 0) * sum(colSums(interaction_matrix) > 0)),
            # Get modularity using the Beckett method
            observed_modularity = observed_modularity@likelihood,
            # Compute expected modularity with the free null model
            expected_modularity = mean(free_null_modularity),
            # Compute modularity z-score
            z_score_modularity = (observed_modularity - expected_modularity) / sd(free_null_modularity),
            # Compute modularity p-value
            p_value_modularity = 1 / rank(c(observed_modularity, free_null_modularity))[1],
            # Get full matrix nestedness (WNODA)
            observed_wnoda_full = observed_nestedness$WNODAmatrix,
            # Get same-module (SM) and different-module (DM) nestedness
            observed_wnoda_sm = NA,
            observed_wnoda_dm = NA,
            # Get expected nestedness for  SM, and DM using the equiprobable and proportional null models
            expected_wnoda_sm_equiprobable = NA,
            expected_wnoda_dm_equiprobable = NA,
            expected_wnoda_sm_proportional = NA,
            expected_wnoda_dm_proportional = NA,
            # Compute SM and DM z-scores
            z_score_wnoda_sm_equiprobable = NA,
            z_score_wnoda_dm_equiprobable = NA,
            z_score_wnoda_sm_proportional = NA,
            z_score_wnoda_dm_proportional = NA,
            # Calculate p- values for SM nestedness with the equiprobable and proportional null models
            p_value_wnoda_sm_equiprobable = NA,
            p_value_wnoda_dm_equiprobable = NA,
            p_value_wnoda_sm_proportional = NA,
            p_value_wnoda_dm_proportional = NA
        )
    }    
    return(network_metrics)
}


##### ID DATA TO INTERACTION MATRICES #####

# Load ABMI dataset
abmi_data <- read_csv("data/tables/abmi_id_data.csv") #%>% 
    #select(1:Collector) %>%
    # Rename column names as lower case and replace spaces with underscores
    #rename_all(tolower) %>%
    #rename_all(~str_replace_all(., " ", "_")) %>%
    # When species complex is "3.10a", replace section with "3.10"
    #mutate(section = case_when(
    #    species_complex == "3.10a" ~ "3.10",
    #    .default = section
    #)) %>%
    # If the section column is not NA, append the string "section" to the section column
    #mutate(section = case_when(
    #    !is.na(section) ~ paste0("section_", section),
    #    .default = section
    #)) %>%
    # If the species_complex column is not NA, append the string "species_complex" to the species_complex column
    #mutate(species_complex = case_when(
    #    !is.na(species_complex) ~ paste0("species_complex_", species_complex),
    #    .default = species_complex
    #)) %>%
    # Coalesce the columns phylogroup and species complex into a new column called "nostoc_otu"
    #mutate(nostoc_otu = coalesce(phylogroup, species_complex))

# Get full interaction matrix including all species
abmi_matrix_full <- abmi_data %>%
    # Filter to records with IDs for both mycobiont and nostoc
    filter(!is.na(mycobiont_molecular_id), !is.na(nostoc_otu)) %>%
    # Select only the columns mycobiont_molecular_id and nostoc_otu
    select(mycobiont_molecular_id, nostoc_otu) %>%
    # Count the number of times each pair of nostoc_otu and mycobiont_molecular_id occurs
    count(nostoc_otu, mycobiont_molecular_id) %>%
    # Convert the data to a matrix with the columns nostoc_otu and mycobiont_molecular_id as rows and the count as the value
    spread(nostoc_otu, n, fill = 0) %>%
    # Make mycobiont_molecular_id the row names
    column_to_rownames("mycobiont_molecular_id") %>%
    as.matrix()

# Get interaction matrix for Peltigera species
abmi_matrix_peltigera <- abmi_data %>%
    # Filter the data to include mycobiont_molecular_id that contains "Peltigera" in the string and remove rows with NA in the nostoc_otu column
    filter(str_detect(mycobiont_molecular_id, "Peltigera") & !is.na(nostoc_otu)) %>%
    # Select only the columns mycobiont_molecular_id and nostoc_otu
    select(mycobiont_molecular_id, nostoc_otu) %>%
    # Count the number of times each pair of nostoc_otu and mycobiont_molecular_id occurs
    count(nostoc_otu, mycobiont_molecular_id) %>%
    # Convert the data to a matrix with the columns nostoc_otu and mycobiont_molecular_id as rows and the count as the value
    spread(nostoc_otu, n, fill = 0) %>%
    # Make mycobiont_molecular_id the row names
    column_to_rownames("mycobiont_molecular_id") %>%
    as.matrix()

# Save the interaction matrix objects to a file
if(!dir.exists("analyses/ecology")) {
    dir.create("analyses/ecology", recursive = TRUE)
}
save(abmi_matrix_full, file = "analyses/ecology/full_interaction_matrix.RData")
save(abmi_matrix_peltigera, file = "analyses/ecology/peltigera_interaction_matrix.RData")


##### REGIONAL NETWORK STRUCTURE ANALYSES #####

# Compute metrics for the Peltigera network and the full network
peltigera_network_metrics <- quantify_network_topology(abmi_matrix_peltigera, "peltigera")
full_network_metrics <- quantify_network_topology(abmi_matrix_full, "full")

# Save the network metric objects to a file
if(!dir.exists("analyses/ecology")) {
    dir.create("analyses/ecology", recursive = TRUE)
}
save(peltigera_network_metrics, file = "analyses/ecology/peltigera_regional_network_metrics.RData")
save(full_network_metrics, file = "analyses/ecology/full_regional_network_metrics.RData")

# Table with reported metrics
regional_network_metrics <- bind_rows(peltigera_network_metrics, full_network_metrics) %>%
    select(-contains("wnoda_sm_proportional"),
            -contains("wnoda_dm_proportional"),
            -observed_wnoda_full) %>%
    rename_with(~ str_replace(., "_equiprobable", ""), contains("wnoda")) %>%
    pivot_longer(
        cols = -network,
        names_to = c(".value", "metric"),
        names_pattern = "(observed|expected|z_score|p_value)_(.*)")


# Save the network metrics table to a CSV file
write_csv(regional_network_metrics, "documents/tables/regional_network_metrics.csv")
