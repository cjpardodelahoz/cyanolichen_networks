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
    free_null_matrices <- bipartite::nullmodel(interaction_matrix, method = "r2d", N = 999)
    # Compute modularity for the free null matrices
    free_null_modularity <- sapply(free_null_matrices, function(x) bipartite::metaComputeModules(x, method = "Beckett")@likelihood)
    # Compute modularity for the network and get module partitions
    observed_modularity <- bipartite::metaComputeModules(interaction_matrix, N = 5, method = "Beckett")
    # Comute expected modularity with the free null model
    expected_modularity <- mean(free_null_modularity)
    # Compute modularity z-score
    modularity_z_score <- (observed_modularity@likelihood - expected_modularity) / sd(free_null_modularity)
    # Compute modularity p-value
    modularity_p_value <- 1 / rank(c(observed_modularity@likelihood, free_null_modularity))[1]
    # If the interaction matrix is modular, continue with low-level nestedness analyses
    if (modularity_p_value < 0.01 & modularity_z_score > 2.0) {
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
            modularity_z_score = (observed_modularity - expected_modularity) / sd(free_null_modularity),
            # Compute modularity p-value
            modularity_p_value = 1 / rank(c(observed_modularity, free_null_modularity))[1],
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
            z_score_sm_equiprobable = (observed_wnoda_sm - expected_wnoda_sm_equiprobable) / sd(equiprobable_null_nestedness["WNODA_SM_matrix", ] %>% unlist()),
            z_score_dm_equiprobable = (observed_wnoda_dm - expected_wnoda_dm_equiprobable) / sd(equiprobable_null_nestedness["WNODA_DM_matrix", ] %>% unlist()),
            z_score_sm_proportional = (observed_wnoda_sm - expected_wnoda_sm_proportional) / sd(proportional_null_nestedness["WNODA_SM_matrix", ] %>% unlist()),
            z_score_dm_proportional = (observed_wnoda_dm - expected_wnoda_dm_proportional) / sd(proportional_null_nestedness["WNODA_DM_matrix", ] %>% unlist())
        )
    # If the interaction matrix is not modular, calculate only observed nestedness and return a table with observed_wnoda_full and the rest columns with NA values
    } else {
        # Compute nestedness for the network with the WNODA method
        observed_nestedness <- bipartite::nest.smdm(interaction_matrix, 
            constraints = module_partitions, 
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
            modularity_z_score = (observed_modularity - expected_modularity) / sd(free_null_modularity),
            # Compute modularity p-value
            modularity_p_value = 1 / rank(c(observed_modularity, free_null_modularity))[1],
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
            z_score_sm_equiprobable = NA,
            z_score_dm_equiprobable = NA,
            z_score_sm_proportional = NA,
            z_score_dm_proportional = NA
        )
    }    
    return(network_metrics)
}
##### ID DATA TO INTERACTION MATRICES #####

# Load ABMI dataset from Pardo-De la Hoz et al. (2024)-Data S1C
abmi_data <- read_csv("data/tables/nostoc_datas1c_includes_ined.csv") %>% 
    select(1:Collector) %>%
    # Rename column names as lower case and replace spaces with underscores
    rename_all(tolower) %>%
    rename_all(~str_replace_all(., " ", "_")) %>%
    # When species complex is "3.10a", replace section with "3.10"
    mutate(section = case_when(
        species_complex == "3.10a" ~ "3.10",
        .default = section
    )) %>%
    # If the section column is not NA, append the string "section" to the section column
    mutate(section = case_when(
        !is.na(section) ~ paste0("section_", section),
        .default = section
    )) %>%
    # If the species_complex column is not NA, append the string "species_complex" to the species_complex column
    mutate(species_complex = case_when(
        !is.na(species_complex) ~ paste0("species_complex_", species_complex),
        .default = species_complex
    )) %>%
    # Coalesce the columns phylogroup and species complex into a new column called "nostoc_otu"
    mutate(nostoc_otu = coalesce(phylogroup, species_complex))

# Get full interaction matrix including all species
abmi_matrix_all <- abmi_data %>%
    # Filter the data to include only rows where mycobiont_molecular_id is not NA
    filter(!is.na(mycobiont_molecular_id)) %>%
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

##### REGIONAL NETWORK TOPOLOGY ANALYSES #####

# Compute metrics for the Peltigera network and the full network
peltigera_network_metrics <- quantify_network_topology(abmi_matrix_peltigera, "peltigera")
all_network_metrics <- quantify_network_topology(abmi_matrix_all, "all")
# Merge the network metrics into a single table
regional_network_metrics <- bind_rows(peltigera_network_metrics, all_network_metrics)
# Save the network metrics table to a CSV file
write_csv(network_metrics, "documents/tables/regional_network_metrics.csv")

##### PLOT THE NETWORKS HIGHLIGHTING COMPOUND STRUCTURE #####

# Get module objects for the networks
peltigera_modules <- bipartite::metaComputeModules(abmi_matrix_peltigera, N = 5, method = "Beckett")
all_modules <- bipartite::metaComputeModules(abmi_matrix_all, N = 5, method = "Beckett")
# Get module partitions 
peltigera_module_partitions <- bipartite::module2constraints(peltigera_module_partitions)
all_module_partitions <- bipartite::module2constraints(all_module_partitions)
# Get row and column partitions
peltigera_row_partitions <- peltigera_module_partitions[1:nrow(abmi_matrix_peltigera)]
peltigera_col_partitions <- peltigera_module_partitions[(nrow(abmi_matrix_peltigera)+1):(nrow(abmi_matrix_peltigera)+ncol(abmi_matrix_peltigera))]
all_row_partitions <- all_module_partitions[1:nrow(abmi_matrix_all)]
all_col_partitions <- all_module_partitions[(nrow(abmi_matrix_all)+1):(nrow(abmi_matrix_all)+ncol(abmi_matrix_all))]
# Sort matrices in compound topology
peltigera_matrix_compound <- bipartite::sortmatrix(abmi_matrix_peltigera, 
    topology = "compound", sort_by = "weights", 
    row_partitions = peltigera_row_partitions, col_partitions = peltigera_col_partitions,
    mod_similarity = TRUE)
# Get row and column partitions for the sorted matrix
peltigera_x_partitions <- peltigera_matrix_compound$row_partitions
peltigera_y_partitions <- rev(peltigera_matrix_compound$col_partitions)

# Reshape matrix data for ggplot
abmi_peltigera_matrix_df <- abmi_matrix_peltigera %>%
    as.data.frame() %>%
    rownames_to_column(var = "mycobiont_molecular_id") %>%
    pivot_longer(cols = -mycobiont_molecular_id, names_to = "nostoc_otu", values_to = "interaction_frequency") %>%
    left_join(distinct(abmi_data, mycobiont_molecular_id, peltigera_section), by = "mycobiont_molecular_id") %>%
    left_join(distinct(abmi_data, nostoc_otu, section), by = "nostoc_otu") %>%
    select(mycobiont_molecular_id, nostoc_otu, interaction_frequency, section,peltigera_section)
# Set order for Nostoc section legend
abmi_peltigera_matrix_df$section <- factor(abmi_peltigera_matrix_df$section, 
                                        levels = c("section_2.4", "section_3.1", "section_3.2", "section_3.4", "section_3.5", "section_3.6",
                                            "section_3.7", "section_3.9", "section_3.10", "section_3.11", "section_3.12"))
# Convert Rows and Columns to factors with levels matching the orded of the sorted matrix
abmi_peltigera_matrix_df$mycobiont_molecular_id <- factor(abmi_peltigera_matrix_df$mycobiont_molecular_id, levels = rownames(peltigera_matrix_compound$matrix))
abmi_peltigera_matrix_df$nostoc_otu <- factor(abmi_peltigera_matrix_df$nostoc_otu, levels = rev(colnames(peltigera_matrix_compound$matrix)))


# Plot the matrix with compound topology
# Start with the base plot
peltigera_compound_plot <- ggplot(abmi_peltigera_matrix_df, aes(x = mycobiont_molecular_id, y = nostoc_otu, fill = interaction_frequency)) +
            geom_tile(height = 0.8, width = 0.8) +
            scale_fill_gradient(low = "white", high = "black",
                guide = guide_colorbar(title = "Interaction frequency",
                                    title.position = "top",
                                    order = 1)) +
            new_scale_fill() +
            geom_tile(data = abmi_peltigera_matrix_df, 
                aes(x = -0.3, y = nostoc_otu, fill = section), height = 1, width = 0.4) +
            scale_fill_brewer(palette = "Set3",
                guide = guide_legend(title = "Nostoc section",
                                    title.position = "top",
                                    order = 2)) +
            new_scale_fill() +
            geom_tile(data = abmi_peltigera_matrix_df, 
                aes(x = mycobiont_molecular_id, y = -0.3, fill = peltigera_section), height = 0.4) +
            scale_fill_brewer(palette = "Set2",
                guide = guide_legend(title = "Peltigera section",
                                    title.position = "top",
                                    order = 3)) +
            labs(x = "Peltigera species", y = "Nostoc OTU") +
            theme_minimal() +
            annotate("rect",
                                xmin = 0,
                                xmax = length(peltigera_x_partitions) + 1,
                                ymin = 0,
                                ymax = length(peltigera_y_partitions) + 1,
                                color = "black", fill = NA, linetype = 1, size = 0.2) +
            theme(axis.text = element_blank(),
                panel.grid = element_blank(),
                legend.position = "top",
                legend.direction = "horizontal",
                legend.key.height = unit(2, "mm"),
                legend.key.width = unit(2, "mm"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 6),
                axis.title = element_text(size = 8))
# Get unique partition values
unique_partitions <- unique(c(peltigera_x_partitions, peltigera_y_partitions))
# Loop over unique partition values
for (i in unique_partitions) {
    # Get indices for rows and columns
    x_indices <- which(peltigera_x_partitions == i)
    y_indices <- which(peltigera_y_partitions == i)       
    # Add annotations for each partition
    peltigera_compound_plot <- peltigera_compound_plot+ annotate("rect", 
                                        xmin = min(x_indices) - 0.5, 
                                        xmax = max(x_indices) + 0.5, 
                                        ymin = min(y_indices) - 0.5, 
                                        ymax = max(y_indices) + 0.5, 
                                        color = "black", fill = NA, linetype = 1, size = 0.2) +
                                    annotate("segment",
                                        x = min(x_indices) - 0.5,
                                        xend = min(x_indices) - 0.5,
                                        y = 0,
                                        yend = min(y_indices) - 0.5,
                                        color = "gray70", linewidth = 0.2, linetype = 2) +
                                    annotate("segment",
                                        x = max(x_indices) + 0.5,
                                        xend = max(x_indices) + 0.5,
                                        y = max(y_indices) + 0.5,
                                        yend = length(peltigera_y_partitions) + 1,
                                        color = "gray70", linewidth = 0.2, linetype = 2) +
                                    annotate("segment",
                                        x = 0,
                                        xend = min(x_indices) - 0.5,
                                        y = min(y_indices) - 0.5,
                                        yend = min(y_indices) - 0.5,
                                        color = "gray70", linewidth = 0.2, linetype = 2) +
                                    annotate("segment",
                                        x = max(x_indices) + 0.5,
                                        xend = length(peltigera_x_partitions) + 1,
                                        y = min(y_indices) - 0.5,
                                        yend = min(y_indices) - 0.5,
                                        color = "gray70", linewidth = 0.2, linetype = 2)
}
# Save the plot
ggsave(plot = peltigera_compound_plot, filename = "documents/plots/peltigera_compound.pdf", 
    width = 17, height = 10, units = "cm")
  

#plotmatrix(peltigera_matrix_compound, binary = FALSE, border = TRUE, plot_labels = TRUE)




