#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load required libraries
library(tidyverse)
library(ggnewscale)


##### PLOT THE NETWORKS HIGHLIGHTING COMPOUND STRUCTURE #####

# Load ABMI data and interaction matrix
abmi_data <- read_csv("data/tables/abmi_id_data.csv") 
load("analyses/ecology/peltigera_interaction_matrix.RData")

# Get module objects for the networks
set.seed(12345)
peltigera_modules <- bipartite::metaComputeModules(abmi_matrix_peltigera, N = 5, method = "Beckett")
#full_modules <- bipartite::metaComputeModules(abmi_matrix_full, N = 5, method = "Beckett")

# Get module partitions 
peltigera_module_partitions <- bipartite::module2constraints(peltigera_modules)
#full_module_partitions <- bipartite::module2constraints(full_module_partitions)

# Get row and column partitions
peltigera_row_partitions <- peltigera_module_partitions[1:nrow(abmi_matrix_peltigera)]
peltigera_col_partitions <- peltigera_module_partitions[(nrow(abmi_matrix_peltigera)+1):(nrow(abmi_matrix_peltigera)+ncol(abmi_matrix_peltigera))]
#full_row_partitions <- full_module_partitions[1:nrow(abmi_matrix_full)]
#full_col_partitions <- full_module_partitions[(nrow(abmi_matrix_full)+1):(nrow(abmi_matrix_full)+ncol(abmi_matrix_full))]

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
                                        levels = c("section24", "section31", "section32", "section34", "section35", "section36",
                                            "section37", "section39", "section310", "section311", "section312"))

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


##### ALTERNATIVE PLOT #####

# Load ABMI data and interaction matrix
abmi_data <- read_csv("data/tables/abmi_id_data.csv") 
load("analyses/ecology/peltigera_interaction_matrix.RData")
abmi_matrix_peltigera <- t(abmi_matrix_peltigera)

# Get module objects for the networks
set.seed(12345)
peltigera_modules <- bipartite::metaComputeModules(abmi_matrix_peltigera, N = 5, method = "Beckett")
#full_modules <- bipartite::metaComputeModules(abmi_matrix_full, N = 5, method = "Beckett")

# Get module partitions 
peltigera_module_partitions <- bipartite::module2constraints(peltigera_modules)
#full_module_partitions <- bipartite::module2constraints(full_module_partitions)

# Get row and column partitions
peltigera_row_partitions <- peltigera_module_partitions[1:nrow(abmi_matrix_peltigera)]
peltigera_col_partitions <- peltigera_module_partitions[(nrow(abmi_matrix_peltigera)+1):(nrow(abmi_matrix_peltigera)+ncol(abmi_matrix_peltigera))]
#full_row_partitions <- full_module_partitions[1:nrow(abmi_matrix_full)]
#full_col_partitions <- full_module_partitions[(nrow(abmi_matrix_full)+1):(nrow(abmi_matrix_full)+ncol(abmi_matrix_full))]

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
    rownames_to_column(var = "nostoc_otu") %>%
    pivot_longer(cols = -nostoc_otu, names_to = "mycobiont_molecular_id", values_to = "interaction_frequency") %>%
    left_join(distinct(abmi_data, mycobiont_molecular_id, peltigera_section), by = "mycobiont_molecular_id") %>%
    left_join(distinct(abmi_data, nostoc_otu, section), by = "nostoc_otu") %>%
    select(mycobiont_molecular_id, nostoc_otu, interaction_frequency, section,peltigera_section)

# Set order for Nostoc section legend
abmi_peltigera_matrix_df$section <- factor(abmi_peltigera_matrix_df$section, 
                                        levels = c("section24", "section31", "section32", "section34", "section35", "section36",
                                            "section37", "section39", "section310", "section311", "section312"))

# Convert Rows and Columns to factors with levels matching the orded of the sorted matrix
abmi_peltigera_matrix_df$mycobiont_molecular_id <- factor(abmi_peltigera_matrix_df$mycobiont_molecular_id, levels = rev(colnames(peltigera_matrix_compound$matrix)))
abmi_peltigera_matrix_df$nostoc_otu <- factor(abmi_peltigera_matrix_df$nostoc_otu, levels = rownames(peltigera_matrix_compound$matrix))

# Plot the matrix with compound topology
# Start with the base plot
peltigera_compound_plot <- ggplot(abmi_peltigera_matrix_df, aes(x = nostoc_otu, y = mycobiont_molecular_id, fill = interaction_frequency)) +
            geom_tile(height = 0.8, width = 0.8) +
            scale_fill_gradient(low = "white", high = "black",
                guide = guide_colorbar(title = "Interaction frequency",
                                    title.position = "top",
                                    order = 1)) +
            new_scale_fill() +
            geom_tile(data = abmi_peltigera_matrix_df, 
                aes(y = 66.7, x = nostoc_otu, fill = section), height = 0.8) +
            scale_fill_brewer(palette = "Set3", guide = "none") +
            new_scale_fill() +
            geom_tile(data = abmi_peltigera_matrix_df, 
                aes(y = mycobiont_molecular_id, x = -0.7, fill = peltigera_section), height = 1, width = 0.8) +
            scale_fill_brewer(palette = "Set2", guide = "none") +
            labs(y = "Peltigera species", x = "Nostoc OTU") +
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
ggsave(plot = peltigera_compound_plot, filename = "documents/plots/peltigera_compound_alternative.pdf", 
    width = 8, height = 17, units = "cm")


##### MAKE TABLES FOR MAPPING MODULES #####

# Make tables with the module assignments for mapping
peltigera_module_assignments <- data.frame(mycobiont_molecular_id = rev(colnames(peltigera_matrix_compound$matrix)), 
    module = peltigera_y_partitions)
nostoc_module_assignments <- data.frame(nostoc_otu = rownames(peltigera_matrix_compound$matrix), 
    module = peltigera_x_partitions)

# Save the tables as RData files
save(peltigera_module_assignments, nostoc_module_assignments, file = "analyses/ecology/peltigera_module_assignments.RData")
