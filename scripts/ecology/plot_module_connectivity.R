#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load required libraries
library(tidyverse)
library(bipartite)

# Function to calculate normalized P_ij
 compute_Pij_normalized <- function(adj_matrix, row_modules, col_modules) {
  # Total weight in the matrix
  w <- sum(adj_matrix)
  
  # Compute row and column sums
  row_sums <- rowSums(adj_matrix)
  col_sums <- colSums(adj_matrix)
  
  # Initialize a matrix to store Pij values
  n_rows <- length(unique(row_modules))
  n_cols <- length(unique(col_modules))
  P_ij_normalized <- matrix(0, nrow = n_rows, ncol = n_cols)
  
  # Loop over row and column modules
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      if (i != j) {
        # Find rows and columns in the respective modules
        rows_in_module_a <- which(row_modules == i)
        cols_in_module_a <- which(col_modules == j)
        rows_in_module_b <- which(row_modules == j)
        cols_in_module_b <- which(col_modules == i)
        
        # Calculate the observed sum for this pair of modules
        observed_sum_a <- sum(adj_matrix[rows_in_module_a, cols_in_module_a])
        observed_sum_b <- sum(adj_matrix[rows_in_module_b, cols_in_module_b])
        
        # Calculate the maximum possible interactions for this pair of modules
        #max_interactions_a <- sum(row_sums[rows_in_module_a]) * sum(col_sums[cols_in_module_a]) / w
        #max_interactions_b <- sum(row_sums[rows_in_module_b]) * sum(col_sums[cols_in_module_b]) / w

        #P_ij_normalized_a <- observed_sum_a / max_interactions_a
        #P_ij_normalized_b <- observed_sum_b / max_interactions_b

        #P_ij_normalized[i, j] <- (P_ij_normalized_a + P_ij_normalized_b)

        # Calculate the normalized Pij for this pair of modules
        P_ij_normalized[i, j] <- sum(observed_sum_a, observed_sum_b) / w
      } else {
        # Find rows and columns in the respective modules
        rows_in_module <- which(row_modules == i)
        cols_in_module <- which(col_modules == j)
      
        # Calculate the observed sum for this pair of modules
        observed_sum <- sum(adj_matrix[rows_in_module, cols_in_module])
      
        # Calculate the maximum possible interactions for this pair of modules
        max_interactions <- sum(row_sums[rows_in_module]) * sum(col_sums[cols_in_module]) / w
      
        # Calculate the normalized Pij for this pair of modules
        P_ij_normalized[i, j] <- observed_sum / max_interactions
      }
    }
  }
  
  return(P_ij_normalized)
}



##### CALCULATE MODULE CONNECTIVITY TEST #####

adj_matrix <- matrix(c(1, 1, 0, 0, 
                       1, 1, 0, 0, 
                       0, 0, 1, 1,
                       0, 0, 1, 1), 
                     nrow = 4, byrow = TRUE)
rownames(adj_matrix) <- c("A", "B", "C", "D")  # Row names (e.g., species)
colnames(adj_matrix) <- c("X", "Y", "Z", "W") # Column names (e.g., samples)

row_modules <- c(A = 1, B = 1, C = 2, D = 2) # Module 1 for rows "A" and "B", Module 2 for "C".

col_modules <- c(X = 1, Y = 1, Z = 2, W = 2) # Module 1 for columns "X" and "Y", Module 2 for "Z" and "W".

# Example usage with simulated data
# adj_matrix: bipartite adjacency matrix
# row_modules and col_modules: module assignments for rows and columns
P_ij_normalized <- compute_Pij_normalized(adj_matrix, row_modules, col_modules)


##### CALCULATE MODULE CONNECTIVITY #####

# Load Peltigera-Nostoc interaction matrix
load("analyses/ecology/peltigera_interaction_matrix.RData")

# Get module objects for the networks
set.seed(12345) # To match modules in Fig. 1
peltigera_modules <- bipartite::metaComputeModules(abmi_matrix_peltigera, N = 5, method = "Beckett")

# Get module partitions 
peltigera_module_partitions <- bipartite::module2constraints(peltigera_modules)

# Get row and column partitions
peltigera_row_partitions <- peltigera_module_partitions[1:nrow(abmi_matrix_peltigera)]
peltigera_col_partitions <- peltigera_module_partitions[(nrow(abmi_matrix_peltigera)+1):(nrow(abmi_matrix_peltigera)+ncol(abmi_matrix_peltigera))]
names(peltigera_row_partitions) <- rownames(abmi_matrix_peltigera)
names(peltigera_col_partitions) <- colnames(abmi_matrix_peltigera)

# Calculate module connectivity
module_connectivity <- compute_Pij_normalized(abmi_matrix_peltigera, peltigera_row_partitions, peltigera_col_partitions)


##### CALCULATE MODULE COOCCURRENCE #####

# Load ABMI data
abmi_data <- read_csv("data/tables/abmi_id_data.csv")

# Load modeule assignments
load("analyses/ecology/peltigera_module_assignments.RData")

# Module records
module_sites <- abmi_data %>%
  filter(str_detect(mycobiont_molecular_id, "Peltigera") & !is.na(nostoc_otu)) %>%
  left_join(peltigera_module_assignments, by = "mycobiont_molecular_id") %>%
  left_join(nostoc_module_assignments, by = "nostoc_otu") %>%
  mutate(module = case_when(module.x == module.y ~ paste0("module", module.x),
                                .default = NA)) %>%
  filter(!is.na(module)) %>%
  distinct(site_year, module)

# Create a cross-tabulation of site_year and module
cross_tab <- module_sites %>%
  count(site_year, module) %>%
  spread(key = module, value = n, fill = 0)

# Generate the symmetric matrix
module_names <- colnames(cross_tab)[-1]  # Exclude site_year column
num_modules <- length(module_names)
symmetric_matrix <- matrix(0, nrow = num_modules, ncol = num_modules, dimnames = list(module_names, module_names))

for (i in 1:num_modules) {
  for (j in i:num_modules) {
    module_i <- module_names[i]
    module_j <- module_names[j]
    cooccurrences <- sum(cross_tab[[module_i]] > 0 & cross_tab[[module_j]] > 0)
    symmetric_matrix[i, j] <- cooccurrences
    symmetric_matrix[j, i] <- cooccurrences  # Ensure symmetry
  }
}

# Convert the upper triangle of the symmetric matrix to a data frame
module_cooccurrence <- data.frame(
  module_row = character(),
  module_column = character(),
  cooccurrences = integer()
)

for (i in 1:(num_modules - 1)) {
  for (j in (i + 1):num_modules) {
    module_cooccurrence <- rbind(module_cooccurrence, data.frame(
      module_row = module_names[i],
      module_column = module_names[j],
      cooccurrences = symmetric_matrix[i, j]
    ))
  }
}

module_cooccurrence <- module_cooccurrence %>%
  mutate(
    module_row_index = as.numeric(str_replace(module_row, "module", "")),
    module_column_index = as.numeric(str_replace(module_column, "module", ""))
  ) %>%
  mutate(connectivity = map2_dbl(module_row_index, module_column_index, ~ module_connectivity[.x, .y])) %>%
  select(-module_row_index, -module_column_index) %>%
  mutate(
    process = case_when(
      str_detect(module_row, "2|3|4") & str_detect(module_column, "2|3|4") ~ "env",
      .default = "other"
    )
  )


##### PLOT MODULE CONNECTIVITY VS COOCCURRENCE #####

# Define a custom theme for plots
custom_theme <- theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), 
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 12)
    )

# Make a scatter plot of module connectivity and co-occurrence and color by process
module_connectivity <- module_cooccurrence %>%
  ggplot(aes(x = cooccurrences, y = connectivity, color = process)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("env" = "#892b46", "other" = "gray50"),
                      labels = c("env" = "Modules C, D, and E", "other" = "Others")) +
  geom_line(data = module_cooccurrence %>% 
    filter(process == "env") %>% 
    arrange(cooccurrences), aes(group = 1), color = "#892b46") +
  labs(x = "No. of cooccurring sites", y = "Inter-module connectivity") +
  custom_theme +
  theme(legend.position = "top")  +
  guides(color = guide_legend(title = ""))

 ggsave(module_connectivity, filename = "documents/plots/module_connectivity.pdf",
    width = 10, height = 9, units = "cm")
