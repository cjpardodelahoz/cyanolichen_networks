#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

library(tidyverse)

##### PREP OCCURRENCE AND SPECIFICITY DATA FOR PLOTTING #####

# Load ABMI ID data and filter to Peltigera-Nostoc
abmi_id_data <- read_csv("data/tables/abmi_id_data.csv") %>%
  filter(str_detect(mycobiont_molecular_id, "Peltigera") & !is.na(nostoc_otu))

# Load module assignments
load("analyses/ecology/peltigera_module_assignments.RData")

# Define custom colors for the modules
module_colors <- c("1" = "#8a708a", "2" = "#7e937b", "3" = "#3277b0", "4" = "#be8551", "5" = "gray70", "6" = "gray20")

# Summarize the table for mycobiont (Peltigera)
peltigera_summary <- abmi_id_data %>%
  group_by(mycobiont_molecular_id) %>%
  summarise(
    n_sites = n_distinct(site),
    n_partners = n_distinct(nostoc_otu),
    symbiont_type = "peltigera"
  ) %>%
  rename(symbiont_id = mycobiont_molecular_id) %>%
  left_join(peltigera_module_assignments, by = c("symbiont_id" = "mycobiont_molecular_id")) %>%
  mutate(color = module_colors[as.character(module)])

# Summarize the table for nostoc
nostoc_summary <- abmi_id_data %>%
  group_by(nostoc_otu) %>%
  summarise(
    n_sites = n_distinct(site),
    n_partners = n_distinct(mycobiont_molecular_id),
    symbiont_type = "nostoc"
  ) %>%
  rename(symbiont_id = nostoc_otu) %>%
  left_join(nostoc_module_assignments, by = c("symbiont_id" = "nostoc_otu")) %>%
  mutate(color = module_colors[as.character(module)])

##### PLOT OCCURRENCE AND SPECIFICITY DATA #####

# Define a custom theme for plots
custom_theme <- theme(
  panel.grid.major = element_blank(),  # Hide major grid lines
  panel.grid.minor = element_blank(),  # Hide minor grid lines
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), 
  axis.text = element_text(color = "black", size = 12),
  legend.text = element_text(color = "black", size = 12),
  legend.title = element_text(color = "black", size = 12),
  legend.position = "right",  # Position legend on the right
  legend.direction = "vertical"  # Make legend vertical
)

# Plot for Peltigera
peltigera_plot <- ggplot(peltigera_summary, aes(x = n_sites, y = n_partners)) +
  geom_point(aes(color = color), shape = 17) +
  geom_smooth(method = "lm", se = FALSE, color = "gray50", linetype = "dashed", size = 0.5) +  # Add dashed trend line
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +  # Set y scale breaks
  scale_color_identity() +  # Use the color column for point colors
  labs(x = "Number of sites",
       y = "Number of partners") +
  custom_theme

# Plot for Nostoc
nostoc_plot <- ggplot(nostoc_summary, aes(x = n_sites, y = n_partners)) +
  geom_point(aes(color = color)) +
  geom_smooth(method = "lm", se = FALSE, color = "gray50", linetype = "dashed", size = 0.5) +  # Add dashed trend line
  scale_color_identity() +  # Use the color column for point colors
  labs(x = "Number of sites",
       y = "Number of partners") +
  custom_theme

# Save the plots
ggsave("documents/plots/pelt_occurrence_vs_specif.pdf",
       plot = peltigera_plot,
       width = 8, height = 6, units = "cm")
ggsave("documents/plots/nostoc_occurrence_vs_specif.pdf",
         plot = nostoc_plot,
         width = 8, height = 6, units = "cm")

