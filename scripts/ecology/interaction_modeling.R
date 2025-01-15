#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load the required libraries
library(tidyverse)
library(ggnewscale)

# Load model functions
source("scripts/ecology/gravel2019_functions/species_models.r")
source("scripts/ecology/gravel2019_functions/interactions_models.r")
source("scripts/ecology/gravel2019_functions/get_LL.r")
source("scripts/ecology/gravel2019_functions/get_probs.r")
source("scripts/ecology/gravel2019_functions/fit_models.r")
source("scripts/ecology/gravel2019_functions/fit_all_pairs.r")


##### PREP REGIONAL DATA FOR MODEL FITTING #####

# Load Peltigera-Nostoc interaction matrix
load("analyses/ecology/peltigera_interaction_matrix.RData")
total_interactions <- sum(abmi_matrix_peltigera)

# Load ABMI dataset and filter to include only Peltigera with associated nostoc
abmi_id_data <- read_csv("data/tables/abmi_id_data.csv") %>% 
    filter(!is.na(nostoc_otu) & !is.na(mycobiont_molecular_id) &  str_detect(mycobiont_molecular_id, "Peltigera"))

# Occurrence and interaction combinations in regional data
mycobiont_molecular_id <- abmi_id_data$mycobiont_molecular_id
nostoc_otu <- abmi_id_data$nostoc_otu
site_year <- abmi_id_data$site_year
mycobiont_occurrence <- paste0(mycobiont_molecular_id, site_year) %>% unique()
nostoc_occurrence <- paste0(nostoc_otu, site_year) %>% unique()
interaction_occurrence <- paste0(mycobiont_molecular_id, nostoc_otu, site_year) %>% unique()

# Table with cooccurrence and interaction data for all mycobiont-nostoc pairs
pairs_data <- expand_grid(IDi = unique(abmi_id_data$mycobiont_molecular_id),
                                    IDj = unique(abmi_id_data$nostoc_otu),
                                    site_year = unique(abmi_id_data$site_year)) %>%
    # Remove NAs and Saskatchewan sites
    filter(!is.na(IDi), !is.na(IDj)) %>%
    filter(!str_detect(site_year, "SK")) %>%
    # Score mycobiont occurrence
    mutate(Xi = if_else(paste0(IDi, site_year) %in% mycobiont_occurrence, 1, 0)) %>%
    # Score nostoc occurrence
    mutate(Xj = if_else(paste0(IDj, site_year) %in% nostoc_occurrence, 1, 0)) %>%
    # Score cooccurrence
    mutate(Xij = if_else(Xi == 1 & Xj == 1, 1, 0)) %>%
    # Score interactions
    mutate(Lij = if_else(paste0(IDi, IDj, site_year) %in% interaction_occurrence, 1, 0)) %>%
    # Filter for Peltigera pairs
    filter(str_detect(IDi, "Peltigera")) #%>%
    # Filter for pairs IDi and IDj that cooccur (Xij) at least once
    #group_by(IDi, IDj) %>%
    #filter(sum(Xij) > 0) %>%
    #ungroup()

# Pair data for pairs that cooccur at least 10 times
pairs_data_cooccur <- pairs_data %>%
    group_by(IDi, IDj) %>%
    filter(sum(Xij) > 10) %>%
    ungroup()

# Pairs that cooccur at least 20 times and interact at least once
pairs_cooccur_20 <- pairs_data %>%
    group_by(IDi, IDj) %>%
    filter(sum(Xij) > 20 & sum(Lij) > 1) %>%
    ungroup() %>%
    select(IDi, IDj) %>%
    distinct()

# Pairs that interact every time they cooccur
model_force <- pairs_data_cooccur %>%
    group_by(IDi, IDj) %>%
    summarise(coo = sum(Xij),
                inter = sum(Lij)) %>%
    filter(coo == inter)

# Summarise properties for all pairs
pairs_data_summary <- pairs_data %>%
    group_by(IDi, IDj) %>%
    summarize(observed = sum(Lij) > 0,
            cooccurring_sites = sum(Xij)) %>%
    ungroup() %>%
    mutate(regional_frequency = map2_dbl(IDi, IDj, ~ abmi_matrix_peltigera[.x, .y]) / total_interactions) %>%
    group_by(IDj) %>%
    mutate(nostoc_degree = n_distinct(IDi)) %>%
        mutate(nostoc_degree = sum(observed)) %>%
    ungroup() %>%
    group_by(IDi) %>%
    mutate(mycobiont_degree = n_distinct(IDj)) %>%
        mutate(mycobiont_degree = sum(observed)) %>%
    ungroup() %>%
    mutate(assymetry = abs(nostoc_degree - mycobiont_degree) / (nostoc_degree + mycobiont_degree))

# Summarise properties for each cooccurring pair
pairs_data_cooccur_summary <- pairs_data_cooccur %>%
    group_by(IDi, IDj) %>%
    summarize(observed = sum(Lij) > 0,
            cooccurring_sites = sum(Xij)) %>%
    ungroup() %>%
    mutate(regional_frequency = map2_dbl(IDi, IDj, ~ abmi_matrix_peltigera[.x, .y]) / total_interactions) %>%
    group_by(IDj) %>%
    mutate(nostoc_degree = n_distinct(IDi)) %>%
        mutate(nostoc_degree = sum(observed)) %>%
    ungroup() %>%
    group_by(IDi) %>%
    mutate(mycobiont_degree = n_distinct(IDj)) %>%
        mutate(mycobiont_degree = sum(observed)) %>%
    ungroup() %>%
    mutate(assymetry = abs(nostoc_degree - mycobiont_degree) / (nostoc_degree + mycobiont_degree))


##### PREP REGIONAL SITE DATA #####

# Load ABMI site data
abmi_site_data <- read_csv("data/tables/abmi_site_data.csv")

# Environmental data for included sites
edata <- pairs_data_cooccur %>%
    select(site_year) %>%
    left_join(abmi_site_data, by = c("site_year" = "SiteYear")) %>%
    select(MAT, 
            MAP,
            mean_canopy_closure,
            Elevation,
            proportion_conifer
            ) %>%
    mutate(MAT2 = MAT^2,
        MAP2 = MAP^2)

# are there sites with missing data?
#edata %>%
#    filter(is.na(T) | is.na(PP))


##### JOIN PAIR AND ENVIRONMENTAL DATA AND SPLIT BY PAIRS #####

# Join pair and climate data
DF <- as.data.frame(pairs_data_cooccur)
DF$E <- as.data.frame(edata)

# Split the DF into a list of pairs of species
pairs_ID <- as.factor(paste(DF$IDi,DF$IDj))
DF_split <- split(DF,pairs_ID)

# Examine specifc pairs in DF_split
DF_split$`Peltigeraextenuata1 sppcomplex31a`[DF_split$`Peltigeraextenuata1 sppcomplex31a`$Xi == 1,] 


##### FIT MODELS FOR ALL PAIRS AND SUMMARIZE RESULTS #####

#DF_split <- DF_split[c(806, 913)]

# Specify the environmental variables
#Enames = c(#"MAT"#,
            #"MAT2",
            #"MAP"#,
            #"MAP2"#, 
            #"mean_canopy_closure"#, 
            #"Elevation"#, 
            #"proportion_conifer"
            #)
Enames_MAT <- c("MAT")
Enames_MAP <- c("MAP")
Enames_proportion_conifer <- c("proportion_conifer")
Enames_elevation <- c("Elevation")

# Fit all models to all pairs
#fit_result <- fit_all_pairs(DF_split, Enames)
fit_result_MAT <- fit_all_pairs(DF_split, Enames_MAT)
fit_result_MAP <- fit_all_pairs(DF_split, Enames_MAP)
fit_result_proportion_conifer <- fit_all_pairs(DF_split, Enames_proportion_conifer)
fit_result_elevation <- fit_all_pairs(DF_split, Enames_elevation)

# Summarize the results across all pairs
#full_fit_summary <- fit_result %>%
#    group_by(model) %>%
#    summarize(LL = sum(LL),
#              npars = sum(npars),
#              AIC = -2 * LL + 2 * npars)
full_fit_summary_MAT <- fit_result_MAT %>%
    group_by(model) %>%
    summarize(LL = sum(LL),
              npars = sum(npars),
              AIC = -2 * LL + 2 * npars)
full_fit_summary_MAP <- fit_result_MAP %>%
    group_by(model) %>%
    summarize(LL = sum(LL),
              npars = sum(npars),
              AIC = -2 * LL + 2 * npars)
full_fit_summary_proportion_conifer <- fit_result_proportion_conifer %>%
    group_by(model) %>%
    summarize(LL = sum(LL),
              npars = sum(npars),
              AIC = -2 * LL + 2 * npars)
full_fit_summary_elevation <- fit_result_elevation %>%
    group_by(model) %>%
    summarize(LL = sum(LL),
              npars = sum(npars),
              AIC = -2 * LL + 2 * npars)

# Write the full fit summary to a CSV file
#write_csv(full_fit_summary, "documents/tables/model_fit_summary.csv")
write_csv(full_fit_summary_MAT, "documents/tables/model_fit_summary_MAT.csv")
write_csv(full_fit_summary_MAP, "documents/tables/model_fit_summary_MAP.csv")
write_csv(full_fit_summary_proportion_conifer, "documents/tables/model_fit_summary_proportion_conifer.csv")
write_csv(full_fit_summary_elevation, "documents/tables/model_fit_summary_elevation.csv")

# Rank models by AIC - Force ties between C2_L0 and C2_L1 to to score the best model as C2_L0
#model_ranks <- fit_result %>%
#    group_by(mycobiont, nostoc) %>%
#    mutate(rank = min_rank(AIC)) %>%
#    ungroup() %>%
#    left_join(model_force, by = c("mycobiont" = "IDi", "nostoc" = "IDj"), suffix = c("", "_tied")) %>%
#    mutate(rank = if_else(!is.na(coo) & model == "C2_L0", 1, rank),
#            rank = if_else(!is.na(coo) & model == "C2_L1", 2, rank)) %>%
#    select(-coo, -inter) %>%
#    left_join(pairs_data_cooccur_summary, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))
model_ranks_MAT <- fit_result_MAT %>%
    group_by(mycobiont, nostoc) %>%
    mutate(rank = min_rank(AIC)) %>%
    ungroup() %>%
    left_join(model_force, by = c("mycobiont" = "IDi", "nostoc" = "IDj"), suffix = c("", "_tied")) %>%
    mutate(rank = if_else(!is.na(coo) & model == "C2_L0", 1, rank),
            rank = if_else(!is.na(coo) & model == "C2_L1", 2, rank)) %>%
    select(-coo, -inter) %>%
    left_join(pairs_data_cooccur_summary, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))
model_ranks_MAP <- fit_result_MAP %>%
    group_by(mycobiont, nostoc) %>%
    mutate(rank = min_rank(AIC)) %>%
    ungroup() %>%
    left_join(model_force, by = c("mycobiont" = "IDi", "nostoc" = "IDj"), suffix = c("", "_tied")) %>%
    mutate(rank = if_else(!is.na(coo) & model == "C2_L0", 1, rank),
            rank = if_else(!is.na(coo) & model == "C2_L1", 2, rank)) %>%
    select(-coo, -inter) %>%
    left_join(pairs_data_cooccur_summary, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))
model_ranks_proportion_conifer <- fit_result_proportion_conifer %>%
    group_by(mycobiont, nostoc) %>%
    mutate(rank = min_rank(AIC)) %>%
    ungroup() %>%
    left_join(model_force, by = c("mycobiont" = "IDi", "nostoc" = "IDj"), suffix = c("", "_tied")) %>%
    mutate(rank = if_else(!is.na(coo) & model == "C2_L0", 1, rank),
            rank = if_else(!is.na(coo) & model == "C2_L1", 2, rank)) %>%
    select(-coo, -inter) %>%
    left_join(pairs_data_cooccur_summary, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))
model_ranks_elevation <- fit_result_elevation %>%
    group_by(mycobiont, nostoc) %>%
    mutate(rank = min_rank(AIC)) %>%
    ungroup() %>%
    left_join(model_force, by = c("mycobiont" = "IDi", "nostoc" = "IDj"), suffix = c("", "_tied")) %>%
    mutate(rank = if_else(!is.na(coo) & model == "C2_L0", 1, rank),
            rank = if_else(!is.na(coo) & model == "C2_L1", 2, rank)) %>%
    select(-coo, -inter) %>%
    left_join(pairs_data_cooccur_summary, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))

# Load module assignments
load("analyses/ecology/peltigera_module_assignments.RData")
peltigera_module_assignments <- mutate(peltigera_module_assignments, 
    module = factor(as.character(module)))
nostoc_module_assignments <- mutate(nostoc_module_assignments, 
    module = factor(as.character(module)))

# Table for delta AIC plot comparing C2_L1 to C2_L2
#delta_aic_table <- model_ranks %>%
#    filter(model %in% c("C2_L1", "C2_L2")) %>%
#    group_by(mycobiont, nostoc) %>%
#    summarize(delta_aic = AIC[model == "C2_L2"] - AIC[model == "C2_L1"]) %>%
#    ungroup() %>%
#    right_join(pairs_cooccur_20, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))
delta_aic_table_MAT <- model_ranks_MAT %>%
    filter(model %in% c("C2_L1", "C2_L2")) %>%
    group_by(mycobiont, nostoc) %>%
    summarize(delta_aic = AIC[model == "C2_L2"] - AIC[model == "C2_L1"]) %>%
    ungroup() %>%
    right_join(pairs_cooccur_20, by = c("mycobiont" = "IDi", "nostoc" = "IDj")) %>%
    left_join(peltigera_module_assignments, by = c("mycobiont" = "mycobiont_molecular_id")) %>%
    left_join(nostoc_module_assignments, by = c("nostoc" = "nostoc_otu"))
delta_aic_table_MAP <- model_ranks_MAP %>%
    filter(model %in% c("C2_L1", "C2_L2")) %>%
    group_by(mycobiont, nostoc) %>%
    summarize(delta_aic = AIC[model == "C2_L2"] - AIC[model == "C2_L1"]) %>%
    ungroup() %>%
    right_join(pairs_cooccur_20, by = c("mycobiont" = "IDi", "nostoc" = "IDj")) %>%
    left_join(peltigera_module_assignments, by = c("mycobiont" = "mycobiont_molecular_id")) %>%
    left_join(nostoc_module_assignments, by = c("nostoc" = "nostoc_otu"))
delta_aic_table_proportion_conifer <- model_ranks_proportion_conifer %>%
    filter(model %in% c("C2_L1", "C2_L2")) %>%
    group_by(mycobiont, nostoc) %>%
    summarize(delta_aic = AIC[model == "C2_L2"] - AIC[model == "C2_L1"]) %>%
    ungroup() %>%
    right_join(pairs_cooccur_20, by = c("mycobiont" = "IDi", "nostoc" = "IDj")) %>%
    left_join(peltigera_module_assignments, by = c("mycobiont" = "mycobiont_molecular_id")) %>%
    left_join(nostoc_module_assignments, by = c("nostoc" = "nostoc_otu"))
delta_aic_table_elevation <- model_ranks_elevation %>%
    filter(model %in% c("C2_L1", "C2_L2")) %>%
    group_by(mycobiont, nostoc) %>%
    summarize(delta_aic = AIC[model == "C2_L2"] - AIC[model == "C2_L1"]) %>%
    ungroup() %>%
    right_join(pairs_cooccur_20, by = c("mycobiont" = "IDi", "nostoc" = "IDj")) %>%
    left_join(peltigera_module_assignments, by = c("mycobiont" = "mycobiont_molecular_id")) %>%
    left_join(nostoc_module_assignments, by = c("nostoc" = "nostoc_otu"))

# What fraction of pairs for which the best model is C2_L0 never interact? - 0.9069767
#C2_L0_observed <- model_ranks %>%
#    filter(model == "C2_L0" & rank == 1 & observed) %>%
#    nrow()
#C2_L0_unobserved <- model_ranks %>%
#    filter(model == "C2_L0" & rank == 1 & !observed) %>%
#    nrow()
#C2_L0_unobserved / (C2_L0_observed + C2_L0_unobserved)

# What are the best models for pairs that never interact?
#model_ranks %>%
#    filter(!observed) %>%
#    group_by(mycobiont, nostoc) %>%
#    summarize(best_model = model[rank == 1]) %>%
#    pull(best_model) %>%
#    unique()

# Calculate the percent likelihood gain for the best model compared to the second-best model
#likelihood_gain <- fit_result %>%
#    group_by(mycobiont, nostoc) %>%
#    mutate(rank = rank(AIC)) %>%
#    filter(rank %in% c(1, 2 ,6)) %>%
#    group_by(mycobiont, nostoc) %>%
#    summarize(percent_likelihood_gain = -(LL[rank == 2] - LL[rank == 1]) / (LL[rank == 6] - LL[rank == 1])*100)
    

##### PLOT MODEL FITTING RESULTS #####

# Define a custom theme for plots
custom_theme <- theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), 
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 12)
    )
custom_theme_1 <- theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(color = "black", size = 7),
    legend.text = element_text(color = "black", size = 7),
    legend.title = element_text(color = "black", size = 7)
    )

# Define custom colors for the modules
module_colors <- c("1" = "#8a708a", "2" = "#7e937b", "3" = "#3277b0", "4" = "#be8551", "6" = "gray20")

# Delta AIC plots comparing C2_L1 to C2_L2

# Define the order of pairs based on delta_aic_plot_MAT
pair_order <- delta_aic_table_MAT %>%
    arrange(desc(mycobiont), desc(nostoc)) %>%
    mutate(pair = paste(mycobiont, nostoc)) %>%
    pull(pair)

# Add the pair column to all tables
delta_aic_table_MAT <- delta_aic_table_MAT %>%
  mutate(pair = factor(paste(mycobiont, nostoc), levels = pair_order))

delta_aic_table_MAP <- delta_aic_table_MAP %>%
  mutate(pair = factor(paste(mycobiont, nostoc), levels = pair_order))

delta_aic_table_proportion_conifer <- delta_aic_table_proportion_conifer %>%
  mutate(pair = factor(paste(mycobiont, nostoc), levels = pair_order))

delta_aic_table_elevation <- delta_aic_table_elevation %>%
  mutate(pair = factor(paste(mycobiont, nostoc), levels = pair_order))

# Plot for MAT
delta_aic_plot_MAT <- delta_aic_table_MAT %>%
  ggplot(aes(x = pair, y = delta_aic, fill = delta_aic > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#c7d19f", "FALSE" = "#aaa58d"), guide = "none") +
  labs(x = "Pair", y = "Delta AIC (C2_L2 - C2_L1)") +
  new_scale_color() +
  geom_point(data = delta_aic_table_MAT,
             aes(y = -10, x = pair, color = module.x), alpha = 0.85, shape = 17, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  new_scale_color() +
  geom_point(data = delta_aic_table_MAT,
             aes(y = -9.5, x = pair, color = module.y), alpha = 0.85, shape = 16, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_y_continuous(limits = c(-10.5, 4), breaks = seq(-8, 4, 2)) +
  coord_flip() +
  custom_theme_1

# Plot for MAP
delta_aic_plot_MAP <- delta_aic_table_MAP %>%
  ggplot(aes(x = pair, y = delta_aic, fill = delta_aic > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#c7d19f", "FALSE" = "#aaa58d"), guide = "none") +
  labs(x = "Pair", y = "Delta AIC (C2_L2 - C2_L1)") +
  new_scale_color() +
  geom_point(data = delta_aic_table_MAP,
             aes(y = -4, x = pair, color = module.x), alpha = 0.85, shape = 17, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  new_scale_color() +
  geom_point(data = delta_aic_table_MAP,
             aes(y = -3.5, x = pair, color = module.y), alpha = 0.85, shape = 16, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-2, 4, 2)) +
  coord_flip() +
  custom_theme_1

# Plot for Proportion Conifer
delta_aic_plot_proportion_conifer <- delta_aic_table_proportion_conifer %>%
  ggplot(aes(x = pair, y = delta_aic, fill = delta_aic > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#c7d19f", "FALSE" = "#aaa58d"), guide = "none") +
  labs(x = "Pair", y = "Delta AIC (C2_L2 - C2_L1)") +
  new_scale_color() +
  geom_point(data = delta_aic_table_proportion_conifer,
             aes(y = -4, x = pair, color = module.x), alpha = 0.85, shape = 17, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  new_scale_color() +
  geom_point(data = delta_aic_table_proportion_conifer,
             aes(y = -3.5, x = pair, color = module.y), alpha = 0.85, shape = 16, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-2, 4, 2)) +
  coord_flip() +
  custom_theme_1

# Plot for Elevation
delta_aic_plot_elevation <- delta_aic_table_elevation %>%
  ggplot(aes(x = pair, y = delta_aic, fill = delta_aic > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#c7d19f", "FALSE" = "#aaa58d"), guide = "none") +
  labs(x = "Pair", y = "Delta AIC (C2_L2 - C2_L1)") +
  new_scale_color() +
  geom_point(data = delta_aic_table_elevation,
             aes(y = -4, x = pair, color = module.x), alpha = 0.85, shape = 17, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  new_scale_color() +
  geom_point(data = delta_aic_table_elevation,
             aes(y = -3.5, x = pair, color = module.y), alpha = 0.85, shape = 16, size = 1.9) +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-2, 4, 2)) +
  coord_flip() +
  custom_theme_1

# Histogram with specialization assymetry for top two models
assymetry_hist <- model_ranks %>%
    filter(model %in% c("C2_L0", "C2_L1") & rank == 1 & observed) %>%
    mutate(model = factor(model, levels = c("C2_L1", "C2_L0"))) %>%
    ggplot(aes(x = assymetry, fill = model)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 15) +
    labs(x = "Specialization assymetry", y = "No. of pairs") +
    scale_fill_manual(
        values = c("C2_L0" = "#8A708A", "C2_L1" = "#C3B8C3"),
        labels = c("C2_L0" = "Model C2_L0", "C2_L1" = "Model C2_L1")
    ) +
    guides(fill = guide_legend(title = "Best model")) +
    geom_vline(data = model_ranks %>%
                       filter(model %in% c("C2_L0", "C2_L1") & rank == 1 & observed) %>%
                       group_by(model) %>%
                       summarize(median_assymetry = median(assymetry)),
                   aes(xintercept = median_assymetry, color = model),
                   linetype = "dashed", size = 0.75) +
        scale_color_manual(
            values = c("C2_L0" = "#8A708A", "C2_L1" = "#C3B8C3"),
            labels = c("C2_L0" = "Model C2_L0", "C2_L1" = "Model C2_L1")
        ) +
    guides(color = "none") +
    custom_theme

# Most pairs that do not interact never cooccur
cooccurring_sites_hist <- pairs_data_summary %>%
    ggplot(aes(x = cooccurring_sites, fill = observed)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 45) +
    labs(x = "Number of cooccurring sites", y = "No. of pairs") +
    scale_fill_manual(values = c("TRUE" = "gray25", "FALSE" = "gray70"), labels = c("TRUE" = "Observed", "FALSE" = "Not observed")) +
    guides(fill = guide_legend(title = "Interaction")) +
    custom_theme

# Custom colors for each model
#model_colors <- c("C0_L2" = "#8F6D4A", 
#                    "C1_L2" = "#CAB79E", 
#                    "C2_L0" = "#8A708A", 
#                    "C2_L1" = "#C3B8C3", 
#                    "C2_L2" = "#7D927A", 
#                    "C3_L2" = "#BFC8BC")

# Custom labels for the legend
#model_labels <- c("C0_L2" = "Model C0_L2",
#                   "C1_L2" = "Model C1_L2",
#                   "C2_L0" = "Model C2_L0",
#                   "C2_L1" = "Model C2_L1",
#                   "C2_L2" = "Model C2_L2",
#                   "C3_L2" = "Model C3_L2")

# Model ranks plot including all pairs
#model_ranks_plot <- model_ranks %>%
#    ggplot(aes(x = rank, fill = model)) +
#    geom_bar(position = "fill") +
#    labs(x = "Rank (best to worst)", y = "Proportion of pairs") +
#    scale_fill_manual(values = model_colors, labels = model_labels) +
#    guides(fill = guide_legend(title = "Model")) +
#    custom_theme

# Model ranks plot weighted by pair regional frequency
#model_ranks_plot_weighted <- model_ranks %>%
#    ggplot(aes(x = rank, fill = model, weight = regional_frequency)) +
#    geom_bar(position = "fill") +
#    labs(x = "Rank (best to worst)", y = "Proportion of pairs weighted by regional freq.") +
#    scale_fill_manual(values = model_colors, labels = model_labels) +
#    guides(fill = guide_legend(title = "Model")) +
#    custom_theme

# Model ranks plot including only pairs with at least one interaction
#model_ranks_plot_interacting <- model_ranks %>%
#    filter(observed) %>%
#    ggplot(aes(x = rank, fill = model)) +
#    geom_bar(position = "fill") +
#    labs(x = "Rank", y = "Proportion") +
#    scale_fill_manual(values = model_colors, labels = model_labels) +
#    guides(fill = guide_legend(title = "Model")) +
#    custom_theme

# Save the plots as a PDF
#ggsave(model_ranks_plot, filename = "documents/plots/model_ranks_plot.pdf",
#    width = 12, height = 9, units = "cm")
#ggsave(model_ranks_plot_weighted, filename = "documents/plots/model_ranks_plot_weighted.pdf",
#    width = 12, height = 9, units = "cm")
ggsave(delta_aic_plot_MAT, filename = "documents/plots/delta_aic_mat.pdf",
    width = 12, height = 12, units = "cm")
ggsave(delta_aic_plot_MAP, filename = "documents/plots/delta_aic_map.pdf",
    width = 12, height = 12, units = "cm")
ggsave(delta_aic_plot_proportion_conifer, filename = "documents/plots/delta_aic_proportion_conifer.pdf",
    width = 12, height = 12, units = "cm")
ggsave(delta_aic_plot_elevation, filename = "documents/plots/delta_aic_elevation.pdf",
    width = 12, height = 12, units = "cm")
ggsave(assymetry_hist, filename = "documents/plots/assymetry_hist.pdf",
    width = 12, height = 9, units = "cm")
ggsave(cooccurring_sites_hist, filename = "documents/plots/cooccurring_sites_hist.pdf",
    width = 12, height = 9, units = "cm")


##### PLOT COOCCURRENCE VS INTERACTIONS #####

# Pair data for pairs that cooccur at least 10 times
coo_vs_inter_data <- pairs_data %>%
    group_by(IDi, IDj) %>%
    summarise(coo = sum(Xij),
                inter = sum(Lij)) %>%
    left_join(peltigera_module_assignments, by = c("IDi" = "mycobiont_molecular_id")) %>%
    left_join(nostoc_module_assignments, by = c("IDj" = "nostoc_otu")) %>%
    mutate(same_module = module.x == module.y)

# Make a plot of inter vs coo with geom point witht the shape 18 and color by same_module (use gray50 if TRUE and #892b46 if FALSE) with alpha 0.5
coo_vs_inter_plot <- coo_vs_inter_data %>%
    ggplot(aes(x = coo, y = inter, color = same_module)) +
    geom_point(shape = 16, alpha = 0.6, size = 2) +
    scale_color_manual(values = c("TRUE" = "#892b46", "FALSE" = "gray50"), guide = "none") +
    labs(x = "No. of sites where symbionts cooccur", y = "No. of sites where symbionts interact") +
    # Add a one to one line
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 0.35) +
    # Add a horizontal line at y = 0.99
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", linewidth = 0.35) +
    custom_theme

# Save the plot as a PDF
ggsave(coo_vs_inter_plot, filename = "documents/plots/coo_vs_inter_plot.pdf",
    width = 12, height = 9, units = "cm")

##### EXAMPLE PAIR #####

# Find the nb of links per pair
nL  = numeric(length(DF_split))
nXi = numeric(length(DF_split))
nXj = numeric(length(DF_split))
nXij = numeric(length(DF_split))
for(x in 1:length(DF_split)) {
	nL[x] = sum(DF_split[[x]]$Lij)
	nXi[x] = sum(DF_split[[x]]$Xi)
	nXj[x] = sum(DF_split[[x]]$Xj)
	nXij[x] = sum(DF_split[[x]]$Xij)
}

cbind(nL,nXi,nXj,nXij,nL,nL/nXij)[nL/nXij<0.7 & nL > 10,]

# Subset the data
pair_index = which(nL == 36 & nXij == 54) 
data = DF_split[[pair_index]]
#data$E = data.frame(T = data$E$T/12, PP = data$E$PP/1000, T2 = data$E$T2/12^2, PP2 = data$E$PP2/1000^2)
sum(data$Lij)
sum(data$Xij)

# Specify the environmental variables
Enames = c("MAT")

# Pick the model
models_C0_L1 = fit_models(data, selection = FALSE, funC = C0, funL = L1, Enames)
models_C2_L0 = fit_models(data, selection = FALSE, funC = C2, funL = L0, Enames)
models_C2_L1 = fit_models(data, selection = FALSE, funC = C2, funL = L1, Enames)
models_C2_L2 = fit_models(data, selection = FALSE, funC = C2, funL = L2, Enames)
models_C0_L2 = fit_models(data, selection = FALSE, funC = C0, funL = L2, Enames)
models_C1_L2 = fit_models(data, selection = FALSE, funC = C1, funL = L2, Enames)
models_C3_L2 = fit_models(data, selection = FALSE, funC = C3, funL = L2, Enames)

# Compute the LL
LL_C0_L1 = get_LL(models_C0_L1, data)
LL_C2_L0 = get_LL(models_C2_L0, data)
LL_C2_L1 = get_LL(models_C2_L1, data)
LL_C2_L2 = get_LL(models_C2_L2, data)
LL_C0_L2 = get_LL(models_C0_L2, data)
LL_C1_L2 = get_LL(models_C1_L2, data)
LL_C3_L2 = get_LL(models_C3_L2, data)

# Collect the results
LL = c(
    LL_C0_L1[1],
	LL_C2_L0[1],
	LL_C2_L1[1],
	LL_C2_L2[1],
	LL_C0_L2[1],
	LL_C1_L2[1],
	LL_C3_L2[1]	
	)

npars = c(
    LL_C0_L1[2],
	LL_C2_L0[2],
	LL_C2_L1[2],
	LL_C2_L2[2],
	LL_C0_L2[2],
	LL_C1_L2[2],
	LL_C3_L2[2]	
	)

AIC = -2*LL + 2*npars

cbind(LL,npars,AIC)


##### DUMMY DATASET FOR FORMATTING TROUBLESHOOTING #####

# Get rows 7 to 9 of the ABMI ID data
dummy <- abmi_id_data[7:9, ]

# Mycobiont molecular ID and site_year values in dummy
mycobiont_molecular_id <- dummy$mycobiont_molecular_id
nostoc_otu <- dummy$nostoc_otu
site_year <- dummy$site_year
mycobiont_occurrence <- paste0(mycobiont_molecular_id, site_year) %>% unique()
nostoc_occurrence <- paste0(nostoc_otu, site_year) %>% unique()
interaction_occurrence <- paste0(mycobiont_molecular_id, nostoc_otu, site_year) %>% unique()

# Create a table where each row has on of all possible unique combination of mycobiont_molecular_id, nostoc_otu, and site_year values in dummy
dummy_combinations <- expand_grid(IDi = unique(dummy$mycobiont_molecular_id),
                                  IDj = unique(dummy$nostoc_otu),
                                  site_year = unique(dummy$site_year)) %>%
    # Remove NAs from the table
    filter(!is.na(IDi), !is.na(IDj)) %>%
    # Add a column called "Xi" with a value of 1 if the combination of mycobiont_molecular_id and site_year is present in the dummy data and 0 otherwise
    mutate(Xi = if_else(paste0(IDi, site_year) %in% mycobiont_occurrence, 1, 0)) %>%
    # Add a column called "Xj" with a value of 1 if the combination of nostoc_otu and site_year is present in the dummy data and 0 otherwise
    mutate(Xj = if_else(paste0(IDj, site_year) %in% nostoc_occurrence, 1, 0)) %>%
    # Add a column called "Xij" with value 1 if both Xi and Xj are 1 and 0 otherwise
    mutate(Xij = if_else(Xi == 1 & Xj == 1, 1, 0)) %>%
    # Add a column called "Lij" with a value of 1 if the combination of mycobiont_molecular_id, nostoc_otu, and site_year is present in the dummy data and 0 otherwise
    mutate(Lij = if_else(paste0(IDi, IDj, site_year) %in% interaction_occurrence, 1, 0))
