#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load the required libraries
library(tidyverse)
# Load some functions
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

# Pair data for pairs that cooccur at least once
pairs_data_cooccur <- pairs_data %>%
    group_by(IDi, IDj) %>%
    filter(sum(Xij) > 10) %>%
    ungroup()

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


##### JOIN E AND SPLIT BY PAIRS #####

# Join pair and climate data
DF <- as.data.frame(pairs_data_cooccur)
DF$E <- as.data.frame(edata)

# Split the DF into a list of pairs of species
pairs_ID <- as.factor(paste(DF$IDi,DF$IDj))
DF_split <- split(DF,pairs_ID)

# Examine specifc pairs in DF_split
#DF_split$`Peltigeraleucophlebia6 XXXIII`[DF_split$`Peltigeraleucophlebia6 XXXIII`$Xij == 1,] 



##### FIT MODELS FOR ALL PAIRS #####

#DF_split <- DF_split[c(806, 913)]
# Specify the environmental variables
Enames = c("MAT"#,
            #"MAT2",
            #"MAP"#,
            #"MAP2"#, 
            #"mean_canopy_closure"#, 
            #"Elevation"#, 
            #"proportion_conifer"
            )

# Fit all models to all pairs
fit_result <- fit_all_pairs(DF_split, Enames)

# Summarize the results across all pairs
full_fit_summary <- fit_result %>%
    group_by(model) %>%
    summarize(LL = sum(LL),
              npars = sum(npars),
              AIC = -2 * LL + 2 * npars)

# Write the full fit summary to a CSV file
write_csv(fit_summary, "documents/tables/model_fit_summary.csv")

# Rank models by AIC - Force ties between C2_L0 and C2_L1 to to score the best model as C2_L0
model_ranks <- fit_result %>%
    group_by(mycobiont, nostoc) %>%
    mutate(rank = min_rank(AIC)) %>%
    ungroup() %>%
    left_join(model_force, by = c("mycobiont" = "IDi", "nostoc" = "IDj"), suffix = c("", "_tied")) %>%
    mutate(rank = if_else(!is.na(coo) & model == "C2_L0", 1, rank),
            rank = if_else(!is.na(coo) & model == "C2_L1", 2, rank)) %>%
    select(-coo, -inter) %>%
    left_join(pairs_data_cooccur_summary, by = c("mycobiont" = "IDi", "nostoc" = "IDj"))

# What are the best models for pairs that never interact?
#model_ranks %>%
#    filter(!observed) %>%
#    group_by(mycobiont, nostoc) %>%
#    summarize(best_model = model[rank == 1]) %>%
#    pull(best_model) %>%
#    unique()

# What fraction of pairs for which the best model is C2_L0 never interact? - 0.9376147
C2_L0_observed <- model_ranks %>%
    filter(model == "C2_L0" & rank == 1 & observed) %>%
    nrow()
C2_L0_unobserved <- model_ranks %>%
    filter(model == "C2_L0" & rank == 1 & !observed) %>%
    nrow()
C2_L0_unobserved / (C2_L0_observed + C2_L0_unobserved)
    

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

# Custom colors for each model
model_colors <- c("C0_L2" = "#8F6D4A", 
                    "C1_L2" = "#CAB79E", 
                    "C2_L0" = "#8A708A", 
                    "C2_L1" = "#C3B8C3", 
                    "C2_L2" = "#7D927A", 
                    "C3_L2" = "#BFC8BC")

# Custom labels for the legend
model_labels <- c("C0_L2" = "Model C0_L2",
                   "C1_L2" = "Model C1_L2",
                   "C2_L0" = "Model C2_L0",
                   "C2_L1" = "Model C2_L1",
                   "C2_L2" = "Model C2_L2",
                   "C3_L2" = "Model C3_L2")

# Model ranks plot including all pairs
model_ranks_plot <- model_ranks %>%
    ggplot(aes(x = rank, fill = model)) +
    geom_bar(position = "fill") +
    labs(x = "Rank (best to worst)", y = "Proportion of pairs") +
    scale_fill_manual(values = model_colors, labels = model_labels) +
    guides(fill = guide_legend(title = "Model")) +
    custom_theme

# Model ranks plot weighted by pair regional frequency
model_ranks_plot_weighted <- model_ranks %>%
    ggplot(aes(x = rank, fill = model, weight = regional_frequency)) +
    geom_bar(position = "fill") +
    labs(x = "Rank (best to worst)", y = "Proportion of pairs weighted by regional freq.") +
    scale_fill_manual(values = model_colors, labels = model_labels) +
    guides(fill = guide_legend(title = "Model")) +
    custom_theme

# Model ranks plot including only pairs with at least one interaction
#model_ranks_plot_interacting <- model_ranks %>%
#    filter(observed) %>%
#    ggplot(aes(x = rank, fill = model)) +
#    geom_bar(position = "fill") +
#    labs(x = "Rank", y = "Proportion") +
#    scale_fill_manual(values = model_colors, labels = model_labels) +
#    guides(fill = guide_legend(title = "Model")) +
#    custom_theme

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

# Save the plots as a PDF
ggsave(model_ranks_plot, filename = "documents/plots/model_ranks_plot.pdf",
    width = 12, height = 9, units = "cm")
ggsave(model_ranks_plot_weighted, filename = "documents/plots/model_ranks_plot_weighted.pdf",
    width = 12, height = 9, units = "cm")
#ggsave(model_ranks_plot_interacting, filename = "documents/plots/model_ranks_plot_interacting.pdf",
#    width = 12, height = 9, units = "cm")
ggsave(assymetry_hist, filename = "documents/plots/assymetry_hist.pdf",
    width = 12, height = 9, units = "cm")
ggsave(cooccurring_sites_hist, filename = "documents/plots/cooccurring_sites_hist.pdf",
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
