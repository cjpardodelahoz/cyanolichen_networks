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


##### PREP REGIONAL DATA FOR MODEL FITTING #####

# Load ABMI extra dataset
abmi_extra_data <- read_csv("documents/tables/abmi_extra_voucher.csv") %>%
    # Add ABMI_ prefix to the Site column
    mutate(site = paste0("ABMI_", site),
            abmi_id = as.numeric(abmi_id)) %>%
    # Remove underscores and periods from the mycobiont_molecular_id, subclade, section, and species_complex columns
    mutate(mycobiont_molecular_id = str_replace_all(mycobiont_molecular_id, "_", "")) %>%
    mutate(mycobiont_molecular_id = str_replace_all(mycobiont_molecular_id, "\\.", "")) %>%
    mutate(subclade = str_replace_all(subclade, "_", "")) %>%
    mutate(section = str_replace_all(section, "_", "")) %>%
    mutate(section = str_replace_all(section, "\\.", "")) %>%
    mutate(species_complex = str_replace_all(species_complex, "_", "")) %>%
    mutate(species_complex = str_replace_all(species_complex, "\\.", ""))

# Load ABMI ID dataset from Pardo-De la Hoz et al. (2024)-Data S1C and merge with ABMI extra data
abmi_id_data <- read_csv("data/tables/nostoc_datas1c.csv") %>% 
    select(1:Collector) %>%
    # Rename column names as lower case and replace spaces with underscores
    rename_all(tolower) %>%
    rename_all(~str_replace_all(., " ", "_")) %>%
    # Merge ABMI ID data with ABMI extra data
    bind_rows(abmi_extra_data) %>%
    # Replace spaces and remove periods from Mycobiont_molecular_ID column with underscores
    mutate(mycobiont_molecular_id = str_replace_all(mycobiont_molecular_id, " ", "")) %>%
    mutate(mycobiont_molecular_id = str_replace_all(mycobiont_molecular_id, "\\.", "")) %>%
    # Remove spaces and periods from section and species_complex columns with underscores
    mutate(section = str_replace_all(section, " ", "")) %>%
    mutate(section = str_replace_all(section, "\\.", "")) %>%
    mutate(species_complex = str_replace_all(species_complex, " ", "")) %>%
    mutate(species_complex = str_replace_all(species_complex, "\\.", "")) %>%
    # Coalesce the columns phylogroup and species complex into a new column called "nostoc_otu"
    mutate(nostoc_otu = coalesce(phylogroup, species_complex)) %>%
    # Generate a site_year column by concatenating the site and year columns
    mutate(site_year = paste(str_replace(site, " ", "_"), year, sep = "_"))

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
    mutate(Lij = if_else(paste0(IDi, IDj, site_year) %in% interaction_occurrence, 1, 0))


##### PREP REGIONAL SITE DATA #####

# Load ABMI site data
load("data/r_datasets/abmi_north_south_quad.Rdata")
rownames(d_qha.all) <- NULL

# Format site data to match pairs data
abmi_site_data <- d_qha.all %>%
    # Add "ABMI_" prefix to the SiteYear column
    mutate(SiteYear = paste0("ABMI_", SiteYear)) %>%
    # Filter QUAD to NE
    filter(QUAD == "NE") %>%
    # Rename DD18 as DDA18 and DD_18 as DDB18
    rename(DDA18 = DD18, DDB18 = DD_18) %>%
    # Remove underscores from column names
    rename_all(~str_replace_all(., "_", ""))

# Environmental data for included sites
edata <- pairs_data %>%
    select(site_year) %>%
    left_join(abmi_site_data, by = c("site_year" = "SiteYear")) %>%
    select(MAT, MAP) %>%
    mutate(T2 = MAT^2,
        PP2 = MAP^2) %>%
    rename(T = MAT,
        PP = MAP)

# are there sites with missing data?
#edata %>%
#    filter(is.na(MAT) | is.na(MAP))


##### JOIN E AND SPLIT BY PAIRS #####

# Join pair and climate data
DF <- as.data.frame(pairs_data)
DF$E <- as.data.frame(edata)

# Split the DF into a list of pairs of species
pairs_ID <- as.factor(paste(DF$IDi,DF$IDj))
DF_split <- split(DF,pairs_ID)


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
pair_index = which(nL == 25 & nXij == 43) 
data = DF_split[[pair_index]]
data$E = data.frame(T = data$E$T/12, PP = data$E$PP/1000, T2 = data$E$T2/12^2, PP2 = data$E$PP2/1000^2)
sum(data$Lij)
sum(data$Xij)

# Pick the model
models_C2_L0 = fit_models(data, selection = FALSE, funC = C2, funL = L0)
models_C2_L1 = fit_models(data, selection = FALSE, funC = C2, funL = L1)
models_C2_L2 = fit_models(data, selection = FALSE, funC = C2, funL = L2)
models_C0_L2 = fit_models(data, selection = FALSE, funC = C0, funL = L2)
models_C1_L2 = fit_models(data, selection = FALSE, funC = C1, funL = L2)
models_C3_L2 = fit_models(data, selection = FALSE, funC = C3, funL = L2)

# Compute the LL
LL_C2_L0 = get_LL(models_C2_L0, data)
LL_C2_L1 = get_LL(models_C2_L1, data)
LL_C2_L2 = get_LL(models_C2_L2, data)
LL_C0_L2 = get_LL(models_C0_L2, data)
LL_C1_L2 = get_LL(models_C1_L2, data)
LL_C3_L2 = get_LL(models_C3_L2, data)

# Collect the results
LL = c(
	LL_C2_L0[1],
	LL_C2_L1[1],
	LL_C2_L2[1],
	LL_C0_L2[1],
	LL_C1_L2[1],
	LL_C3_L2[1]	
	)

npars = c(
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
