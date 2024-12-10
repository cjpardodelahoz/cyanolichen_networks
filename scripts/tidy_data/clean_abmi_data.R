#!/usr/bin/env Rscript

# Load the required libraries
library(tidyverse)


###### MERGE ABMI CYANOLICHEN ID DATASETS (Pardo-De la Hoz et al + extra) ######




###### SUMMARIZE ABMI SITE DATA ######

# Canopy closure data
canopy_closure <- read_csv("data/tables/publicABMIData12052024/A_T11_Canopy_Closure_AllYear_AB_8348838353460696594.csv") %>%
    mutate(SiteYear = paste(`ABMI Site`, Year, sep = "_")) %>%
    filter(`Canopy Closure (Open dots)` != "VNA", `Canopy Closure (Open dots)` != "DNC") %>%
    mutate(canopy_closure = 96 - as.numeric(`Canopy Closure (Open dots)`)) %>%
    group_by(SiteYear) %>%
    summarise(mean_canopy_closure = mean(canopy_closure, na.rm = TRUE))

# Define vegetation type groups
conifer <- c("(08a) Sb", "(10a) SbLt", "(02c) Sb", "(06c) Spruce", "(04e) Spruce", "(01a) Pine", 
            "(05c) Sb", "(09a) SbLt", "(05b) Spruce", "(04a) Pine", "(02a) Pine", "(03d) Spruce", 
            "(03b) Pine", "(06a) Pine")
deciduous <- c("(04b) PjMix", "(06b) Poplar", "(04c) Aw", "(04d) AwMix", "(05a) Poplar", "(03c) AwMix",
                "(12a) Tree", "(10.5a) Tree")
non_treed <- c("(10c) None", "(11a) None", "(07b) Flood", "(14b) River", "(14a) Lake", "(10b) Shrub", 
                "(07e) Geo", "(09b) Shrub", "(02b) Other", "(07d) Dry", "(07f) Human", "(08b) Shrub", 
                "(03a) None", "(07a) Alpine", "(12b) Shrub", "(07c) Ice", "(10.5b) Shrub", "(10.5c) None",
                "(13a) None")

# Load ecosite data and condense vegetation types
ecosite_data <- read_csv("data/tables/publicABMIData12052024/A_T01C_Site_Capability_AllYear_AB_16861735982408730727.csv") %>%
    filter(`Point Count Station` == "1", Subpoint == "P", `Time Period` == "C") %>%
    mutate(SiteYear = paste(`ABMI Site`, Year, sep = "_"),
        vegetation_type = case_when(`Ecosite - Tree Species Modifier` %in% conifer ~ "conifer",
                                `Ecosite - Tree Species Modifier` %in% deciduous ~ "deciduous",
                                `Ecosite - Tree Species Modifier` %in% non_treed ~ "non-treed", 
                                .default = `Ecosite - Tree Species Modifier`)) %>%
    select(SiteYear, `Ecosite - Tree Species Modifier`, vegetation_type, `Percent Area of Ecological Site Classification`)

# Load ABMI site data
load("data/r_datasets/abmi_north_south_quad.Rdata")
rownames(d_qha.all) <- NULL

# Summarize ABMI site data
abmi_site_data <- d_qha.all %>%
    group_by(SiteYear) %>%
    summarise_all(list(~ if(is.numeric(.)) 
                            mean(., na.rm = TRUE)
                        else 
                            first(.))) %>%
    left_join(canopy_closure, by = "SiteYear") %>%
    left_join(ecosite_data, by = "SiteYear") %>%
    mutate(deciduous_natural = rowSums(select(., 48:56), na.rm = TRUE),
            deciduous_CC = rowSums(select(., 103:107)),
            conifer_natural = rowSums(select(., 66:94)),
            conifer_CC = rowSums(select(., 113:122, 190)),
            mixedwood_natural = rowSums(select(., 57:65)),
            mixedwood_CC = rowSums(select(., 108:112)),
            all_forest_stands = deciduous_natural + deciduous_CC + conifer_natural + conifer_CC + mixedwood_natural + mixedwood_CC,
            proportion_conifer = (conifer_natural + conifer_CC) / all_forest_stands,
            SiteYear = paste0("ABMI_", SiteYear)) %>%
    select(1:3, 6:13, 16:46, deciduous_natural, deciduous_CC, conifer_natural, 
            conifer_CC, mixedwood_natural, mixedwood_CC, all_forest_stands, proportion_conifer,
            mean_canopy_closure, `Ecosite - Tree Species Modifier`, 
            vegetation_type, `Percent Area of Ecological Site Classification`)