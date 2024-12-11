#!/usr/bin/env Rscript

# Load the required libraries
library(tidyverse)


###### MERGE ABMI CYANOLICHEN ID DATASETS - ANALYSIS VERSION ######

# Load ABMI extra dataset
abmi_extra_data <- read_csv("data/tables/abmi_extra_voucher.csv") %>%
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

# Save the ABMI ID dataset
write_csv(abmi_id_data, "data/tables/abmi_id_data.csv")


###### MERGE ABMI CYANOLICHEN ID DATASETS - MANUSCRIPT VERSION ######

# Load ABMI extra dataset
abmi_extra_voucher <- read_csv("data/tables/abmi_extra_voucher.csv") %>%
    # Add ABMI_ prefix to the Site column
    mutate(site = paste0("ABMI_", site),
            abmi_id = as.numeric(abmi_id)) %>%
    # Remove underscores from the mycobiont_molecular_id, subclade, section, and species_complex columns
    mutate(mycobiont_molecular_id = str_replace_all(mycobiont_molecular_id, "_", " "),
           subclade = str_replace_all(subclade, "_", " "),
           section = str_replace_all(section, "_", " "),
           species_complex = str_replace_all(species_complex, "_", " "),
           site = str_replace_all(site, "_", " "),
           Reference = "This study")

# Load ABMI ID dataset from Pardo-De la Hoz et al. (2024)-Data S1C and merge with ABMI extra data
abmi_id_voucher <- read_csv("data/tables/nostoc_datas1c.csv") %>% 
    select(1:Collector) %>%
    # Rename column names as lower case and replace spaces with underscores and add reference column
    rename_all(tolower) %>%
    rename_all(~str_replace_all(., " ", "_")) %>%
    mutate(Reference = "Pardo-De la Hoz et al. (2025)") %>%
    # Merge ABMI ID data with ABMI extra data
    bind_rows(abmi_extra_voucher) %>%
    # Coalesce the columns phylogroup and species complex into a new column called "nostoc_otu"
    mutate(nostoc_otu = coalesce(phylogroup, species_complex)) %>%
    # Rename columns by replacing underscores with spaces and making it sentence case
    rename_all(~str_replace_all(., "_", " ")) %>%
    rename_all(~str_to_title(.)) %>%
    rename(`DNA ID` = `Dna Id`,
            `ABMI ID` = `Abmi Id`,
            `Mycobiont molecular ID` = `Mycobiont Molecular Id`,
            `Nostoc OTU` = `Nostoc Otu`) %>%
    # Reorder columns
    select(`DNA ID`, `ABMI ID`, `Mycobiont molecular ID`, `Nostoc OTU`, Subclade, Section, `Species Complex`, Phylogroup, Site, Year, Collector, Reference)

# Save the ABMI ID voucher
write_csv(abmi_id_voucher, "documents/tables/abmi_id_voucher.csv")


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
    group_by(`ABMI Site`, Year) %>%
  filter(
    case_when(
      any(`Collection Methodology` == 1) ~ `Collection Methodology` == 1,
      TRUE ~ `Collection Methodology` >= 2)) %>%
  ungroup() %>%
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

# Sites included in ABMI data
included_sites <- unique(abmi_id_data$site_year)

# Summarize ABMI site data
abmi_site_data <- d_qha.all %>%
    group_by(SiteYear) %>%
    summarise_all(list(~ if(is.numeric(.)) 
                            mean(., na.rm = TRUE)
                        else 
                            first(.)))  %>%
    # Add dummy SiteYear columns to update the year of data collection for problematic sites
    mutate(SiteYear_dummy_1 = case_when(SiteYear == "1070_2008" ~ "1070_2007",
                                SiteYear == "120_2016" ~ "120_2017",
                                #SiteYear == "1521_2016" ~ "1521_2015",
                                SiteYear == "15_2015" ~ "15_2016",
                                SiteYear == "1_2015" ~ "1_2016",
                                SiteYear == "236_2015" ~ "236_2016",
                                SiteYear == "264_2015" ~ "264_2016",
                                SiteYear == "2_2015" ~ "2_2016",
                                SiteYear == "32_2015" ~ "32_2016",
                                SiteYear == "34_2015" ~ "34_2016",
                                SiteYear == "3_2015" ~ "3_2016",
                                SiteYear == "497_2016" ~ "497_2015",
                                SiteYear == "498_2016" ~ "498_2015",
                                SiteYear == "499_2016" ~ "499_2015",
                                SiteYear == "528_2016" ~ "528_2015",
                                SiteYear == "529_2016" ~ "529_2015",
                                SiteYear == "559_2016" ~ "559_2015",
                                SiteYear == "560_2016" ~ "560_2015",
                                SiteYear == "584_2016" ~ "584_2015",
                                SiteYear == "585_2016" ~ "585_2015",
                                SiteYear == "586_2016" ~ "586_2015",
                                SiteYear == "616_2016" ~ "616_2015",
                                .default = SiteYear)) %>%
    mutate(SiteYear_dummy_2 = case_when(SiteYear == "1070_2008" ~ "1070_2007",
                                SiteYear == "120_2016" ~ "120_2017",
                                #SiteYear == "1521_2016" ~ "1521_2015",
                                SiteYear == "112_2013" ~ "112_2012",
                                SiteYear == "113_2013" ~ "113_2012",
                                SiteYear == "134_2013" ~ "134_2012",
                                SiteYear == "56_2013" ~ "56_2012",
                                SiteYear == "605_2012" ~ "605_2013",
                                SiteYear == "607_2012" ~ "607_2013",
                                SiteYear == "637_2012" ~ "637_2013",
                                SiteYear == "639_2012" ~ "639_2013",
                                SiteYear == "82_2013" ~ "82_2012",
                                SiteYear == "83_2013" ~ "83_2012",
                                .default = SiteYear)) %>%
    left_join(canopy_closure, by = c("SiteYear_dummy_1" = "SiteYear")) %>%
    left_join(ecosite_data, by = c("SiteYear_dummy_2" = "SiteYear")) %>%
    mutate(deciduous_natural = rowSums(select(., 48:56), na.rm = TRUE),
            deciduous_CC = rowSums(select(., 103:107)),
            conifer_natural = rowSums(select(., 66:94)),
            conifer_CC = rowSums(select(., 113:122, 190)),
            mixedwood_natural = rowSums(select(., 57:65)),
            mixedwood_CC = rowSums(select(., 108:112)),
            all_forest_stands = deciduous_natural + deciduous_CC + conifer_natural + conifer_CC + mixedwood_natural + mixedwood_CC,
            proportion_conifer = case_when(all_forest_stands > 0 ~ (conifer_natural + conifer_CC) / all_forest_stands,
                                            .default = 0),
            SiteYear = paste0("ABMI_", SiteYear)) %>%
    select(1:3, 6:13, 16:46, deciduous_natural, deciduous_CC, conifer_natural, 
            conifer_CC, mixedwood_natural, mixedwood_CC, all_forest_stands, proportion_conifer,
            mean_canopy_closure, `Ecosite - Tree Species Modifier`, 
            vegetation_type, `Percent Area of Ecological Site Classification`) %>%
    filter(SiteYear %in% included_sites)

# Assuming ABMI_1521_2016 is grassland, update mean_canopy_closure to 0 and vegetation_type to non_treed
abmi_site_data[abmi_site_data$SiteYear == "ABMI_1521_2016",]$mean_canopy_closure <- 0
abmi_site_data[abmi_site_data$SiteYear == "ABMI_1521_2016",]$vegetation_type <- "non-treed"

# Which records have NA in the mean_canopy_closure column? convert to character before filtering
#abmi_site_data[is.na(as.character(abmi_site_data$mean_canopy_closure)),]$SiteYear

# Which records have NA in the vegetation_type column?
#abmi_site_data[is.na(abmi_site_data$vegetation_type),]$SiteYear

# Which sites are duplicated in the ABMI site data?
#abmi_site_data %>% 
#    group_by(SiteYear) %>%
#    filter(n() > 1)

# Which records have VNA or DNC in the vegetation_type column?
#abmi_site_data[abmi_site_data$vegetation_type == "VNA" | abmi_site_data$vegetation_type == "DNC",]

# Which record has NA in the site column?
#abmi_site_data[is.na(abmi_site_data$Site),]

# Which sites are in included_sites but not in abmi_site_data?
#setdiff(included_sites, abmi_site_data$SiteYear)
#[1] "ABMI_OG-ALPAC-SK-1_2012" "ABMI_OG-ALPAC-SK-2_2012"
#[3] "ABMI_OG-ALPAC-SK-9_2012"

# Save the ABMI site data
write_csv(abmi_site_data, "data/tables/abmi_site_data.csv")