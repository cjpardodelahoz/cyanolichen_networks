#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load the required libraries
library(gjam)
library(tidyverse)
library(ggpubr)
library(knitr)
library(cowplot)

##### PREP ABMI CDM #####

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

# Load ABMI ID dataset from Pardo-De la Hoz et al. (2024)-Data S1C
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

# Convert ABMI ID data to Nostoc community data matrix
abmi_nostoc_cdm <- abmi_id_data %>%
    # Group by site_year and nostoc_otu and count the number of occurrences
    group_by(site_year, nostoc_otu) %>%
    summarise(count = n()) %>%
    # Spread the data to create a community matrix
    spread(nostoc_otu, count, fill = 0) %>%
    # Remove the <NA column
    select(-c("<NA>"))

# Convert ABMI ID data to Mycobiont community data matrix
abmi_mycobiont_cdm <- abmi_id_data %>%
    # Group by site_year and mycobiont_molecular_id and count the number of occurrences
    group_by(site_year, mycobiont_molecular_id) %>%
    summarise(count = n()) %>%
    # Spread the data to create a community matrix
    spread(mycobiont_molecular_id, count, fill = 0) %>%
    # Remove the <NA column
    select(-c("<NA>"))

# Are there rows in abmi_nostoc_cdm with all zeros?
empty_sites_nostoc <- abmi_nostoc_cdm %>%
    column_to_rownames("site_year") %>%
    filter_all(all_vars(. == 0)) %>%
    rownames()

# Are there rows in abmi_mycobiont_cdm with all zeros?
empty_sites_mycobiont <- abmi_mycobiont_cdm %>%
    column_to_rownames("site_year") %>%
    filter_all(all_vars(. == 0)) %>%
    rownames()

# Are there sites empty in both the Nostoc and Mycobiont community data matrices?
empty_sites_both <- intersect(empty_sites_nostoc, empty_sites_mycobiont)

# Join the Nostoc and Mycobiont community data matrices by site_year
abmi_cdm <- abmi_nostoc_cdm %>%
    left_join(abmi_mycobiont_cdm, by = "site_year") %>%
    # Remove sites empty in both the Nostoc and Mycobiont community data matrices
    filter(!site_year %in% empty_sites_both) %>%
    # Arrange the data by site_year
    arrange(site_year) %>%
    # Filter Saskatchewan sites
    filter(!str_detect(site_year, "SK")) %>%
    # Make the site_year column the row names
    column_to_rownames("site_year")

# Get list of sites included in the CDM
cdm_sites <- rownames(abmi_cdm)

# Trim the CDM
abmi_cdm_trimmed <- as.data.frame(gjamTrimY(abmi_cdm, minObs = 40)$y)
abmi_cdm_trimmed <- sapply(abmi_cdm_trimmed, as.numeric)
abmi_cdm_trimmed <- as.data.frame(abmi_cdm_trimmed)


##### PREP ABMI SITE DATA #####

load("data/r_datasets/abmi_north_south_quad.Rdata")

rownames(d_qha.all) <- NULL

abmi_site_data <- d_qha.all %>%
    # Add "ABMI_" prefix to the SiteYear column
    mutate(SiteYear = paste0("ABMI_", SiteYear)) %>%
    # Filter QUAD to NE
    filter(QUAD == "NE") %>%
    # Filter to sites included in the CDM
    filter(SiteYear %in% cdm_sites) %>%
    # Ensure SiteYear is ordered the same as in the CDM
    arrange(SiteYear) %>%
    # Make the SiteYear column the row names
    column_to_rownames("SiteYear") %>%
    # Rename DD18 as DDA18 and DD_18 as DDB18
    rename(DDA18 = DD18, DDB18 = DD_18) %>%
    # Remove underscores from column names
    rename_all(~str_replace_all(., "_", ""))

# Which sites are included in the CDM but not in the site data?
#sites_not_in_site_data <- setdiff(cdm_sites, rownames(abmi_site_data))


##### SET UP MODEL AND RUN gJAM #####

# Set up the formula for the model
formula1 <- as.formula(~ Elevation + MWMT + MCMT + MSP)

# Set 𝑁 and 𝑟 for dimension reduction. 𝑁 is the potential number of response groups, and 𝑟 is the dimensionality (AKA flexibility) of those groups. 
# 𝑟 must be smaller than 𝑁 , which should be much smaller than the number of species 𝑆
rlist   <- list(r = 5, N = 8)

# Set up priors

spLo <- "XXXIII"
sp <- length(spLo)
lo <- vector("list", sp)
# add names to the list
names(lo) <- paste0("X_", spLo)
# add values to the list
lo[1] <- Inf
spHi <- c("XXXIII")
sp <- length(spHi)
hi <- vector("list", sp)
# add names to the list
names(hi) <- paste0("X_", spHi)
# add values to the list
hi[1:length(hi)] <- Inf

prior <- gjamPriorTemplate(formula = formula1, 
                           xdata = abmi_site_data, ydata = abmi_cdm_trimmed,
                           lo = lo, hi = hi)

# Effort
#abmi_effort <- list(columns = 1:ncol(abmi_cdm_trimmed), 
#                    values = rep(1, nrow(abmi_cdm_trimmed)))

# Set up the model and run parameters
mlist <- list(ng=40000, burnin=20000, typeNames = 'DA', betaPrior = prior, reductList = rlist)

# Run the model
out <- gjam(formula1, xdata = abmi_site_data, ydata = abmi_cdm_trimmed, modelList = mlist)

# Plot the results
gjamPlot(output = out, plotPars = list(PLOTALLY = T, SAVEPLOTS = T, outfolder = "./"))

out$fit$rmspeAll
out$fit$DIC


# create base for plots
base <- ggplot() +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

base +
  geom_point(aes(x = out$inputs$y[,16],
                 y = out$prediction$ypredMu[,16])) +
  labs(x = "observed", y = "predicted", title = paste0(colnames(out$inputs$y)[16], ", common"))



xMu  <- out$prediction$xpredMu        # inverse prediction of x
xSd  <- out$prediction$xpredSd
yMu  <- out$prediction$ypredMu        # predicted y
hold <- out$modelList$holdoutIndex    # holdout observations (rows)

plot(out$inputs$xUnstand[hold,-1],xMu[hold,-1], cex=.2, xlab='True', ylab='Predictive mean')
title('holdouts in x')
abline( 0, 1, lty=2 )
