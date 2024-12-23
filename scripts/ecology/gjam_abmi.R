#!/usr/bin/env Rscript

##### LOAD REQUIRED LIBRARIES AND CUSTOM FUNCTIONS #####

# Load the required libraries
library(gjam)
library(tidyverse)
library(ggpubr)
library(knitr)
library(cowplot)
library(ggnewscale)

# Function to plot the environmental responses
plot_env_response <- function(gjam_beta, env_variable, symb, xlab, ylab) {
    # Define a custom theme for plots
    custom_theme <- theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), 
        axis.text = element_text(color = "black", size = 12)
    )
    
    # Create the boxplot
    gjam_beta %>%
        filter(variable == env_variable & symbiont == symb) %>%
        ggplot(aes(x = species, y = Estimate)) +
        geom_pointrange(aes(ymin = CI_025, ymax = CI_975, color = sig95)) +
        scale_color_manual(values = c("YES" = "black", "NO" = "gray50")) +
        guides(color = "none") +
        new_scale_color() +
        geom_errorbar(aes(ymin = CI_025, ymax = CI_975, color = sig95), width = 0.2) +
        scale_color_manual(values = c("YES" = "black", "NO" = "gray50")) +
        scale_y_continuous(limits = c(-0.75, 1.35)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        labs(x = xlab, y = ylab) +
        custom_theme +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        guides(color = "none")
}


##### PREP ABMI CDM #####

# Load ABMI ID data
abmi_id_data <- read_csv("data/tables/abmi_id_data.csv")

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
    column_to_rownames("site_year") #%>%
    #mutate_all(~ if_else(. > 0, 1, 0))

# Get list of sites included in the CDM
cdm_sites <- rownames(abmi_cdm)

# Trim the CDM
abmi_cdm_trimmed <- as.data.frame(gjamTrimY(abmi_cdm, minObs = 30)$y)
abmi_cdm_trimmed <- sapply(abmi_cdm_trimmed, as.numeric)
abmi_cdm_trimmed <- as.data.frame(abmi_cdm_trimmed)


##### PREP ABMI SITE DATA #####

# Load ABMI site data and format for GJAM
abmi_site_data <- read_csv("data/tables/abmi_site_data.csv") %>%
    # Convert to missing for vegetation type
    mutate(vegetation_type = as.factor(if_else(
                str_detect(vegetation_type, "deciduous|non-treed|conifer"), 
                    vegetation_type, NA))) %>%
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
    

##### SET UP MODEL AND RUN GJAM #####

# Set up the formula for the model
formula1 <- as.formula(~ Elevation + MAT + MAP + proportionconifer)

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

# Set up the model and run parameters
mlist <- list(ng = 80000, burnin = 40000, typeNames = 'DA', betaPrior = prior, reductList = rlist)

# Run 3 iterations of the model
set.seed(10403)
gjam1 <- gjam(formula1, xdata = abmi_site_data, ydata = abmi_cdm_trimmed, modelList = mlist)
set.seed(63490)
gjam2 <- gjam(formula1, xdata = abmi_site_data, ydata = abmi_cdm_trimmed, modelList = mlist)
set.seed(78477)
gjam3 <- gjam(formula1, xdata = abmi_site_data, ydata = abmi_cdm_trimmed, modelList = mlist)

# Save the GJAM outputs
save(gjam1, gjam2, gjam3, file = "analyses/ecology/gjam_output.RData")

# Plot the results
gjamPlot(output = gjam1, plotPars = list(PLOTALLY = T, SAVEPLOTS = T, outfolder = "./")) # FOR SOME REASON THIS DOESN'T LIKE SAVING TO OTHER DIRECTORIES


##### PREPARE BETAS FOR PLOTTING #####

# Load the GJAM outputs
#load("analyses/ecology/gjam_output.RData")

# Load module assignments
load("analyses/ecology/peltigera_module_assignments.RData")

# Merge and reorder module assigments for both symbionts
peltigera_module_assignments <- peltigera_module_assignments %>%
    mutate(species = mycobiont_molecular_id,
           symbiont = "peltigera") %>%
    select(species, symbiont, module)
nostoc_module_assignments <- nostoc_module_assignments %>%
    mutate(species = nostoc_otu,
           symbiont = "nostoc") %>%
    select(species, symbiont, module)
module_assignments <- bind_rows(peltigera_module_assignments, nostoc_module_assignments)

# Get beta coefficients from the GJAM output
module_order <- c(1, 5, 2, 4, 3, 6)
gjam_beta <- gjam1$parameters$betaStandXWTable %>%
    rownames_to_column("key") %>%
    separate(key, into = c("species", "variable"), sep = "_") %>%
    left_join(module_assignments, by = "species") %>%
    mutate(sig95 = if_else(sig95 == "*", "YES", "NO"),
        module = factor(module, levels = module_order)) %>%
        arrange(module, species) %>%
        mutate(species = factor(species, levels = unique(species)))


##### PLOT ENVIRONMENTAL RESPONSES #####

# Plot the MAT response
mat_response_pelt <- plot_env_response(gjam_beta, "MAT", "peltigera", "Peltigera species", "Mean annual temperature effect")
mat_response_nostoc <- plot_env_response(gjam_beta, "MAT", "nostoc", "Nostoc species", "Mean annual temperature effect")

# Plot the MAP response
map_response_pelt <- plot_env_response(gjam_beta, "MAP", "peltigera", "Peltigera species", "Mean annual precipitation effect")
map_response_nostoc <- plot_env_response(gjam_beta, "MAP", "nostoc", "Nostoc species", "Mean annual precipitation effect")

# Plot the elevation response
elevation_response_pelt <- plot_env_response(gjam_beta, "Elevation", "peltigera", "Peltigera species", "Elevation effect")
elevation_response_nostoc <- plot_env_response(gjam_beta, "Elevation", "nostoc", "Nostoc species", "Elevation effect")

# Plot the proportion conifer response
conifer_response_pelt <- plot_env_response(gjam_beta, "proportionconifer", "peltigera", "Peltigera species", "Proportion conifer effect")
conifer_response_nostoc <- plot_env_response(gjam_beta, "proportionconifer", "nostoc", "Nostoc species", "Proportion conifer effect")

# Save the plots
ggsave(mat_response_pelt, filename = "documents/plots/mat_response_pelt.pdf",
    width = 10, height = 12, units = "cm")
ggsave(mat_response_nostoc, filename = "documents/plots/mat_response_nostoc.pdf",
    width = 10, height = 12, units = "cm")
ggsave(map_response_pelt, filename = "documents/plots/map_response_pelt.pdf",
    width = 10, height = 12, units = "cm")
ggsave(map_response_nostoc, filename = "documents/plots/map_response_nostoc.pdf",
    width = 10, height = 12, units = "cm")
ggsave(elevation_response_pelt, filename = "documents/plots/elevation_response_pelt.pdf",
    width = 10, height = 12, units = "cm")
ggsave(elevation_response_nostoc, filename = "documents/plots/elevation_response_nostoc.pdf",
    width = 10, height = 12, units = "cm")
ggsave(conifer_response_pelt, filename = "documents/plots/conifer_response_pelt.pdf",
    width = 10, height = 12, units = "cm")
ggsave(conifer_response_nostoc, filename = "documents/plots/conifer_response_nostoc.pdf",
    width = 10, height = 12, units = "cm")


##### COMPUTE AND PLOT PCA OF ENVIRONMENTAL RESPONSES #####

# Prepare the data for clustering, filtering out species from module 6
gjam_beta_wide <- gjam_beta %>%
  filter(module != 6) %>%
  select(species, variable, Estimate) %>%
  pivot_wider(names_from = variable, values_from = Estimate)

# Remove species with missing values
gjam_beta_wide <- gjam_beta_wide %>%
  drop_na()

# Extract the species and module information
species_info <- gjam_beta_wide %>%
  select(species) %>%
  left_join(module_assignments, by = "species")

# Convert module to character and then to factor
species_info <- species_info %>%
  mutate(module = as.character(module),
         module = factor(module))

# Extract the data for clustering
clustering_data <- gjam_beta_wide %>%
  select(-species)

# Perform PCA
pca_result <- prcomp(clustering_data, scale. = TRUE)

# Calculate the percentage of variation explained by PC1 and PC2
pca_var <- summary(pca_result)$importance[2, 1:2] * 100
pc1_var <- round(pca_var[1], 2)
pc2_var <- round(pca_var[2], 2)

# Create a data frame with PCA results
pca_df <- as.data.frame(pca_result$x) %>%
  mutate(species = gjam_beta_wide$species) %>%
  left_join(species_info, by = "species")

# Create a data frame for the variable loadings
loadings_df <- as.data.frame(pca_result$rotation) %>%
  rownames_to_column("variable")

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

# Define custom colors for the modules
module_colors <- c("1" = "#8a708a", "2" = "#7e937b", "3" = "#3277b0", "4" = "#be8551")

# Plot PCA results without species names and with custom colors and shapes
pca_plot <- ggplot(pca_df) +
  geom_point(size = 3, aes(x = PC1, y = PC2, color = module, shape = symbiont)) +
  scale_color_manual(values = module_colors) +
  scale_shape_manual(values = c("nostoc" = 16, "peltigera" = 17)) +  # Circles for Nostoc, diamonds for Peltigera
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "gray50") +
  geom_text(data = loadings_df, aes(x = PC1, y = PC2, label = variable), 
            color = "gray50", vjust = 1.5) +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)")) +
  custom_theme

# Save the plot as a PDF
ggsave(pca_plot, filename = "documents/plots/env_response_pca.pdf",
       width = 13, height = 10, units = "cm")


###### EXPLORE GJAM OUTPUT ######

# Plot the results
#gjamPlot(output = out1, plotPars = list(PLOTALLY = T, SAVEPLOTS = T, outfolder = "./"))
#out$fit$rmspeAll
#out$fit$DIC

# create base for plots
#base <- ggplot() +
#  theme_bw() +
#  theme(panel.background = element_blank(),
#        panel.grid = element_blank(),
#        plot.title = element_text(hjust = 0.5),
#        plot.subtitle = element_text(hjust = 0.5))

#base +
#  geom_point(aes(x = out$inputs$y[,16],
#                 y = out$prediction$ypredMu[,16])) +
#  labs(x = "observed", y = "predicted", title = paste0(colnames(out$inputs$y)[16], ", common"))

#xMu  <- out$prediction$xpredMu        # inverse prediction of x
#xSd  <- out$prediction$xpredSd
#yMu  <- out$prediction$ypredMu        # predicted y
#hold <- out$modelList$holdoutIndex    # holdout observations (rows)

#plot(out$inputs$xUnstand[hold,-1],xMu[hold,-1], cex=.2, xlab='True', ylab='Predictive mean')
#title('holdouts in x')
#abline( 0, 1, lty=2 )


###### VARIANCE PARTITIONING ######

# Extract parameters from GJAM output
residual_covariance <- out1[["parameters"]][["sigMu"]]  # Residual covariance matrix (S x S)
predictor_coefficients <- out1[["parameters"]][["betaMu"]]  # Beta coefficients (Q x S)
predictor_covariance <- cov(out1[["inputs"]][["xUnstand"]])  # Covariance of predictors (Q x Q)

# Number of species
species_names <- colnames(out1[["prediction"]][["ypredMu"]])
num_species <- length(species_names)

# Initialize an empty results table
variance_partitioning_table <- data.frame(
  Species = species_names,
  Total_Variance = numeric(num_species),
  Fraction_Explained_By_Environment = numeric(num_species),
  Fraction_Residual_Covariance = numeric(num_species)
)

# Compute variance components for each species
for (s in seq_len(num_species)) {
  # Environmental variance for species s: b_s^T * V * b_s
  beta_vector <- predictor_coefficients[, s]  # Coefficients for species s (Q x 1)
  environmental_variance <- sum(beta_vector %*% predictor_covariance %*% beta_vector)
  
  # Residual variance for species s: diagonal element of residual covariance
  residual_variance <- residual_covariance[s, s]
  
  # Total variance
  total_variance <- environmental_variance + residual_variance
  
  # Update the results table
  variance_partitioning_table$Total_Variance[s] <- total_variance
  variance_partitioning_table$Fraction_Explained_By_Environment[s] <- environmental_variance / total_variance
  variance_partitioning_table$Fraction_Residual_Covariance[s] <- residual_variance / total_variance
}


###### SARAH CODE ######

load("out_full_CPUE_fall.RData")

out <- out1

sigMu <- out[["parameters"]][["sigMu"]]
sensBeta <- out[["parameters"]][["sensBeta"]]

sensMu <- cbind(sensMu, sensBeta[,1])
sensSe <- cbind(sensSe, sensBeta[,1])
colnames(sensMu)[ncol(sensMu)] <- colnames(sensSe)[ncol(sensMu)] <- 'beta'

        sigma <- sqrt( diag( sigMu )) #sensBeta on sd scale
        sens  <- cbind(sensMu, sigma)
        
        sens <- sens^2
        
        sprop <- sweep( sens, 1, rowSums(sens), '/')
        ord <- order(sprop[,2], decreasing = T)
        
        smu <- t(sprop[ord,])
        
names <- colnames(out[["prediction"]][["ypredMu"]])
smu2 <- t(smu)  
    smu2 <- as.data.frame(smu2) 
    smu2$snames <- rownames(smu2)
    specs <- as.data.frame(names)
    colnames(specs) <- c("snames")
    smu2 <- left_join(smu2, specs, by = c("snames"))
    
#Plot 
    smu3 <- t(smu2)
    colnames(smu3) <- smu2$specs
        smax <- max( colSums(smu) )
        tmp0 <- barplot( smu, beside = F, xaxt = 'n',
                        ylim = c(-.5, smax), ylab = 'Proportion of total variance' )
        tmp0
        text( tmp0 - .2*diff(tmp0)[1], -.1, colnames(smu), srt = 90, pos = 1, cex=.6, col = as.vector(smu2$specColor))