#!/usr/bin/env Rscript

##### LOAD PACKAGES AND CUSTOM FUNCTIONS #####

# Load required packages
library(sf)
library(terra)
library(ggspatial)
library(ggplot2)
library(dplyr)
library(readr)
library(magrittr)
library(ggnewscale)
library(stringr)

# Function to transform coordinates to AB shape file CRS
transform_crs <- function(lat_long_df) {
  lat_long_df %>%
  # Convert to sf data
  st_as_sf(., coords = c("Long", "Lat"), crs = 4326) %>%
  # Transform to match CRS of AB shapefile
  st_transform(., crs = st_crs(ab_nr_shp)) # transforming specifically to the ABM shapefile with natural subregion boundaries removed
}

##### BASE AB MAP #####

# Load hillshade file and convert to dataframe for ggplot 2 compatibility
hillshade <- rast("data/gis/Hillshade.tif")
hillshade <- stretch(hillshade, minq = 0.02, maxq = 0.98) # To maximize contrast in hillshade
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
colnames(hillshade_df) <- c("x", "y", "hillshade")

# Load AB shapefile and remove subregion boundaries
ab <- st_read("data/gis/Natural_Regions_Subregions_of_Alberta.shp")
ab_nr_shp <- ab %>%
  group_by(NRNAME) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")

#print(ab)
#print(ab_nr_shp)

# Colors for natural regions
nr_colors <- c("Boreal" = "#93C192",
               "Canadian Shield" = "#A38682",
               "Foothills" = "#9C97CE",
               "Grassland" = "#E2C68D",
               "Parkland" = "#D3965F",
               "Rocky Mountain" = "#9FB4C1")

# Base map plot
base_map <- ggplot() +
  # Hillshade as background
  geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade)) +
  scale_fill_gradient(low = "gray30", high = "white", guide = "none") +
  # Reset the fill scale for shapefile regions
  new_scale_fill() +
  # Overlay shapefile with transparent fill
  geom_sf(data = ab_nr_shp, aes(fill = NRNAME), color = NA, linewidth = 0.1, alpha = 0.5) +
  labs(x = "Lat", y = "Lon") +
  scale_fill_manual(values = nr_colors, name = "Natural region") +
  # Add scale bar
  annotation_scale(location = "bl", width_hint = 0.3) + # "bl" = bottom left
  # Add north arrow
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_minimal) +
  theme_minimal() + 
  theme(axis.text = element_text(colour = "black"),
        legend.text = element_text(colour = "black"))

##### PLOT ALL AB SITES INCLUDED #####

# Set seed for reproducibility
set.seed(42) 

# Load ABMI site data (regional) and transform to reference CRS
all_abmi_sites <- read_csv("data/tables/abmi_site_data.csv") %>%
  # Remove duplicated sites
  filter(!duplicated(Site)) %>%
  select(Lat, Long)

# Identify duplicate coordinates
is_duplicate <- duplicated(all_abmi_sites) | duplicated(all_abmi_sites, fromLast = TRUE)

# Add noise to duplicated rows
all_abmi_sites[is_duplicate, ] <- all_abmi_sites[is_duplicate, ] + 
  matrix(runif(sum(is_duplicate) * 2, -0.1, 0.1), ncol = 2)

# Transform crs
all_abmi_sites <- transform_crs(all_abmi_sites)

# 14 sites from regional also included in local dataset
local_sites <- c("1116",
                 "1077",
                 "1169",
                 "1293",
                 "1378",
                 "1427",
                 "1529",
                 "1631",
                 "1499",
                 "1307",
                 "OG-ABMI-1122-1",
                 "877",
                 "OG-DH-751-1",
                 "OG-DH-785-1")

# Get coordinates for local sites
local_sites_df <- read_csv("data/tables/abmi_site_data.csv") %>%
  filter(Site %in% local_sites) %>%
  select(Lat, Long)

# Add coordinates for site 1321, which was not part of regional dataset
local_sites_df[15,1] <- 52.3611
local_sites_df[15,2] <- -116.3646

# Trnasform local sites to AB CRS
local_sites_sf <- transform_crs(local_sites_df)

# Plot all sites
map_with_all_sites <- base_map +
  # Add points
  geom_sf(data = all_abmi_sites, fill ="white", color = "black", size = 1.5,
          stroke = 0.2, shape = 21) +
  geom_sf(data = local_sites_sf, fill = "black", color = "black", size =1.5,
          stroke = 0.2, shape = 21)

# Save the map
ggsave(map_with_all_sites, filename = "documents/plots/map_with_all_sites.pdf",
       width = 13, height = 14, units = "cm")