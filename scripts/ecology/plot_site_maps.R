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

# Load ABMI site data (regional) and transform to reference CRS
all_abmi_sites <- read_csv("data/tables/abmi_site_data.csv") %>%
  # Remove duplicated sites
  filter(!duplicated(Site)) %>%
  select(Lat, Long) %>%
  # Transform crs
  transform_crs()

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

##### PLOT SITES BY INTERACTION MODULE #####




# Category 1 (6 colors)
c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

# Category 2 (11 colors)
c("#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3")

# Category 3 (6 colors)
c("#FFED6F", "#D9D9D9", "#BC80BD", "#CCEBC5", "#DECBE4", "#FED9A6")



# Category 1 (6 colors)
c("#003F5C", "#58508D", "#BC5090", "#FF6361", "#FFA600", "#5C8001")

# Category 2 (11 colors)
c("#D73027", "#FC8D59", "#FEE08B", "#91CF60", "#1A9850", "#4575B4", 
  "#313695", "#A6D96A", "#F46D43", "#FEE08B", "#D53E4F")

# Category 3 (6 colors)
c("#6A51A3", "#FEC44F", "#FC9272", "#FF9896", "#E41A1C", "#C7E9B4")




# Category 1 (6 colors, earthy greens and blues)
c("#8C510A", "#BF812D", "#DFC27D", "#80CDC1", "#35978F", "#01665E")

# Category 2 (11 colors, bright rainbow colors)
c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62")

# Category 3 (6 colors, muted warm tones)
c("#FEE0D2", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#67000D")


# Category 1 (6 colors, greens and blues)
c("#4CAF50", "#81C784", "#A5D6A7", "#E8F5E9", "#388E3C", "#2E7D32")

# Category 2 (11 colors, oranges and purples)
c("#FF9800", "#FFB74D", "#FFE0B2", "#FFCCBC", "#E64A19", "#BF360C", 
  "#CE93D8", "#BA68C8", "#9C27B0", "#6A1B9A", "#4A148C")

# Category 3 (6 colors, reds and yellows)
c("#FFEB3B", "#FFF59D", "#FDD835", "#FFCA28", "#FFB300", "#FF6F00")

  
