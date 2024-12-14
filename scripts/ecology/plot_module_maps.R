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

##### TEMPERATURE BASE MAP #####

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

# Load the climatic rasters
mat_raster <- terra::rast("data/gis/climate_variables_2024/MAT.tif")

# Rescale for maximum constrast
mat_raster_stretched <- stretch(mat_raster, minq = 0.001, maxq = 0.999)

# Map rescaled T values to orginal ones and convert to DF for ggplot
mat_df <- as.data.frame(mat_raster, xy = TRUE, na.rm = TRUE)
mat_raster_stretched_df <- as.data.frame(mat_raster_stretched, xy = T, na.rm = T)
mat_df$stretched_mat <- mat_raster_stretched_df$MAT
colnames(mat_df) <- c("x", "y", "temperature", "stretched_mat")

# Base plot with MAT
mat_base_map <- ggplot() +
  # Plot hillshade as a grayscale raster
  geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade), alpha = 0.8) +
  scale_fill_gradient(low = "white", high = "black", name = "Hillshade", na.value = "transparent") +
  # Overlay precipitation raster
  geom_raster(data = mat_df, aes(x = x, y = y, fill = stretched_mat), alpha = 0.85) +
  scale_fill_viridis_c(option = "C", name = "Temperature", na.value = "transparent",
                       breaks = seq(0, 255, length.out = 5),  # Breaks in stretched scale (0 to 1)
                       labels = round(seq(min(mat_df$temperature), max(mat_df$temperature), length.out = 5), 1)) +
  # Add AB NRs shapefile overlay
  geom_sf(data = ab_nr_shp, fill = NA, color = "black", size = 0.7) +
  # Add labels
  labs(x = "Lon", y = "Lat") +
  # Add scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_minimal) +
  theme_minimal() + 
  theme(axis.text = element_text(colour = "black"),
        legend.text = element_text(colour = "black"))

##### ANNUAL PRECIPITATION BASE MAP #####

# Load the climatic rasters
precip_raster <- terra::rast("data/gis/climate_variables_2024/MAP.tif")

# Rescale for maximum constrast
#cmi_raster_stretched <- stretch(cmi_raster, minq = 0.001, maxq = 0.999)

# Map rescaled T values to orginal ones and convert to DF for ggplot
precip_df <- as.data.frame(precip_raster, xy = TRUE, na.rm = TRUE)
#cmi_raster_stretched_df <- as.data.frame(cmi_raster_stretched, xy = T, na.rm = T)
#cmi_df$stretched_cmi <- cmi_raster_stretched_df$CMI
colnames(precip_df) <- c("x", "y", "map")

# Base plot with MAT
precip_base_map <- ggplot() +
  # Plot hillshade as a grayscale raster
  geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade), alpha = 0.8) +
  scale_fill_gradient(low = "white", high = "black", name = "Hillshade",
                      guide = "none") +
  # Overlay precipitation raster
  new_scale_fill() +
  geom_raster(data = precip_df, aes(x = x, y = y, fill = map), alpha = 0.85) +
  scale_fill_gradientn(colors = c("#BF812C", "#66CCFF", "#0090B7", "#03526D", "#003C30"),
                       limits = c(278, 1000),
                       oob = scales::squish,) +
                       #values = scales::rescale(c(-42.34, 0, 10, 30, 100, 256.59))) +
  #scale_fill_viridis_c(option = "C", name = "CMI", na.value = "transparent") +
  # Add AB NRs shapefile overlay
  geom_sf(data = ab_nr_shp, fill = NA, color = "black", size = 0.7) +
  # Add labels
  labs(x = "Lon", y = "Lat") +
  # Add scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_minimal) +
  theme_minimal() + 
  theme(axis.text = element_text(colour = "black"),
        legend.text = element_text(colour = "black"))


##### CMI BASE MAP #####

# Load the climatic rasters
cmi_raster <- terra::rast("data/gis/climate_variables_2024/CMI.tif")

# Rescale for maximum constrast
#cmi_raster_stretched <- stretch(cmi_raster, minq = 0.001, maxq = 0.999)

# Map rescaled T values to orginal ones and convert to DF for ggplot
cmi_df <- as.data.frame(cmi_raster, xy = TRUE, na.rm = TRUE)
#cmi_raster_stretched_df <- as.data.frame(cmi_raster_stretched, xy = T, na.rm = T)
#cmi_df$stretched_cmi <- cmi_raster_stretched_df$CMI
colnames(cmi_df) <- c("x", "y", "cmi")

# Base plot with MAT
cmi_base_map <- ggplot() +
  # Plot hillshade as a grayscale raster
  geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade), alpha = 0.8) +
  scale_fill_gradient(low = "white", high = "black", name = "Hillshade", na.value = "transparent") +
  # Overlay precipitation raster
  new_scale_fill() +
  geom_raster(data = cmi_df, aes(x = x, y = y, fill = cmi), alpha = 0.85) +
  scale_fill_gradientn(colors = c("#BF812C", "#F5E7C3", "#0090B7", "#03526D", "#003C30"),
                       limits = c(-43, 40),
                       oob = scales::squish,) +
  #values = scales::rescale(c(-42.34, 0, 10, 30, 100, 256.59))) +
  #scale_fill_viridis_c(option = "C", name = "CMI", na.value = "transparent") +
  # Add AB NRs shapefile overlay
  geom_sf(data = ab_nr_shp, fill = NA, color = "black", size = 0.7) +
  # Add labels
  labs(x = "Lon", y = "Lat") +
  # Add scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_minimal) +
  theme_minimal() + 
  theme(axis.text = element_text(colour = "black"),
        legend.text = element_text(colour = "black"))

##### PLOT SITES BY INTERACTION MODULE #####

# Load ABMI site coordinates
abmi_sites_coord <- read_csv("data/tables/abmi_site_data.csv") %>%
  distinct(Site, .keep_all = T) %>%
  select(Site, Lat, Long)

# Load module assignments data
load("analyses/ecology/peltigera_module_assignments.RData")

# Load ABMI ID data
module_data <- read_csv("data/tables/abmi_id_data.csv") %>%
  filter(!is.na(mycobiont_molecular_id), !is.na(nostoc_otu),
         str_detect(mycobiont_molecular_id, "Peltigera.*")) %>%
  group_by(site, mycobiont_molecular_id, nostoc_otu) %>%
  summarise(frequency = n()) %>%
  ungroup() %>%
  left_join(peltigera_module_assignments, by = "mycobiont_molecular_id") %>%
  left_join(nostoc_module_assignments, by = "nostoc_otu") %>%
  rename(peltigera_module = module.x,
         nostoc_module = module.y)

# Get within-module interactions counts per site
within_module_counts <- module_data %>%
  filter(peltigera_module == 4 & nostoc_module == 4) %>%
  group_by(site) %>%
  summarise(n_specimens = sum(frequency)) %>%
  ungroup() %>%
  mutate(Site = str_remove(site, "ABMI ")) %>%
  filter(!str_detect(Site, "SK")) %>%                # Removing SK sites
  left_join(abmi_sites_coord, by = "Site") %>%
  select(Lat, Long, n_specimens) %>%
  transform_crs()

# Peltigera switch module
peltigera_switch_counts <- module_data %>%
  filter(peltigera_module == 2 & nostoc_module != 2) %>%
  group_by(site) %>%
  summarise(n_specimens = sum(frequency)) %>%
  ungroup() %>%
  mutate(Site = str_remove(site, "ABMI ")) %>%
  filter(!str_detect(Site, "SK")) %>%                # Removing SK sites
  left_join(abmi_sites_coord, by = "Site") %>%
  select(Lat, Long, n_specimens) %>%
  transform_crs()

# Nostoc switch module
peltigera_switch_counts <- module_data %>%
  filter(peltigera_module == 2 & nostoc_module != 2) %>%
  group_by(site) %>%
  summarise(n_specimens = sum(frequency)) %>%
  ungroup() %>%
  mutate(Site = str_remove(site, "ABMI ")) %>%
  filter(!str_detect(Site, "SK")) %>%                # Removing SK sites
  left_join(abmi_sites_coord, by = "Site") %>%
  select(Lat, Long, n_specimens) %>%
  transform_crs()

mat_base_map +
geom_sf(data = within_module_counts, aes(size = n_specimens), 
        fill = "white", color = "black", shape = 21, stroke = 0.2, alpha = 1) +
  scale_size_continuous(name = "Number of Specimens", limits = c(1, 32))


ggsave(filename = "documents/plots/map_module4_mat.pdf",
       width = 11, height = 10, units = "cm")
