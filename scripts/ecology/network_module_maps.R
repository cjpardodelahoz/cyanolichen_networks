# Load required libraries
library(sf)
library(ggplot2)

# Use the example shapefile included in the sf package (North Carolina dataset)
nc <- st_read(system.file("shape/nc.shp", package = "sf"))

# Print a summary of the shapefile
print(nc)

# Create a dataframe with points (example locations in North Carolina)
points_df <- data.frame(
  longitude = c(-81.0, -79.5, -76.0),
  latitude = c(35.5, 36.0, 34.8)
)

# Plot the shapefile and overlay the points
ggplot() +
  geom_sf(data = nc, fill = "lightblue", color = "black", alpha = 0.5) + # Shapefile
  geom_point(data = points_df, aes(x = longitude, y = latitude), 
             color = "red", size = 3) + # Points
  labs(title = "North Carolina Map with Points",
       x = "Longitude", y = "Latitude") +
  theme_minimal()




ggplot(data = ab) +
  geom_sf(aes(fill = NRNAME))


library(sf)
library(terra)
library(ggspatial)
library(ggplot2)
library(dplyr)
library(ggnewscale)

# Load hillshade file (GeoTIFF format)
hillshade <- rast("data/gis/Hillshade.tif")

ab <- st_read("data/gis/Natural_Regions_Subregions_of_Alberta.shp")

print(ab)

# Convert hillshade to a dataframe for ggplot2 compatibility
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
colnames(hillshade_df) <- c("x", "y", "hillshade") # Ensure column names are clear

# Dissolve the shapefile by NRNAME (removing subregion boundaries)
ab_nrname <- ab %>%
  group_by(NRNAME) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")

# Plot with boundaries only for NRNAME
ggplot(data = ab_nrname) +
  geom_sf(aes(fill = NRNAME), color = "black", size = 0.5) + # Boundaries only for NRNAME
  scale_fill_viridis_d(option = "C", name = "NRNAME") +
  labs(title = "Map Filled by NRNAME with NRNAME Boundaries Only", 
       x = "Longitude", y = "Latitude") +
  theme_minimal()


# Step 3: Plot hillshade and overlay shapefile
ggplot() +
  # Hillshade as background
  geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade)) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  # Reset the fill scale for shapefile regions
  new_scale_fill() +
  # Overlay shapefile with transparent fill
  geom_sf(data = ab_nrname, aes(fill = NRNAME), color = NA, linewidth = 0.1, alpha = 0.3) +
  labs(x = "Lat", y = "Lon") +
  # Adjust alpha for transparency
  # Add a color scale for the shapefile's regions
  scale_fill_viridis_d(option = "C", name = "Natural Regions") +
  # Add scale bar
  annotation_scale(location = "bl", width_hint = 0.3) + # "bl" = bottom left
  
  # Add north arrow
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_minimal) +
  theme_minimal()

