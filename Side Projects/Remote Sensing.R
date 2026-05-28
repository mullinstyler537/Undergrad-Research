# Project: Satellite Wildfire Data

# --- Load in the libraries ---
library(rstac)
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(leaflet)

# --- Define Area of Interest ---
bbox_co <- c(
  xmin = -105.70, 
  ymin = 40.55, 
  xmax = -105.40, 
  ymax = 40.75
)

# --- Query Satellite API for a verified crisp, clear summer day ---
stac_client <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")

stac_query <- stac_client %>%
  stac_search(
    collections = "sentinel-2-l2a",
    bbox = bbox_co,
    datetime = "2020-08-20/2020-08-25", # Targeted window with perfect visibility
    limit = 5
  ) %>%
  post_request()

signed_items <- items_sign(stac_query, sign_fn = sign_planetary_computer())
selected_scene <- signed_items$features[[1]]

# --- Fetch & Crop Raster Bands ---
red_url <- selected_scene$assets$B04$href
nir_url <- selected_scene$assets$B08$href

red_raster <- rast(paste0("/vsicurl/", red_url))
nir_raster <- rast(paste0("/vsicurl/", nir_url))

bbox_sf <- st_bbox(bbox_co, crs = st_crs(4326)) %>% 
  st_as_sfc()

bbox_projected <- st_transform(bbox_sf, crs(red_raster))

red_crop <- crop(red_raster, bbox_projected)
nir_crop <- crop(nir_raster, bbox_projected)

# --- Calculate NDVI with float normalization ---
ndvi_pre <- (nir_crop - red_crop) / (nir_crop + red_crop)

# --- Dynamic Map Render (Auto-scales colors to available pixel values) ---
ggplot() +
  geom_spatraster(data = ndvi_pre) +
  scale_fill_gradientn(
    colors = c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571"),
    na.value = "transparent",
    name = "NDVI"
  ) +
  theme_minimal() +
  labs(
    title = "Pre-Fire Vegetation Health (August 2020)",
    subtitle = "Area: Arapaho/Roosevelt National Forests, CO",
    x = "Longitude",
    y = "Latitude"
  )

# --- Query Satellite API for Post-Fire Imagery (August 2021) ---
stac_query_post <- stac_client %>%
  stac_search(
    collections = "sentinel-2-l2a",
    bbox = bbox_co,
    datetime = "2021-08-20/2021-08-25", 
    limit = 5
  ) %>%
  post_request()

signed_items_post <- items_sign(stac_query_post, sign_fn = sign_planetary_computer())
selected_scene_post <- signed_items_post$features[[1]]

# --- Fetch, Crop, and Calculate Post-Fire NDVI ---
red_url_post <- selected_scene_post$assets$B04$href
nir_url_post <- selected_scene_post$assets$B08$href

red_raster_post <- rast(paste0("/vsicurl/", red_url_post))
nir_raster_post <- rast(paste0("/vsicurl/", nir_url_post))

red_post <- crop(red_raster_post, bbox_projected)
nir_post <- crop(nir_raster_post, bbox_projected)

ndvi_post <- (nir_post - red_post) / (nir_post + red_post)

# --- Calculate Burn Severity (NDVI Change Detection) ---
ndvi_difference <- ndvi_post - ndvi_pre

# --- Plot the Final Burn Severity Disturbance Map ---
ggplot() +
  geom_spatraster(data = ndvi_difference) +
  scale_fill_gradientn(
    colors = c("#990000", "#d7301f", "#feb24c", "#f7f7f7", "#006837"),
    limits = c(-0.5, 0.2), 
    na.value = "transparent",
    name = "NDVI Change"
  ) +
  theme_minimal() +
  labs(
    title = "Wildfire Burn Severity & Forest Disturbance Map",
    subtitle = "Cameron Peak Fire Impact (2020 vs 2021) | Arapaho National Forest",
    x = "Longitude",
    y = "Latitude"
  )  

# --- Combine All Three Layers Into a Single SpatRaster ---
all_layers <- c(ndvi_pre, ndvi_post, ndvi_difference)

# --- Rename the layers so they look clean on the plot labels ---
names(all_layers) <- c("A: Pre-Fire Forest Health (2020)", 
                       "B: Post-Fire Forest Health (2021)", 
                       "C: Calculated Burn Severity Scar")

# --- Plot the Side-by-Side Faceted Map with Text Fixes ---
ggplot() +
  geom_spatraster(data = all_layers) +
  scale_fill_gradientn(
    colors = c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b", 
               "#ffffbf", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"),
    limits = c(-0.6, 0.85),
    na.value = "transparent",
    name = "NDVI / Change"
  ) +
  facet_wrap(~lyr, ncol = 3) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    # FIX: Rotates the X-axis numbers by 45 degrees so they clear each other cleanly
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    title = "Satellite Remote Sensing Disturbance Dashboard",
    subtitle = "Cameron Peak Fire Analysis | Sentinel-2 Multispectral Imagery (10m Resolution)",
    x = "Longitude",
    y = "Latitude"
  )

# --- Extract Raster Data to a Dataframe ---
df_pre  <- as.data.frame(ndvi_pre, xy = FALSE) %>% drop_na() %>% mutate(Period = "2020: Pre-Fire")
df_post <- as.data.frame(ndvi_post, xy = FALSE) %>% drop_na() %>% mutate(Period = "2021: Post-Fire")

# Rename the columns so they match perfectly for combining
colnames(df_pre)[1]  <- "NDVI"
colnames(df_post)[1] <- "NDVI"

# Combine into a single long dataframe
comparison_df <- bind_rows(df_pre, df_post)

# --- Plot an Ecological Distribution Graph ---
ggplot(comparison_df, aes(x = NDVI, fill = Period)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("2020: Pre-Fire" = "#238b45", "2021: Post-Fire" = "#d7301f")) +
  theme_minimal() +
  labs(
    title = "Shift in Forest Canopy Health Distribution",
    subtitle = "Density curve showing the massive shift toward lower NDVI values post-fire",
    x = "NDVI Value (Vegetation Density)",
    y = "Pixel Density"
  )

# --- Reproject and Downsample for the Web Interface ---
# Project to standard web mercator projection
ndvi_diff_web <- project(ndvi_difference, "EPSG:3857")

# Downsample the raster by a factor of 2 (aggregates pixels to keep file size small)
ndvi_diff_leaflet <- aggregate(ndvi_diff_web, fact = 2, fun = mean)

# --- Build the Fixed Color Palette Function ---
# Extract raw numeric values to avoid scale warnings
pixel_values <- values(ndvi_diff_leaflet, mat = FALSE)

pal <- colorNumeric(
  palette = c("#990000", "#d7301f", "#feb24c", "#f7f7f7", "#006837"), 
  domain = c(-0.5, 0.2), 
  na.color = "transparent"
)

# --- Render the Interactive Leaflet Map ---
leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addRasterImage(ndvi_diff_leaflet, colors = pal, opacity = 0.7) %>%
  addLegend(
    pal = pal, 
    values = pixel_values, 
    title = "NDVI Change",
    position = "bottomright"
  )
