#Project: Mapping Bighorn Sheep in the US

#Make sure to install the packages below
#install.packages(c("rgbif", "sf", "tidyverse", "maps", "leaflet))

library(rgbif)
library(sf)
library(tidyverse)
library(maps)
library(leaflet)

#Fetch Data
bighorn_species <- name_backbone(name = "Ovis canadensis")

raw_bighorn_data <- occ_search(
  scientificName = "Ovis canadensis",
  limit = 1000,
  hasCoordinate = TRUE
)

bighorn_df <- raw_bighorn_data$data

#Clean Data
bighorn_clean <- bighorn_df %>%
  select(species, decimalLatitude, decimalLongitude, year, stateProvince, countryCode) %>%
  filter(countryCode == "US") %>%
  drop_na(decimalLatitude, decimalLongitude)

#Spatial Conversion
bighorn_spatial <- st_as_sf(
  bighorn_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

us_states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

#Plot Map
ggplot() +
  geom_sf(data = us_states, fill = "gray95", color = "gray70") +
  geom_sf(data = bighorn_spatial, color = "darkorange", size = 1.5, alpha = 0.6) +
  coord_sf(xlim = c(-125, -100), ylim = c(31, 49)) +
  theme_minimal() +
  labs(
    title = "Bighorn Sheep (Ovis canadensis) Observations",
    subtitle = "Data Source: GBIF Sighting Records",
    x = "Longitude",
    y = "Latitude"
  )

#Install packages 
install.packages(c("geodata", "terra"))
library(geodata)
library(terra)

#Download global climate rasters at a 10 minute spatial resolution
climate_global <- worldclim_global(var = "bio", res = 10, path = tempdir())

#Subset to just our two variables: Bio1 (Mean Temp) and Bio12 (Annual Rain)
#Note: WorldClim stores temperature multiplied by 10 to save file size, so divide by 10
bioclim_layers <- climate_global[[c("wc2.1_10m_bio_1", "wc2.1_10m_bio_12")]]
names(bioclim_layers) <- c("Annual_Mean_Temp", "Annual_Precipitation")

#Extract the underlying pixel data for our spatial points
sheep_climate_values <- terra::extract(bioclim_layers, bighorn_spatial)

# Combine the climate values back with our cleaned bighorn table
bighorn_with_climate <- bighorn_clean %>% 
  mutate(
    Temp_C = sheep_climate_values$Annual_Mean_Temp / 10,
    Precip_mm = sheep_climate_values$Annual_Precipitation
  ) %>% 
  drop_na(Temp_C, Precip_mm)

# Inspect the new dataframe
head(bighorn_with_climate)

#Calculate the 5th and 95th percentiles for temperature and rain
niche_stats <- bighorn_with_climate %>%
  summarize(
    Min_Temp = quantile(Temp_C, 0.05),
    Max_Temp = quantile(Temp_C, 0.95),
    Min_Precip = quantile(Precip_mm, 0.05),
    Max_Precip = quantile(Precip_mm, 0.95)
 )

#View the climate limits of the Bighorn Sheep
print(niche_stats)

# Plot the climate niche envelope
ggplot(data = bighorn_with_climate, aes(x = Temp_C, y = Precip_mm)) +
  # Layer 1: Draw the shaded climate niche box
  geom_rect(
    aes(xmin = niche_stats$Min_Temp, xmax = niche_stats$Max_Temp, 
        ymin = niche_stats$Min_Precip, ymax = niche_stats$Max_Precip),
    fill = "lightblue", alpha = 0.02, color = "blue", linetype = "dashed"
  ) +
  # Layer 2: Plot the individual sheep sightings
  geom_point(color = "darkorange", alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(
    title = "Bighorn Sheep Ecological Climate Niche",
    subtitle = "Shaded box represents the core 90% climate envelope",
    x = "Annual Mean Temperature (°C)",
    y = "Annual Precipitation (mm)"
  )

# Create an interactive web map
interactive_sheep_map <- leaflet(data = bighorn_spatial) %>%
  # Add a topographic terrain background
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  
  # Add our bighorn sheep sighting points (FIX: Added a comma after radius = 4)
  addCircleMarkers(
    radius = 4,
    color = "darkorange",
    stroke = FALSE, 
    fillOpacity = 0.7,
    # Create a popup that displays information when a user clicks a point
    popup = ~paste0(
      "<strong>Species:</strong> ", species, "<br>",
      "<strong>State:</strong> ", stateProvince, "<br>",
      "<strong>Year Observed:</strong> ", year
    )
  )

# View the interactive map in your RStudio 'Viewer' tab
interactive_sheep_map