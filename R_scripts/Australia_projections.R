# Script to apply the species distribution model to Australia

# Load the required libraries
library(tidyverse)
library(sf)
library(biomod2)
library(terra)
library(usdm)
library(rnaturalearth)

# Download Australia
australia <- ne_countries(country = "australia",
                          scale = "medium",
                          returnclass = "sf")

# Bounding box
bbox <- st_bbox(australia)
print(bbox)

# as sf
bbox_sf <- st_as_sfc(bbox)

# plot
ggplot() +
  geom_sf(data = australia, fill = "lightgreen", color = "darkgreen") +
  geom_sf(data = bbox_sf, fill = NA, color = "red", linetype = "dashed") +
  labs(title = "Australia y su bounding box") +
  theme_minimal()


# Load predictors
predictor_list <- list.files("C:/SCIENCE/2025_UGR_Acebuche/datos_proyecciones_australia/", pattern = ".tif$", full.names = T)

# Stack predictors
myExpl <- rast(predictor_list) %>% 
  crop(bbox)
names(myExpl) <- str_replace(names(myExpl), "_", "")
names(myExpl) <- toupper(names(myExpl))


# Load model 
myBiomodModelOut <- get(load("C:/SCIENCE/2025_UGR_Acebuche/SDM_output/olive/olive.1745443608.models.out"))
myBiomodEM <- get(load("C:/SCIENCE/2025_UGR_Acebuche/SDM_output/olive/olive.1745443608.ensemble.models.out"))

#### Project current models in Australia
# Project single models
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "Current_Australia",
  new.env = myExpl,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  build.clamping.mask = TRUE
)

myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj,
  proj.name = "Current_EM_Australia",
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  on_0_1000 = FALSE,
  nb.cpu = 8,
  na.rm = FALSE
)




# Load projection results
proj <- rast("C:/SCIENCE/2025_UGR_Acebuche/SDM_output/olive/proj_Current_EM_Australia/proj_Current_EM_Australia_olive_ensemble.tif")
proj <- proj/1000


sf <- read.csv("C:/SCIENCE/2025_UGR_Acebuche/datos_proyecciones_australia/poblaciones_acebuche_aus_cornuault_2015.csv") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  rename(ID_pop = ID)

hs_values <- terra::extract(proj, sf) %>% 
  cbind(sf) %>% 
  dplyr::select(-ID) %>% 
  mutate(lon = st_coordinates(sf)[,1],
         lat = st_coordinates(sf)[,2]) %>%
  st_drop_geometry() %>% 
  dplyr::select(-geometry) %>%
  rename(ID = ID_pop) %>% 
  dplyr::select(ID, lon, lat, everything())

colnames(hs_values)[4] <- "hs"

# Plot hs values with boxplot
ggplot(hs_values, aes(x = "", y = hs)) +
  geom_violin(fill = "lightgreen") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkred") +
  labs(title = "Distribución de valores de hábitat potencial\nen poblaciones de acebuche en Australia",
       y = "Idoneidad del hábitat potencial",
       x = "") +
  theme_minimal()


write.csv(hs_values, "C:/SCIENCE/2025_UGR_Acebuche/Australia_outputs/hs_values_australia.csv", row.names = FALSE)
