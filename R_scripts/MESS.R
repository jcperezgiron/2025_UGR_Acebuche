# Run MOP

# Load the required libraries
library(tidyverse)
library(sf)
library(biomod2)
library(terra)
library(tidyterra)
library(usdm)
library(mop)
library(ggpubr)

library(raster)

#### 1.- Read presences ####
presences <- read.csv("C:/SCIENCE/2025_UGR_Acebuche/data/presences/maxent oleaster 2-24.csv") %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  mutate(olive = 1) %>%
  select(olive)

#### 2.- Read environmental data ####

predictor_list <- list.files("C:/SCIENCE/2025_UGR_Acebuche//data/spatial_predictors/", pattern = ".asc$", full.names = T)

for (i in 1:length(predictor_list)) {
  if (i == 1) {
    template <- rast(predictor_list[i])
    rasStack <- rast(predictor_list[i])
    names(rasStack) <- str_replace_all(names(rasStack), pattern = "_clip", replacement = "")
  } else {
    raster <- rast(predictor_list[i]) %>%
      crop(template) %>%
      resample(template)
    names(raster) <- str_replace_all(names(raster), pattern = "_clip", replacement = "")
    rasStack <- c(rasStack, raster)
  }
}

# Generate 1000 randon points in rasStack
absences <- spatSample(rasStack, size = 10000, method = "random", na.rm = TRUE, as.points = TRUE, values = FALSE) %>%
  st_as_sf() %>%
  mutate(olive = 0)

#### 3.- Combine presences and absences ####
m_matrix <- terra::extract(rasStack, rbind(presences, absences)) %>%
  dplyr::select(-ID)

#### 7.- Future conditions  and compute mop ####
scenarios <- c("245", "370", "585") # RCP scenarios
years <- c("2060", "2080") # years to project

mop <- rast()
for (year in years) {
  for (scenario in scenarios) {
    # read climate projections and correct names to match rasStack_1
    list.climProj <- list.files(file.path("C:/SCIENCE/2025_UGR_Acebuche/data/spatial_predictors_CMMC-ESM2", paste(year, scenario, sep = "_")), pattern = ".tif$", full.names = T, recursive = T)
    future_predictors <- rast(c(list.climProj))
    names(future_predictors) <- str_split_i(names(future_predictors), "_", 1)
    
    mess <- dismo::mess(stack(future_predictors), m_matrix, full = T)
    
  }
}
