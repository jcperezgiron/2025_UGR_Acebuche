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

    mop_results <- mop(
      m = m_matrix,
      g = future_predictors,
      type = "basic",
      calculate_distance = TRUE,
      rescale_distance = TRUE
    )

    mop_results <- mop_results$mop_distances

    pretty_scenario <- ifelse(scenario == "245", "SSP2-4.5 ",
      ifelse(scenario == "370", "SSP3-7.0 ", "SSP5-8.5 ")
    )

    names(mop_results) <- paste0(pretty_scenario, year)
    mop <- c(mop, mop_results)
  }
}

for (i in 1:nlyr(mop)) {
  writeRaster(mop[[i]], paste0("C:/SCIENCE/2025_UGR_Acebuche/MOP_results/", str_replace(names(mop[[i]]), " ", "_"),".tif"), overwrite = TRUE)
}

#### 4.- Plot MOP ####
mop_plot <- function(lyr) {
  ggplot() +
    geom_spatraster(data = mop[[lyr]]) +
    scale_fill_whitebox_c(
      palette = "muted", direction = 1,
      breaks = seq(0, 1, 0.25), limits = c(0, 1),
      na.value = "transparent"
    ) +
    labs(fill = "MOP distance", title = names(mop[[lyr]]), linetype = "", color = "") +
    guides(
      fill = guide_colourbar(order = 1, title.position = "top", title.hjust = 0.5, keywidth = unit(4, "cm")),
      linetype = guide_legend(ncol = 1),
      color = guide_legend(ncol = 1)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
    )
}

p1 <- mop_plot(1)
p2 <- mop_plot(2)
p3 <- mop_plot(3)
p4 <- mop_plot(4)
p5 <- mop_plot(5)
p6 <- mop_plot(6)

Fig_mop <- ggarrange(p1, p4, p2, p5, p3, p6,
                     ncol = 2, nrow = 3,
                     align = "hv", common.legend = TRUE, legend = "bottom"
)

ggsave("C:/SCIENCE/2025_UGR_Acebuche/Manuscript/Figure_MOP.jpg", Fig_mop, units = "mm", dpi = 300, width = 183, height = 183)
