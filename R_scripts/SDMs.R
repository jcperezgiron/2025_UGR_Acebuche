# Script to perform a species distribution model using the Maxent algorithm

# Load the required libraries
library(tidyverse)
library(sf)
library(biomod2)
library(terra)
library(usdm)

#### Global params ####

# initial time
tmp <- Sys.time()

# define the output directory
output_dir <- file.path("./SDM_output")

# create output directory if it does not exist.
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# make a copy of maxent.jar to the output directory
file.copy("./maxent/maxent.jar", output_dir)

# set the working directory to the output directory
setwd(output_dir)


#### 1.- Read presences ####
presences <- read.csv("../data/presences/maxent oleaster 2-24.csv") %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  mutate(olive = 1) %>%
  select(olive)

#### 2.- Read environmental data ####

predictor_list <- list.files("../data/spatial_predictors/", pattern = ".asc$", full.names = T)

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

#### 3.- Check multicollinearity ####
rasValue <- presences %>%
  terra::extract(rasStack, .) %>%
  as.data.frame() %>% # Convierto a df
  dplyr::select(-ID)

# By spearman correlation of 0.7
v1 <- vifcor(rasValue, th = 0.7, method = "spearman")
v1

rasStack_1 <- exclude(rasStack, v1) # remove the variables with high correlation

#### 4.- Modelling parameters ####
# Select the name of the studied species
myRespName <- "olive"

# Get corresponding presence/absence data
myResp <- presences %>%
  vect()

# Environmental variables selected after correlation
myExpl <- rasStack_1

# PARAMETERS
myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp,
  expl.var = myExpl,
  resp.name = myRespName,
  PA.nb.rep = 10, # 1 to test processing time; 10 to run models
  PA.nb.absences = 10000, # Random number; 10000 to run models
  PA.strategy = "random",
  na.rm = T,
  filter.raster = F,
  dir.name = getwd()
)

myBiomodData

# Create the different validation datasets
# k-fold selection
cv.k <- bm_CrossValidation(
  bm.format = myBiomodData,
  strategy = "kfold",
  nb.rep = 1,
  k = 10
)

# Specify the models to be run
models <- c("MAXENT") # species presence-only methods

# Configure the models
user.MAXENT <- list("for_all_datasets" = list(
  memory_allocated = 1024
))

user.val <- list(MAXENT.binary.MAXENT.MAXENT = user.MAXENT)

opt.u <- bm_ModelingOptions(
  data.type = "binary",
  models = models,
  strategy = "user.defined",
  user.val = user.val,
  bm.format = myBiomodData,
  calib.lines = cv.k
)



#### 5.- Model fitting ####
## Single models
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = models,
  OPT.strategy = "user.defined",
  OPT.user.val = user.val,
  OPT.user.base = "default",
  CV.do.full.models = FALSE,
  CV.strategy = "user.defined",
  CV.user.table = cv.k,
  var.import = 10,
  metric.eval = c("BOYCE"),
  prevalence = NULL,
  seed.val = 123,
  do.progress = F,
  nb.cpu = 8
)

myBiomodModelOut

# Get evaluation scores & variables importance
single_eval <- get_evaluations(myBiomodModelOut)
single_eval

single_varimp <- get_variables_importance(myBiomodModelOut)

# save the evaluation scores& variables importance
write.csv(single_eval, file.path(output_dir, "single_eval.csv"), row.names = F)
write.csv(single_varimp, file.path(output_dir, "single_varimp.csv"), row.names = F)

# Ensemble models
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = "all",
  em.by = "algo",
  em.algo = c("EMmedian"),
  metric.select = c("BOYCE"),
  metric.select.thresh = c(0.7),
  metric.eval = c("BOYCE"),
  var.import = 10,
  seed.val = 123,
  nb.cpu = 8
)
myBiomodEM

# Get evaluation scores & variables importance
ensemble_eval <- get_evaluations(myBiomodEM)
ensemble_eval

ensemble_varimp <- get_variables_importance(myBiomodEM)

# save the evaluation scores
write.csv(ensemble_eval, file.path(output_dir, "ensemble_eval.csv"), row.names = F)
write.csv(ensemble_varimp, file.path(output_dir, "ensemble_varimp.csv"), row.names = F)

#### 6.- Project current models ####
# Project single models
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "Current",
  new.env = myExpl,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  build.clamping.mask = TRUE
)

myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj,
  proj.name = "Current_EM",
  models.chosen = "all",
  metric.binary = "BOYCE",
  metric.filter = "BOYCE",
  on_0_1000 = FALSE,
  nb.cpu = 8,
  na.rm = FALSE
)


myBiomodEMProj


#### 7.- Project onto future conditions ####

##
## This is an example
## This code is not going to work as is
## We need to adapt this code to your specific case
## 

# # Read the future environmental data
# years <- c("2050", "2070")
# models <- c("CCSM4", "GFDL-ESM2M", "HadGEM2-ES", "MIROC5")
# scenarios <- c("rcp26", "rcp45", "rcp85")
# 
# 
# for (year in years) {
#   for (scenario in scenarios) {
#     for (model in models) {
#       predictor_list_future <- list.files("../data/spatial_predictors_future/", pattern = ".asc$", full.names = T)
#       myExplFuture <- rast(predictor_list_future)
# 
#       # Project onto future conditions
#       myBiomodProjectionFuture <- BIOMOD_Projection(
#         bm.mod = myBiomodModelOut,
#         proj.name = paste0(year, "_", model, "_", scenario),
#         new.env = myExplFuture,
#         models.chosen = "all",
#         metric.binary = "all",
#         metric.filter = "all",
#         build.clamping.mask = TRUE
#       )
# 
#       # Project ensemble models
#       myBiomodEMProj <- BIOMOD_EnsembleForecasting(
#         bm.em = myBiomodEM,
#         bm.proj = myBiomodProjectionFuture,
#         proj.name = paste0(year, "_", model, "_", scenario, "_EM"),
#         models.chosen = "all",
#         metric.binary = "BOYCE",
#         metric.filter = "BOYCE",
#         on_0_1000 = FALSE,
#         nb.cpu = 8,
#         na.rm = FALSE
#       )
# 
#       myBiomodEMProj
#     }
#   }
# }

# End time
Sys.time() - tmp
