# ==============================================================================
# SEM MODEL CONFIGURATION & SENSITIVITY ANALYSIS SETTINGS
# Description: Defines the baseline SEM formulas and dynamically generates 
# alternative formula sets for robustness checks (spatial scales, human impact, 
# and non-linear weather relationships) before passing them to the optimizer.
# ==============================================================================

library(dplyr)
library(parallel)
library(piecewiseSEM)
library(DHARMa)
library(stringr)
library(nlme)
library(performance)
library(MASS)
library(pbapply)
library(doParallel)
library(DiagrammeR)
library(tidyr)

# Load custom functions (contains data.prep and model.run)
source("./src/functions.model_fitting.R")

# --- 1. System & Data Preparation ---

# Detect the number of cores for parallel processing downstream
numberOfCores <- 25

# Load and prep the master dataset (imputed)
all.data <- read.csv("./data/alldata.impute.csv") %>% 
  data.prep()

# Quick sanity check of the main data
summary(all.data)

# Load and prep the smaller dataset (used to test imputation robustness)
all.data.small <- read.csv("./data/alldata.small.csv") %>% 
  data.prep()


# --- 2. Baseline Model Specification (Hypothesis 1) ---

# Define the starting formulas for the initially hypothesized model structure.
# (Formatted with line breaks for readability)
start.formulas_1 <- list(
  
  # A. Prevalence Model
  pr.new ~ mean_ndvi + mangrove_cover + coastline.dist + 
    lc_crop + pop_dens_median + 
    anomaly_2t + anomaly_tp + 
    mean_2t + mean_tp + 
    anomaly_2t_6m + anomaly_tp_6m + 
    mean_2t_6m + mean_tp_6m,
  
  # B. NDVI Model
  mean_ndvi ~ mangrove_cover + pop_dens_median + 
    anomaly_2t + anomaly_tp + 
    mean_2t + mean_tp + 
    anomaly_2t_6m + anomaly_tp_6m + 
    mean_2t_6m + mean_tp_6m,
  
  # C. Mangrove Cover Model
  mangrove_cover ~ mangrove.cover.min1 + lc_crop + pop_dens_median + coastline.dist + 
    anomaly_2t + anomaly_tp + 
    mean_2t + mean_tp + 
    anomaly_2t_6m + anomaly_tp_6m + 
    mean_2t_6m + mean_tp_6m
)

# Identify variables that will be swapped out during spatial sensitivity checks
fixed.var <- c("lc_crop", "pop_dens_median", "anomaly_2t", "anomaly_tp",
               "mean_2t", "mean_tp", "anomaly_2t_6m", "anomaly_tp_6m",
               "mean_2t_6m", "mean_tp_6m")


# --- 3. Alternative Formulas for Robustness Checks ---

# Model Set 2: Small Dataset (Formulas remain identical, data swaps in model.run)
start.formulas_1_small <- start.formulas_1



# Model Set 3: 5km Spatial Radius
# Dynamically removes the base 'fixed.var' and adds the '_5' versions
start.formulas_1_5k <- lapply(1:3, function(x) {
  # Constructs a string like: "~ . - lc_crop ... + lc_crop_5 ..."
  update_string <- paste0("~ . - ", paste(fixed.var, collapse = " - "),
                          " + ", paste(str_c(fixed.var, "_5"), collapse = " + "))
  update(start.formulas_1[[x]], as.formula(update_string))
})

# Model Set 4: 20km Spatial Radius
# Dynamically removes the base 'fixed.var' and adds the '_20' versions
start.formulas_1_20k <- lapply(1:3, function(x) {
  update_string <- paste0("~ . - ", paste(fixed.var, collapse = " - "),
                          " + ", paste(str_c(fixed.var, "_20"), collapse = " + "))
  update(start.formulas_1[[x]], as.formula(update_string))
})

# Model Set 5: Additional Malaria Variables
# Adds healthcare and ITN access strictly to the prevalence formula (Index 1)
start.formulas_1_newvars <- list(
  update(start.formulas_1[[1]], ~ . + t_healthcare_motor.median + ITN_access_mean_median),
  start.formulas_1[[2]],
  start.formulas_1[[3]]
)

# Model Set 6: Non-Linear (Quadratic) Weather Effects
# Appends squared weather terms to ALL three sub-models to test for non-linear relationships
weather.formula <- "~ . + anomaly_2t_sqr + anomaly_tp_sqr + mean_2t_sqr + mean_tp_sqr + 
                    anomaly_2t_6m_sqr + anomaly_tp_6m_sqr + mean_2t_6m_sqr + mean_tp_6m_sqr"

start.formulas_1_sqr <- lapply(1:3, function(s) {
  update(start.formulas_1[[s]], as.formula(weather.formula))
})


# --- 4. Model Execution & Optimization ---

## NOTE: The 'model.run' function handles the parallel execution.
## If using a Linux HPC, it likely uses mclapply. If Windows, it should use parLapply.

## Optimisation steps executed inside model.run():
##  - 1: Default baseline model
##  - 2 "_small": Reduced data to test robustness of imputation
##  - 3 "_5k": Fixed variables calculated at alternative 5-km radius
##  - 4 "_20k": Fixed variables calculated at alternative 20-km radius
##  - 5 "_newvars": Malaria-relevant variables added to the prevalence model
##  - 6 "_sqr": Quadratic functions for weather variables

# Run stepwise optimization across all 6 model variations simultaneously
model.run(step = 1, model.select = 1:6)
model.run(step = 2, model.select = 1:6)
model.run(step = 3, model.select = 1:6)
model.run(step = 4, model.select = 1:6)
model.run(step = 5, model.select = 1:6)
model.run(step = 6, model.select = 1:6)