# ==============================================================================
# Script Name: Machine Learning Input Preparation
# Description: Cleans, engineers, and transforms an imputed malaria dataset 
#              into a structured format ready for Machine Learning algorithms.
#              Specifically, it converts aggregated prevalence survey data 
#              into binary target classes (infected/uninfected) with case weights.
# ==============================================================================

## =============================================================================
## 1. LOAD NECESSARY LIBRARIES
## =============================================================================
library(dplyr)          # For core data manipulation and piping (%>%)
library(tidyr)          # For reshaping data (pivot_longer)

## =============================================================================
## 2. DATA LOADING & INITIAL FEATURE ENGINEERING
## =============================================================================

# Load the previously imputed dataset
me_impute_filter = read.csv("./data/alldata.impute.csv") %>%
  # Standardise column names for consistency
  rename(mean_ndvi = mean.ndvi, 
         mangrove_cover = mangrove.cover) %>%
  # Engineer new variables based on the survey data
  mutate(
    # Reconstruct the absolute number of infected individuals.
    # pr (prevalence rate) * examined (total tested) = infected cases.
    # We round it to ensure we have clean integer counts for individuals.
    infected = round(pr * examined, 0),
    
    # Recalculate prevalence as a sanity check (though it gets dropped later)
    pr.new = infected / examined,
    
    # Parse the start year and Julian day/month into a proper R Date object
    date = as.Date(paste(year_start, month_start), '%Y%j')
  )

## =============================================================================
## 3. BINARY TARGET TRANSFORMATION & FEATURE SELECTION
## =============================================================================

# Transform the data from an aggregated summary format into a binary 
# classification format suitable for Machine Learning models.


input = me_impute_filter %>% 
  # Calculate the absolute number of uninfected individuals
  mutate(uninfected = examined - infected) %>%
  
  # Drop the proportion/summary columns, as we are shifting to absolute counts
  dplyr::select(-pr, -pr.new, -examined) %>%
  
  # RESHAPE DATA: 
  # This splits a single survey row into two distinct rows: 
  # One for the 'infected' count and one for the 'uninfected' count.
  # The counts are stored in a new column called 'weight'.
  pivot_longer(cols = c("infected", "uninfected"),
               names_to = "infected", 
               values_to = "weight") %>%
  
  # Convert the text label back into a strictly binary numeric target (1 or 0).
  # 1 = Infected, 0 = Uninfected.
  mutate(infected = ifelse(infected == "infected", 1, 0)) %>%
  
  # Remove rows where the weight is 0. 
  # (e.g., if a survey had 0 infections, we don't need a row saying 'infected=1, weight=0')
  filter(weight != 0) %>% 
  
  # FINAL FEATURE SELECTION: 
  # Keep only the target variable, the predictors (environmental/spatial/climatic), 
  # and the row 'weight' for the ML algorithm.
  dplyr::select(
    # Target Variable
    infected, 
    
    # Mangrove & Land Cover Predictors
    mean_ndvi, mangrove_cover, mangrove.cover.min1, coastline.dist, lc_crop, 
    
    # Population & Weather Predictors
    pop_dens_median, anomaly_2t, anomaly_tp, mean_2t, mean_tp, 
    anomaly_2t_6m, anomaly_tp_6m, mean_2t_6m, mean_tp_6m,
    
    # Spatiotemporal & Methodological Control Predictors
    lat, lon, year_start, month_start, method, 
    
    # Algorithm Case Weight
    weight
  )

## =============================================================================
## 4. EXPORT MACHINE LEARNING DATASET
## =============================================================================

# Export the fully processed dataset to a CSV file.
# Note: 'row.names = F' ensures we don't export an unnecessary index column.
write.table(input, "./data/ML_input_file.csv", row.names = FALSE, sep = ",")