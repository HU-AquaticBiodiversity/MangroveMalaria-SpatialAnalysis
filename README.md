# Mangrove-Malaria Project

Welcome to the repository for the Mangrove-Malaria study by the Aquatic Biodiversity Group at Hasselt University. This repository contains all R and Python code employed in the spatial and structural equation modeling (SEM) analysis of this study.

To explore our findings interactively, you can access the ShinyApp containing all path diagrams of the structural equation model (SEM) analysis by following [this link](https://github.com/HU-AquaticBiodiversity/Mangrove-Malaria_Study/tree/main/src/Mangrove-Malaria_ShinyApp).

## Installation & Setup

This project uses a combination of R (for spatial extraction, data assembly, and SEM) and Python (for Machine Learning pipelines). 

### R Dependencies
Most spatial and data wrangling operations require the following R packages. You can install them using:

```R
install.packages(c("dplyr", "tidyr", "sf", "terra", "exactextractr", "parallel", "doParallel", "ggplot2", "piecewiseSEM", "nlme", "malariaAtlas"))
```

### Python Dependencies
The machine learning optimization scripts require Python 3.8+ and the following libraries. You can install them via pip:

```bash
pip install pandas numpy scikit-learn xgboost tensorflow scikeras shap
```

## Data Description

### Scripts (`src/`)

**`01_data_processing.R`** Downloads malaria prevalence (PR) data from the Malaria Atlas Project and Demographic and Health Surveys (DHS). Filters points to those within 50 km of the African coastline and generates concentric spatial polygons (1-50 km radii) for both public and DHS datasets simultaneously.

**`02_extract_mangrove_cover.R`** Calculates the total surface area of mangrove forests within 1-50 km radii around coastal malaria survey locations using Global Mangrove Watch (GMW) high-resolution raster tiles.  
*Note: This code should ideally be run on a High-Performance Computing (HPC) cluster, as the spatial operations on large rasters are highly memory- and CPU-intensive.*

**`03_extract_ndvi.R`** For specific survey years, this script extracts the spatial extent of mangroves within 1-50 km radii of malaria survey sites. It then extracts the mean vegetation health (NDVI) specifically within those mangrove areas over time.

**`04_data_assembly.R`** Merges point prevalence data with multi-scale spatial rasters (weather, NDVI, human impact, mosquito distributions) using parallelized extraction. Prepares the final imputed dataset for SEM modeling.

**`05a_sem_functions.R`** Contains the core functions to prep data, dynamically update SEM formulas based on d-separation and p-values, and run parallelized iterations.

**`05b_sem_fitting.R`** Defines the baseline SEM formulas and dynamically generates alternative formula sets for robustness checks (spatial scales, human impact, and non-linear weather relationships) before passing them to the optimizer.

**`05c_sem_bootstrapping.R`** Generates a cluster bootstrap for a piecewise SEM across spatial scales accounting for spatial autocorrelation. Filters non-converged runs, calculates CIs, and back-transforms direct, indirect, and total effects into comparative units and Odds Ratios.

**`06a_sem_predictions.R`** Plots raw data against glmmPQL trendlines and bootstrapped SEM effects.

**`06b_sem_interpretation.R`** Generates heatmaps, structural equation model (SEM) path diagrams, effect size bar charts, and spatial maps for the analysis. Evaluates model robustness across spatial resolutions.

**`07a_ml_data_prep.R`** Cleans, engineers, and transforms an imputed malaria dataset into a structured format ready for Machine Learning algorithms. Converts aggregated prevalence survey data into binary target classes (infected/uninfected) with case weights.

**`07b_ml_pipeline.py`** Trains and evaluates Logistic Regression, XGBoost, and a Keras Neural Network using GridSearchCV. Includes data preprocessing, feature scaling, and performance evaluation via ROC AUC.

**`08_plot_spatial_buffers.R`** Generates a high-resolution, 3-panel figure illustrating mangrove land cover, mangrove NDVI (vegetation health), and population density within specific spatial buffer zones (1-50 km) around a selected site in the Saloum Delta, Senegal.

---

### Data (`data/`)

**`country.table.csv`** Metadata for country selection including ISO codes and appropriate coordinate reference system (CRS).

**`Coastal.PR.final.csv`** Malaria prevalence data from MalariaAtlas for 28 coastal countries with data points being 50 km or less off the coastline.

> **⚠️ Note on Health Data Access:** > Malaria infections were accessed via the MalariaAtlas project. For a large portion of these data, access is restricted and needs to be requested via the [Demographic and Health Surveys (DHS)](https://www.dhsprogram.com/) programme of USAID. To download the DHS data through the `malariaAtlas` R package, you must first create an account and request access on their portal.
