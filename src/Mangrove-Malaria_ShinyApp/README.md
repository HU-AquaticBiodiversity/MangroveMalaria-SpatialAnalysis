# Mangrove-Malaria Spatial Resolution Explorer (Shiny App)

This Shiny application provides an interactive environment to explore the effect of spatial resolution on the **Mangrove-Malaria relationship** across coastal Africa from Cruz-Laufer et al. (in review).

By dynamically adjusting the spatial scale, users can visualise how different environmental and geographical factors influence local malaria prevalence.

## Overview

Each path diagram generated in this app represents a **piecewise Structural Equation Model (SEM)**. The 50 models shown here correspond to the main models from Cruz-Laufer et al (in review), with missing data having been imputed and human impact variables (population density, agricultural land cover) and weather variables (mean and anomaly of temperature and precipitation during and 6 months prior to malaria survey period) having been calculated at a fixed radius of 10 km (see the preprint https://ecoevorxiv.org/repository/view/10430/ for more details).

The spatial resolution ($r$) corresponds to the radius (in kilometres) surrounding a malaria survey site. Within this radius, key mangrove variables — specifically **mangrove land cover** and **mangrove NDVI (greenness)** — were calculated. 

### How to interpret the visualisations:
* **Nodes:** Represent the variables in the model (e.g., Malaria Prevalence, Mangrove Cover, Weather anomalies).
* **Paths (Arrows):** Represent the relationships between variables.
* **Path Thickness:** Indicates the **robustness** of the relationship (i.e., the proportion of models that statistically support the relationship).
* **Path Colour:** Indicates the direction of the standardised effect:
    * **Blue:** Positive association
    * **Red:** Negative association
    * **Grey:** Unsupported / non-significant relationship

## Prerequisites

Before running the application, ensure you have **R** (and optionally RStudio) installed on your machine. The application relies on the following R packages:

* `shiny`
* `DiagrammeR`
* `dplyr`
* `tidyr`

*(Note: The app script is designed to automatically detect and install missing packages, but it is good practice to ensure they are available).*

## How to run the Shiny App

**Step 1: Download the repository**
Clone this GitHub repository to your local machine using Git, or download it as a ZIP file and extract it.

**Step 2: Launch the App in R**
Open R or RStudio, and run the script below.

```text
# 1. In your terminal or command prompt:
git clone [https://github.com/HU-AquaticBiodiversity/MangroveMalaria-SpatialAnalysis.git](https://github.com/HU-AquaticBiodiversity/MangroveMalaria-SpatialAnalysis.git)

# 2. In R or RStudio (ensure your working directory is set to where you cloned the repo):
library(shiny)
runApp("MangroveMalaria-SpatialAnalysis/src/Mangrove-Malaria_ShinyApp")
