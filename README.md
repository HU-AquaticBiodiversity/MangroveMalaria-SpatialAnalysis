# Mangrove Malaria project

Welcome to the repository for the Mangrove-Malaria study of Aquatic Biodiversity Group at Hasselt University. This repo contain all R and python codes employed in the study, but still remains under construction.

To access the ShinyApp with all path diagrams of the structural equation model (SEM) analysis please follow <a href="https://github.com/HU-AquaticBiodiversity/Mangrove-Malaria_Study/tree/main/src/Mangrove-Malaria_ShinyApp">this link</a>.

## Data description
### Codes (src)
<b>Data_processing.R</b> - Data download data from repositories. Currently, it only works for the mangrove polygon data (which I am not using anymore) and the prevalence data. Please be awar that in order to do download the DHS data through the MalariaAtlas package you need to create an account and request access.</br></br>
<b>Data_extraction_ME.R</b> - Calculates mangrove land cover in a radius from 1 to 50 km for each disease prevalence data point. The code uses a parallel version of the lapply (mclapply) function to speed up calculations, and should be performed on an high-perfomance computing cluster (HPC).</br></br>
<b>NDVI_extraction.R</b> - Calculates mangrove NDVI cover in a radius from 1 to 50 km for each disease prevalence data point. The code uses a parallel version of the lapply (mclapply) function to speed up calculations, and should be performed on an high-perfomance computing cluster (HPC).</br></br>
<b>Data_assembly.R</b> - Extracts data from all other data layers and table and merges them into a single file together with the mangrove cover data.</br></br>
<b>model_fitting.R</b> - Codes for fitting piecewise structural equation models and model optimisation steps.</br></br>
<b>model_interpretation.R</b> - Codes for analysing model output and producing figures.</br></br>
<b>ML_DataPrep.R</b> - Transformations to prepare data for machine learning analyses.</br></br>
<b>ML_Hyperparameters.py</b> - Optimisation of machine learning models.</br></br>
<b>Mangrove-Malaria</b> - Optimisation of machine learning models.</br></br>

### Data (data)
<b>country.table.csv</b> - Metadata for country selection including ISO codes and appropriate coordinate reference system (CRS).</br>
<b>Coastal.PR.final.csv</b> - Malaria prevalence data from MalariaAtlas for 28 coastal countries with data point being 50 km or less off the coastline.</br></br>

<i>NOTE: Health data (malaria infections) were accessed via the MalariaAtlas project. For a large portion of these data, access is restricted and needs to be requested via the <a href="https://www.dhsprogram.com">Demographic and Health Surveys (DHS)</a> programme of USAID.</br>
