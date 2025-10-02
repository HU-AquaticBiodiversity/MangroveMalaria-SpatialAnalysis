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

<i>NOTE: Access to health data (malaria infections) needs to be requested via the <a href="https://www.dhsprogram.com">Demographic and Health Surveys (DHS)</a> programme of USAID via </br>
Therefore, the geographical coordinates of all the dataset below have been removed.</i>
</br></br>
<b>alldata.impute.csv</b> - .</br></br>
<b>annual.raster.data.csv</b> - Values for raster layers, which are available for a multiple years incl. health management data and population density.</br></br>
<b>const.raster.data.csv</b> - Values for raster layers, which are available for a single year.</br></br>
<b>Coastal.PR.dhs.final.csv</b> - Malaria prevalence data from DHS for 28 coastal countries with data point being 50 km or less off the coastline.</br></br>
<b>cropland.csv</b> - .</br></br>
<b>Polygons_pr_malaria.shp</b> - Shape for buffer zones of 1 to 50 km around each unique (coastal) prevalence data location from the Coastal.PR.final.csv dataset. CRS: WGS84.</br></br>
<b>ML_input_file.csv</b> - .</br></br>
<b>ndvi_data.csv</b> - NDVI within mangroves for each shape in Polygons_pr_malaria.shp.</br></br>
<b>nHDI_data.csv</b> - .</br></br>

