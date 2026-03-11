# ==============================================================================
# Script Name: Malaria Prevalence Data Processing & Spatial Buffering
# Description: Downloads malaria prevalence (PR) data from the Malaria Atlas 
#              Project and Demographic and Health Surveys (DHS). Filters points 
#              to those within 50 km of the African coastline and generates 
#              concentric spatial polygons (1-50 km radii) for both public 
#              and DHS datasets simultaneously.
# ==============================================================================

##-------------------------------------------
## 1. Load Required Libraries
##-------------------------------------------
# Malaria data retrieval
library(malariaAtlas)
library(rdhs)

# Spatial data processing and mapping
library(sf)
library(sp)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(osmdata)

# Data wrangling and general utilities
library(dplyr)
library(tidyverse)
library(zen4R)
library(pbapply)

# Visualisation
library(ggplot2)
library(viridis)

##----------------------------------------------------------
## 2. Create Low-Resolution Shapefile for Mangrove Scenes
##----------------------------------------------------------
# This section creates simplified, low-resolution boundary boxes (polygons) of 
# the high-resolution mangrove raster tiles. This is necessary to locate and 
# request specific scenes from the USGS Earth Explorer without exceeding 
# download file size limits.

# Load names of individual mangrove raster tiles 
tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_mng_union/"),
                   pattern = ".tif", full.names = T)

# Load the raster information using the 'terra' package
tile.list = pblapply(tiles, rast)

# Aggregate (downsample) the raster heavily to simplify the shapefile
# Using fact=4500 means combining 4500x4500 pixels into a single mean pixel
tiles.aggr = pblapply(tile.list, function(x) aggregate(x, fact = 4500, fun = 'mean', na.rm = T))

# Split rasters into three distinct geographical regions (West, Northeast, Southeast)
# This meets the strict file size requirements for the NASA/USGS download portals.

tiles.mosaic.west = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(-20, 20, -37, 30))
tiles.mosaic.northeast = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(20, 55, 5, 30))
tiles.mosaic.southeast = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(20, 55, -37, 5))

# Convert the raster mosaics into vector polygons (shapefiles)
mangrove.shp.west = st_as_sf(as.polygons(tiles.mosaic.west))
mangrove.shp.northeast = st_as_sf(as.polygons(tiles.mosaic.northeast))
mangrove.shp.southeast = st_as_sf(as.polygons(tiles.mosaic.southeast))

# Export the region shapefiles
write_sf(mangrove.shp.west, dsn = './data/mangrove_vcts_aggr/mangrove_vec_west.shp', append = F)
write_sf(mangrove.shp.northeast, dsn = './data/mangrove_vcts_aggr/mangrove_vec.northeast.shp', append = F)
write_sf(mangrove.shp.southeast, dsn = './data/mangrove_vcts_aggr/mangrove_vec.southeast.shp', append = F)

# Code snippet for quick visual check of the raster boundaries
# ggplot() +  
#   geom_sf(data = mangrove.shp.southeast) +
#   scale_fill_viridis()

##----------------------------------------------------------
## 3. Download and Clean Malaria Prevalence Data
##----------------------------------------------------------

# Load reference list of countries relevant to the study
country.list = read.csv("./data/country_table.csv")

# Create an expanded list of target years. 
# We include the year immediately following our target years (year + 1) to 
# capture surveys that span across calendar years.
years.expanded = c(1996:1997, 2007:2011, 2015:2021)

# Download publicly available malaria prevalence (PR) data from the Malaria Atlas Project
malaria.PR.raw = getPR(continent = "Africa", species = "Pf") %>%
  # Filter strictly for countries included in this study
  filter(country_id %in% country.list$ISO3,
         year_start %in% years.expanded | year_end %in% years.expanded)

# Request and append restricted coordinate data from the Demographic and Health Surveys (DHS)
# Note: Data access must be requested and approved through the official DHS website.
malaria.PR.raw.dhs = fillDHSCoordinates(malaria.PR.raw %>% filter(is.na(latitude)),
                                        email = "cruzmamo@googlemail.com",
                                        project = "Malaria prevalence in mangrove forests",
                                        timeout = 60)

# Clean the public dataset (remove rows missing coordinates or PR values)
malaria.PR = malaria.PR.raw %>% 
  filter(!is.na(latitude)) %>%
  filter(!is.na(pr)) %>%
  as.data.frame()

# Clean the restricted DHS dataset and standardise the ID column
malaria.PR.dhs = malaria.PR.raw.dhs %>% 
  filter(!is.na(latitude)) %>%
  filter(!is.na(pr)) %>%
  as.data.frame() %>%
  mutate(site_id = dhs_id)

# Compare data volumes in the console
print(paste("Public PR records:", nrow(malaria.PR)))
print(paste("DHS + Public PR records:", nrow(malaria.PR.dhs)))

##----------------------------------------------------------
## 4 & 5. Filter Coastal Data & Generate Buffer Polygons
##        (Automatically processes both Public & DHS data)
##----------------------------------------------------------

# Load high-resolution global coastline shape data from Natural Earth
coastline = ne_coastline(scale = 10, returnclass = "sf")

# Place both datasets into a named list to process them iteratively
malaria_datasets = list(
  "public" = malaria.PR,
  "dhs" = malaria.PR.dhs
  # ==============================================================================
  # Script Name: Malaria Prevalence Data Processing & Spatial Buffering
  # Description: Downloads malaria prevalence (PR) data from the Malaria Atlas 
  #              Project and Demographic and Health Surveys (DHS). Filters points 
  #              to those within 50 km of the African coastline and generates 
  #              concentric spatial polygons (1-50 km radii) for both public 
  #              and DHS datasets simultaneously.
  # ==============================================================================
  
  ##-------------------------------------------
  ## 1. Load Required Libraries
  ##-------------------------------------------
  # Malaria data retrieval
  library(malariaAtlas)
  library(rdhs)
  
  # Spatial data processing and mapping
  library(sf)
  library(sp)
  library(terra)
  library(tidyterra)
  library(rnaturalearth)
  library(osmdata)
  
  # Data wrangling and general utilities
  library(dplyr)
  library(tidyverse)
  library(zen4R)
  library(pbapply)
  
  # Visualisation
  library(ggplot2)
  library(viridis)
  
  ##----------------------------------------------------------
  ## 2. Create Low-Resolution Shapefile for Mangrove Scenes
  ##----------------------------------------------------------
  # This section creates simplified, low-resolution boundary boxes (polygons) of 
  # the high-resolution mangrove raster tiles. This is necessary to locate and 
  # request specific scenes from the USGS Earth Explorer without exceeding 
  # download file size limits.
  
  # Load names of individual mangrove raster tiles 
  tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_mng_union/"),
                     pattern = ".tif", full.names = T)
  
  # Load the raster information using the 'terra' package
  tile.list = pblapply(tiles, rast)
  
  # Aggregate (downsample) the raster heavily to simplify the shapefile
  # Using fact=4500 means combining 4500x4500 pixels into a single mean pixel
  tiles.aggr = pblapply(tile.list, function(x) aggregate(x, fact = 4500, fun = 'mean', na.rm = T))
  
  # Split rasters into three distinct geographical regions (West, Northeast, Southeast)
  # This meets the strict file size requirements for the NASA/USGS download portals.
  
  tiles.mosaic.west = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(-20, 20, -37, 30))
  tiles.mosaic.northeast = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(20, 55, 5, 30))
  tiles.mosaic.southeast = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(20, 55, -37, 5))
  
  # Convert the raster mosaics into vector polygons (shapefiles)
  mangrove.shp.west = st_as_sf(as.polygons(tiles.mosaic.west))
  mangrove.shp.northeast = st_as_sf(as.polygons(tiles.mosaic.northeast))
  mangrove.shp.southeast = st_as_sf(as.polygons(tiles.mosaic.southeast))
  
  # Export the region shapefiles
  write_sf(mangrove.shp.west, dsn = './data/mangrove_vcts_aggr/mangrove_vec_west.shp', append = F)
  write_sf(mangrove.shp.northeast, dsn = './data/mangrove_vcts_aggr/mangrove_vec.northeast.shp', append = F)
  write_sf(mangrove.shp.southeast, dsn = './data/mangrove_vcts_aggr/mangrove_vec.southeast.shp', append = F)
  
  # Code snippet for quick visual check of the raster boundaries
  # ggplot() +  
  #   geom_sf(data = mangrove.shp.southeast) +
  #   scale_fill_viridis()
  
  ##----------------------------------------------------------
  ## 3. Download and Clean Malaria Prevalence Data
  ##----------------------------------------------------------
  
  # Load reference list of countries relevant to the study
  country.list = read.csv("./data/country_table.csv")
  
  # Create an expanded list of target years. 
  # We include the year immediately following our target years (year + 1) to 
  # capture surveys that span across calendar years.
  years.expanded = c(1996:1997, 2007:2011, 2015:2021)
  
  # Download publicly available malaria prevalence (PR) data from the Malaria Atlas Project
  malaria.PR.raw = getPR(continent = "Africa", species = "Pf") %>%
    # Filter strictly for countries included in this study
    filter(country_id %in% country.list$ISO3,
           year_start %in% years.expanded | year_end %in% years.expanded)
  
  # Request and append restricted coordinate data from the Demographic and Health Surveys (DHS)
  # Note: Data access must be requested and approved through the official DHS website.
  malaria.PR.raw.dhs = fillDHSCoordinates(malaria.PR.raw %>% filter(is.na(latitude)),
                                          email = "cruzmamo@googlemail.com",
                                          project = "Malaria prevalence in mangrove forests",
                                          timeout = 60)
  
  # Clean the public dataset (remove rows missing coordinates or PR values)
  malaria.PR = malaria.PR.raw %>% 
    filter(!is.na(latitude)) %>%
    filter(!is.na(pr)) %>%
    as.data.frame()
  
  # Clean the restricted DHS dataset and standardise the ID column
  malaria.PR.dhs = malaria.PR.raw.dhs %>% 
    filter(!is.na(latitude)) %>%
    filter(!is.na(pr)) %>%
    as.data.frame() %>%
    mutate(site_id = dhs_id)
  
  # Compare data volumes in the console
  print(paste("Public PR records:", nrow(malaria.PR)))
  print(paste("DHS + Public PR records:", nrow(malaria.PR.dhs)))
  
  ##----------------------------------------------------------
  ## 4 & 5. Filter Coastal Data & Generate Buffer Polygons
  ##        (Automatically processes both Public & DHS data)
  ##----------------------------------------------------------
  
  # Load high-resolution global coastline shape data from Natural Earth
  coastline = ne_coastline(scale = 10, returnclass = "sf")
  
  # Place both datasets into a named list to process them iteratively
  malaria_datasets = list(
    "public" = malaria.PR,
    "dhs" = malaria.PR.dhs
  )
  
  # Use lapply to run the entire spatial filtering and buffering pipeline for both datasets
  lapply(names(malaria_datasets), function(data_type) {
    
    print(paste("=========================================="))
    print(paste("Processing spatial buffers for dataset:", toupper(data_type)))
    print(paste("=========================================="))
    
    # Extract the current dataset from the list
    PR = malaria_datasets[[data_type]]
    
    # Isolate only the countries for which we actually have prevalence data in this dataset
    countries.available = country.list %>% filter(ISO3 %in% PR$country_id)
    
    # Filter dataset to only include survey points located up to 50 km inland from the coast
    coastal.PR = pblapply(1:nrow(countries.available), function(c) {
      
      # Extract data for the specific country and convert to spatial (sf) object
      disease.data = PR %>%
        filter(country_id == countries.available$ISO3[c]) %>%
        as.data.frame() %>%
        st_as_sf(coords = c('longitude', 'latitude'), 
                 crs = st_crs(4326), dim = "XY") %>%
        st_transform(countries.available$CRS[c]) 
      
      # Transform the coastline shapefile to match the country's projected CRS
      coastline.crs = coastline %>%
        st_transform(countries.available$CRS[c])
      
      # Filter out any invalid geometries from the coastline
      coastline.filtered = coastline.crs[ifelse(is.na(st_is_valid(coastline.crs)), F, T),]
      
      # Calculate the distance matrix between all survey points and all coastline segments
      coastline.alldist = st_distance(coastline.filtered, disease.data)
      
      # Extract the absolute minimum distance (in kilometres) for each individual survey point
      coastline.dist =  sapply(1:nrow(disease.data), function(k) {          
        min(coastline.alldist[,k], na.rm = T)/1000
      })
      
      # Attach distances and filter out points further than 50 km inland
      disease.filtered = disease.data %>%
        mutate(coastline.dist = coastline.dist) %>%
        filter(coastline.dist <= 50) 
    })
    
    # Identify which countries contain viable coastal PR data after filtering
    data.available = sapply(1:nrow(countries.available), function (c) {
      nrow(coastal.PR[[c]]) != 0
    })
    
    # Merge processed country datasets back together
    coastal.PR.final = bind_rows(
      lapply(coastal.PR, function(x) {
        x %>% 
          st_transform(crs = 4326) %>% # Revert back to WGS84 for tabular output
          dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                        lat = sf::st_coordinates(.)[,2]) %>%
          st_drop_geometry() 
      })
    )
    
    # Dynamically assign the CSV file name based on the dataset type
    csv_file_name = ifelse(data_type == "dhs", 
                           "./data/coastal.PR.dhs.final.csv", 
                           "./data/coastal.PR.final.csv")
    
    # Export tabular coordinate data
    write.table(coastal.PR.final, file = csv_file_name, sep = ",", row.names = F)
    print(paste("Saved CSV:", csv_file_name))
    
    # Generate 1-50 km concentric buffer polygons for downstream raster extraction
    
    coastal.PR.polygons = bind_rows(pblapply(1:50, function(y) {
      bind_rows(
        lapply(coastal.PR[data.available], function(x) {
          
          # Extract raw longitude and latitude points
          points = x %>%
            dplyr::select(geometry) %>%
            st_transform(crs = 4326) %>%
            dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                          lat = sf::st_coordinates(.)[,2]) %>%
            st_drop_geometry()
          
          # Buffer the point by 'y' kilometres (y * 1000 metres), attach coords & radius label
          x %>%
            st_buffer(y * 1000) %>%
            dplyr::select(geometry) %>%
            st_transform(crs = 4326) %>%
            bind_cols(points) %>%
            unique() %>%
            mutate(buffer = y)
        })
      )
    }))
    
    # Dynamically assign the Shapefile name based on the dataset type
    shp_file_name = ifelse(data_type == "dhs", 
                           "./data/polygons_pr_malaria_dhs.shp", 
                           "./data/polygons_pr_malaria.shp")
    
    # Export the concentric buffer polygons as a shapefile
    st_write(coastal.PR.polygons, shp_file_name, append = F)
    print(paste("Saved Shapefile:", shp_file_name))
    
    # Print the final count of 50 km polygons for this dataset
    print(paste("Total 50km polygons generated:", nrow(coastal.PR.polygons %>% filter(buffer == 50))))
    
  })