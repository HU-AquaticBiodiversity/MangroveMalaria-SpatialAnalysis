# ==============================================================================
# Script Name: Mangrove Cover Spatial Extraction 
# Description: Calculates the total surface area of mangrove forests within 
#              1-50 km radii around coastal malaria survey locations. 
#              Uses Global Mangrove Watch (GMW) high-resolution raster tiles.
#
# ### NOTE: This code should ideally be run on a High-Performance Computing 
# ### (HPC) cluster, as the spatial operations on large rasters are highly 
# ### memory- and CPU-intensive.
# ==============================================================================

##-------------------------------------------
## 1. Load Required Libraries
##-------------------------------------------
# Spatial data processing and geometry operations
library(sf)
library(sp)
library(lwgeom)
library(geosphere)
library(osmdata)
library(units)

# Raster data processing and fast extraction
library(terra)
library(exactextractr)
library(tidyterra)

# Data wrangling and tabular operations
library(dplyr)
library(tidyverse)
library(data.table)

# Parallel processing and progress tracking
library(parallel)
library(doParallel)
library(pbapply)

# Plotting and mapping
library(ggplot2)
library(rnaturalearth)
library(ggsflabel)
library(ape)

# Define the years for which mangrove data are available
years = c(1996, 2007:2010, 2015:2020)

##--------------------------------------------
## 2. Download and Unzip Mangrove Raster Data
##--------------------------------------------
# NOTE: This section is commented out. Only run this if the Global Mangrove 
# Watch (GMW) GeoTIFF files are not already present on your local machine/HPC.

# outDir = "./res/mangrove_raster/"
# dir.create(outDir)
# mangrove.files = sapply(years, function(x) paste0("gmw_v3_", x, '_gtiff.zip'))
# download_zenodo(doi = "10.5281/zenodo.6894272", path = outDir, mangrove.files)

# Unzip files. (NOTE: It might be faster to do this manually via the terminal).
# pblapply(list.files(path = outDir, pattern = "*.zip"),
#        function(i) {unzip(paste0(outDir,i), exdir=outDir, overwrite = F)})

##-----------------------------------------------
## 3. Load Vector Data (Survey Polygons & Points)
##-----------------------------------------------
# Load the pre-computed concentric buffer polygons (1 to 50 km) for all survey sites
coastal.PR.polygons = st_read("./data/polygons_pr_malaria.shp")



# Extract only the maximum radius (50 km) polygons. We will use these as bounding 
# boxes to quickly check if a survey site is near ANY mangroves before doing heavy maths.
coastal.PR.polygons.50 = coastal.PR.polygons %>%
  filter(buffer == 50)

# Load the tabular malaria prevalence survey data
coastal.PR = read.csv("./data/coastal.PR.final.csv", sep = ",", header = T)

##----------------------------------------------
## 4. Extracting Mangrove Surfaces for Each Year
##----------------------------------------------
# Loop through each available year of mangrove data
lapply(years, function(y) {
  start.time = Sys.time()
  
  # Identify all GeoTIFF raster tiles for the current year
  tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", y),
                     pattern = ".tif", full.names = T)
  
  # Load the raster tiles into a list using the 'terra' package
  tile.list = lapply(tiles, rast)
  
  # Filter survey locations that match the current raster year (y) 
  # or the subsequent year (y+1) to account for survey overlap.
  # Create a 'dummy' variable combining longitude and latitude to use as a unique ID.
  coord = coastal.PR %>% filter(year_start %in% c(y, y+1)) %>%
    select(lon, lat) %>% unite("dummy", lon, lat, sep = "_", remove = F)
  
  # Subset the 50-km polygons to only include locations surveyed in this specific year
  polygons = coastal.PR.polygons.50 %>% 
    unite("dummy", lon, lat, sep = "_", remove = F) %>%
    filter(dummy %in% coord$dummy) %>%
    select(-dummy)
  
  nrow(polygons) # Print the number of sites to be processed this year
  
  ## CONDITION: Only proceed if there are actually malaria surveys for this year
  if(nrow(polygons) > 0) {
    
    ## Calculate mangrove extent for each survey location using a progress-bar apply loop
    me.table = pblapply(1:nrow(polygons), function(p) {
      
      # Isolate the current survey location's 50-km polygon
      poly = polygons[p,]
      
      # OPTIMISATION STEP: Determine which mangrove raster tiles overlap with the 50-km polygon.
      # st_crop will fail/return empty if the polygon and tile do not intersect.
      tile.overlaps = lapply(1:length(tile.list), function(x) {
        nrow(st_crop(poly, tile.list[[x]])) != 0
      }) %>%  unlist()
      
      # Select all concentric polygons (1 km through 50 km) for this specific location
      poly.1.50 = coastal.PR.polygons %>% filter(lon == poly$lon & lat == poly$lat)
      
      ## CONDITION: Do any raster tiles overlap with the 50-km bounding box?
      ## If NO: Save time by skipping extraction and simply assigning 0 area.
      if(sum(tile.overlaps)==0) {
        mc = rep(0, 50) 
        
        ## If YES: Calculate the exact area of mangrove pixels inside each 1-50 km polygon
      } else {
        # Stitch the overlapping raster tiles together into a single mosaic
        tile.sel = mosaic(sprc(tile.list[tile.overlaps]))
        
        
        
        ## 1. Convert pixel values into actual geographic surface area (km²) using 'cellSize'.
        ##    Non-mangrove pixels are masked out (ignored).
        ## 2. 'exact_extract' calculates the precise sum of the mangrove area falling 
        ##    inside the boundaries of the 1-50 km vector polygons.
        mc = exact_extract(cellSize(tile.sel, unit = 'km', mask = T),
                           poly.1.50, fun = 'sum', progress = F)
      }
      
      # Print progress to the console
      print(paste0(p, "/", nrow(polygons)))
      
      # Attach the calculated mangrove cover (mc) vector to the polygon data
      poly.1.50 %>% mutate(mangrove.cover = mc)
      
    }) %>% bind_rows() %>%
      st_drop_geometry() # Drop the heavy spatial geometry to save as a clean data frame
    
    me.table
    
    ## EXPORT the final table for the current year as a CSV
    write.table(me.table, paste0('./data/MangroveCover/MangroveCover_', y, '.csv'), 
                sep = ",", row.names = F)
  }
  
  # Log the time taken to process the year
  end.time = Sys.time()
  time.taken = round(end.time - start.time, 2)
  print(time.taken)
  print(paste("Year", y, ": Done!", sep = " "))
})