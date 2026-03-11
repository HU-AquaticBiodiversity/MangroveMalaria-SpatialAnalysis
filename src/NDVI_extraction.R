# ==============================================================================
# Script Name: Mangrove Extent and NDVI Extraction pipeline
# Description: For specific survey years, this script extracts the spatial 
#              extent of mangroves within 1-50 km radii of malaria survey 
#              sites. It then extracts the mean vegetation health (NDVI) 
#              specifically within those mangrove areas over time.
# ==============================================================================

## =============================================================================
## 1. LOAD NECESSARY LIBRARIES & SET PARAMETERS
## =============================================================================
# Vector & Raster Spatial processing
library(sf)
library(sp)
library(lwgeom)
library(terra)
library(exactextractr)

# Data manipulation & utilities
library(dplyr)
library(tidyverse)
library(data.table)
library(units)
library(pbapply)
library(geosphere)

# Mapping & Visualization (Though mostly used for backend data prep here)
library(rnaturalearth)
library(ggplot2)
library(viridis)
library(osmdata)

# Parallel computing for faster raster iteration
library(parallel)
library(doParallel)

# Define the survey years to process and detect available CPU cores
years = c(2007:2010, 2015:2020)
numberOfCores = detectCores()

## =============================================================================
## 2. LOADING VECTOR DATA (POLYGONS & POINTS)
## =============================================================================

# Load the buffer polygons (1-50km rings) and combine DHS and public datasets
coastal.PR.polygons = st_read("./data/polygons_pr_malaria_dhs.shp") %>%
  bind_rows(st_read("./data/polygons_pr_malaria.shp"))

# Isolate just the maximum 50km extent. We use this as a "bounding box" to 
# quickly check if a site is near mangroves before processing all 50 inner rings.
coastal.PR.polygons.50 = coastal.PR.polygons %>%
  filter(buffer == 50)

# Load the point prevalence data and combine DHS and public datasets
coastal.PR = read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = T) %>%
  bind_rows(read.csv("./data/coastal.PR.final.csv", sep = ",", header = T) %>%
              mutate(site_id = as.factor(site_id)))



## =============================================================================
## 3. MAIN LOOP: EXTRACTING MANGROVE SURFACES PER YEAR
## =============================================================================

lapply(years, function(y) {
  start.time = Sys.time()
  
  # Filter survey coordinates that match the current loop's year (y)
  # Create a 'dummy' unique ID combining longitude and latitude
  coord = coastal.PR %>% filter(year_start == y) %>%
    select(lon, lat) %>% unite("dummy", lon, lat, sep = "_", remove = F)
  
  # Subset the 50km bounding polygons to only those needed for this specific year
  polygons = coastal.PR.polygons.50 %>% 
    unite("dummy", lon, lat, sep = "_", remove = F) %>%
    filter(dummy %in% coord$dummy) %>%
    select(-dummy)
  
  # CONDITION: Only proceed if there are actually surveys in this year
  if(nrow(polygons) > 0) {
    
    # Load file paths for Global Mangrove Watch (GMW) raster tiles for year 'y'
    tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", y),
                       pattern = ".tif", full.names = T)
    
    # Load the raster tile metadata into memory using terra
    tile.list = lapply(tiles, rast)
    
    ## -------------------------------------------------------------------------
    ## 3A. PARALLEL PROCESSING: ISOLATING MANGROVE EXTENTS
    ## -------------------------------------------------------------------------
    # Iterate through every 50km polygon to find overlapping mangroves
    me.table = mclapply(1:nrow(polygons), function(p) {
      
      # Print progress tracker to console
      print(paste0(p, "/", nrow(polygons)))
      
      # Select the current 50-km polygon
      poly = polygons[p,]
      
      # Optimization: Check which mangrove tiles intersect with this 50-km polygon
      # This prevents the script from loading raster data for the entire continent.
      tile.overlaps = lapply(1:length(tile.list), function(x) {
        nrow(st_crop(poly, tile.list[[x]])) != 0
      }) %>%  unlist()
      
      # CONDITION: Do any raster tiles overlap with this 50-km area?
      if(sum(tile.overlaps)!=0) {
        
        # Mosaic (merge) only the intersecting tiles together
        tile.sel = mosaic(sprc(tile.list[tile.overlaps]))
        
        # Pull all 50 individual concentric rings (1km, 2km... 50km) for this site.
        # Use exact_extract to quickly flag which specific rings actually touch mangroves.
        poly.1.50 = coastal.PR.polygons %>% filter(lon == poly$lon & lat == poly$lat) %>%
          mutate(overlap = exact_extract(tile.sel, ., fun = 'sum', progress = F) > 0) %>%
          filter(overlap)
        
        # CONDITION: Do any of the specific inner buffer rings touch mangroves?
        if(nrow(poly.1.50) > 0) {
          
          # Crop the mangrove raster precisely to the shape of each overlapping ring
          # and convert that raster extent into a spatial vector polygon
          tiles.crop = lapply(1:nrow(poly.1.50), function(x) {
            as.polygons(crop(tile.sel, poly.1.50[x,], mask = T)) %>%
              st_as_sf() %>%
              # Attach spatial metadata back to the new mangrove polygon
              mutate(buffer = poly.1.50[x,]$buffer,
                     lat = poly.1.50[x,]$lat,
                     lon = poly.1.50[x,]$lon)
          }) %>% .[!is.na(.)] %>% bind_rows() %>%
            # Clean up default naming from the GMW raster
            select(-starts_with("GMW"))
          
          tiles.crop # Return the finalized mangrove geometry for this site
          
        } else {NA}
      } else {NA}
    }, mc.cores = numberOfCores) %>% .[!is.na(.)] %>%
      bind_rows() # Combine parallel lists into one large spatial dataframe
    
    
    
    ## -------------------------------------------------------------------------
    ## 3B. NDVI (VEGETATION HEALTH) EXTRACTION
    ## -------------------------------------------------------------------------
    # CONDITION: If we successfully isolated mangrove extents, extract their NDVI
    if(nrow(me.table) > 0) {
      
      # Locate the NDVI raster time-series files for this year across three regions
      ndvi.files.w = list.files(
        "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/West",
        full.names = T, pattern = paste0("*NDVI_doy", y,".*tif"))
      
      ndvi.files.ne = list.files(
        "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/NorthEast",
        full.names = T, pattern = paste0("*NDVI_doy", y,".*tif"))
      
      ndvi.files.se = list.files(
        "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/SouthEast",
        full.names = T, pattern = paste0("*NDVI_doy", y,".*tif"))
      
      # Load the raster stacks
      NDVI.w = rast(ndvi.files.w)
      NDVI.ne = rast(ndvi.files.ne)
      NDVI.se = rast(ndvi.files.se)
      
      # Extract the mean NDVI value *strictly within the mangrove polygons* we generated above
      ndvi.data = lapply(list(NDVI.w, NDVI.ne, NDVI.se), function(x) {
        exact_extract(x, me.table, stack_apply = T, fun = 'mean',
                      append_cols = c('buffer', 'lon', 'lat'), progress = F) %>%
          drop_na() %>%
          # Reshape from wide (columns for each date) to long (date as a row variable)
          pivot_longer(cols = starts_with('mean'), values_to = 'mean', names_to = 'time') %>%
          # String manipulation: extract the Julian date from the raster filename text
          mutate(time = str_sub(time, start = 40, end = 46)) %>%
          mutate(month = as.Date(time, '%Y%j'))
      }) %>%
        bind_rows()
      
      # Export the final NDVI dataset for this specific year
      write.table(ndvi.data, file = paste0("./data/ndvi.data.", y, ".csv"),
                  row.names = F, sep = ",")
    }
  }
  
  # Timing and Progress Logging
  end.time = Sys.time()
  time.taken = round(end.time - start.time, 2)
  print(paste("Year", y, ": Done!", sep = " "))
  print(time.taken)
})