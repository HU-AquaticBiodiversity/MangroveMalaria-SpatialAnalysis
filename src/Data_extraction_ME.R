### NOTE: This could should only be run through high-perfomance computing ###
### as the dataset might be too large ###

##load necessary libraries
library(dplyr)
library(sf)
library(sp)
library(rnaturalearth)
library(ggplot2)
library(lwgeom)
library(tidyverse)
library(geosphere)
library(osmdata)
library(units)
library(parallel)
library(doParallel)
library(data.table)
library(exactextractr)
library(terra)
library(pbapply)
library(ape)
library(tidyterra)
library(ggsflabel)

years = c(1996, 2007:2010, 2015:2020)

##--------------------------------------------
## Download mangrove raster data and unzip
##   - Only use if files aren't already there
##--------------------------------------------
#outDir = "./res/mangrove_raster/"
#dir.create(outDir)
#mangrove.files = sapply(years, function(x) paste0("gmw_v3_", x, '_gtiff.zip'))
#download_zenodo(doi = "10.5281/zenodo.6894272", path = outDir,
#                mangrove.files)
## unzip files. NOTE: might be fast to just do this manually
#pblapply(list.files(path = outDir, pattern = "*.zip"),
#       function(i) {unzip(paste0(outDir,i), exdir=outDir, overwrite = F)})

##-----------------------------------------------
## Loading Data
##-----------------------------------------------
## load polygons files and extract those for maximum radius (50 km)
coastal.PR.polygons = st_read("./data/polygons_pr_malaria.shp")
coastal.PR.polygons.50 = coastal.PR.polygons %>%
  filter(buffer == 50)
## load prevalence data
coastal.PR = read.csv("./data/coastal.PR.final.csv", sep = ",", header = T)

##----------------------------------------------
## Extracting mangrove surfaces for each year (y)
##----------------------------------------------
lapply(years, function(y) {
  start.time = Sys.time()

  # load names of mangrove tiles (part of overall raster)
  tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", y),
                     pattern = ".tif", full.names = T)
  # load information of mangrove raster in year y
  tile.list = lapply(tiles, rast)

  # select coordinates matching the year y or y + 1
  coord = coastal.PR %>% filter(year_start %in% c(y, y+1)) %>%
    select(lon, lat) %>% unite("dummy", lon, lat, sep = "_", remove = F)
  # subset polygons by coordinates to only calculate necessary values
  polygons = coastal.PR.polygons.50 %>% 
    unite("dummy", lon, lat, sep = "_", remove = F) %>%
    filter(dummy %in% coord$dummy) %>%
    select(-dummy)
  nrow(polygons)

  ## CONDITION: if there are location for the selected year (y) to continue. 
  ##            Otherwise, don't do anything
  if(nrow(polygons) > 0) {

    ## for each location calculate mangrove extent
    me.table = pblapply(1:nrow(polygons), function(p) {
      # select location and its 50-km polygon
      poly = polygons[p,]

      # determine which mangrove tiles overlap with the 50-km polygon
      tile.overlaps = lapply(1:length(tile.list), function(x) {
        nrow(st_crop(poly, tile.list[[x]])) != 0
      }) %>%  unlist()

      # select all polygons (1-50 km) in the same location as the 50-km polygon 
      poly.1.50 = coastal.PR.polygons %>% filter(lon == poly$lon & lat == poly$lat)

      ## CONDITION: Are any of the tiles overlapping with the 50-km polygon?
      ##  if yes: produce a vector of 0
      ##  if no: calculate mangrove area inside 1-50 km polygons
      if(sum(tile.overlaps)==0) {
        mc = rep(0, 50)
      } else {
        # select tile(s) within 50-km polygon
        tile.sel = mosaic(sprc(tile.list[tile.overlaps]))
        
        ## 1. surface area for each cell of the raster tile is calculated through
        ##    'cellSize' with NA's being masked
        ## 2. the surface areas (i.e. mangrove cover) of all cells within 
        ##    the 1-50 km polygons are calculated through 'exact_extract'
        mc = exact_extract(cellSize(tile.sel, unit = 'km', mask = T),
                           poly.1.50, fun = 'sum', progress = F)
      }
      
      # print number of p vector to see progress
      print(paste0(p, "/", nrow(polygons)))
      
      # add calculated values to polygon data
      poly.1.50 %>% mutate(mangrove.cover = mc)
    }) %>% bind_rows() %>%
      st_drop_geometry()
    me.table
    ## EXPORT table for each year
    write.table(me.table, paste0('./data/MangroveCover/MangroveCover_', y, '.csv'), sep = ",", row.names = F)

  }
  
  end.time = Sys.time()
  time.taken = round(end.time - start.time,2)
  print(time.taken)
  print(paste("Year", y, ": Done!", sep = " "))
})
