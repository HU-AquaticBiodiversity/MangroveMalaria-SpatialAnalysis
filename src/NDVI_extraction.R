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
library(viridis)


years = c(2007:2010, 2015:2020)
numberOfCores = detectCores()

##-----------------------------------------------
## Loading Data
##-----------------------------------------------
## load polygons files and extract those for maximum radius (50 km)
coastal.PR.polygons = st_read("./data/polygons_pr_malaria_dhs.shp") %>%
  bind_rows(st_read("./data/polygons_pr_malaria.shp"))
coastal.PR.polygons.50 = coastal.PR.polygons %>%
  filter(buffer == 50)
## load prevalence data
coastal.PR = read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = T) %>%
  bind_rows(read.csv("./data/coastal.PR.final.csv", sep = ",", header = T) %>%
              mutate(site_id = as.factor(site_id)))

##----------------------------------------------
## Extracting mangrove surfaces for each year (y)
##----------------------------------------------
lapply(years, function(y) {
  start.time = Sys.time()

  # select coordinates matching the year y
  coord = coastal.PR %>% filter(year_start == y) %>%
    select(lon, lat) %>% unite("dummy", lon, lat, sep = "_", remove = F)
  
  # subset polygons by coordinates to only calculate necessary values
  polygons = coastal.PR.polygons.50 %>% 
    unite("dummy", lon, lat, sep = "_", remove = F) %>%
    filter(dummy %in% coord$dummy) %>%
    select(-dummy)

  ## CONDITION: if there are location for the selected year (y) to continue. 
  ##            Otherwise, don't do anything
  if(nrow(polygons) > 0) {

    # load names of mangrove tiles (part of overall raster)
    tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", y),
                       pattern = ".tif", full.names = T)
    # load information of mangrove raster in year y
    tile.list = lapply(tiles, rast)
    
    ## for each location calculate mangrove extent
    me.table = mclapply(1:nrow(polygons), function(p) {

      # print number of p vector to see progress
      print(paste0(p, "/", nrow(polygons)))
      
      # select location and its 50-km polygon
      poly = polygons[p,]

      # determine which mangrove tiles overlap with the 50-km polygon
      tile.overlaps = lapply(1:length(tile.list), function(x) {
        nrow(st_crop(poly, tile.list[[x]])) != 0
      }) %>%  unlist()
      
      ## CONDITION: Are any of the tiles overlapping with the 50-km polygon?
      ##  if yes: check mangrove overlap
      ##  if no: create empty data frame
      if(sum(tile.overlaps)!=0) {
        
        # select tile(s) within 50-km polygon
        tile.sel = mosaic(sprc(tile.list[tile.overlaps]))
      
        # select all polygons (1-50 km) in the same location as the 50-km polygon 
        poly.1.50 = coastal.PR.polygons %>% filter(lon == poly$lon & lat == poly$lat) %>%
          mutate(overlap = exact_extract(tile.sel, ., fun = 'sum', progress = F) > 0) %>%
          filter(overlap)
        
        ## CONDITION: Are any of the polygons actually overlapping with the mangroves?
        ##  if yes: render overlap polygons (but only for those that overlap)
        ##  if no: create empty data frame
        if(nrow(poly.1.50) > 0) {
          
          tiles.crop = lapply(1:nrow(poly.1.50), function(x) {

            as.polygons(crop(tile.sel, poly.1.50[x,], mask = T)) %>%
              st_as_sf() %>%
              # add metadata
              mutate(buffer = poly.1.50[x,]$buffer,
                     lat = poly.1.50[x,]$lat,
                     lon = poly.1.50[x,]$lon)
          }) %>% .[!is.na(.)] %>% bind_rows() %>%
            # remove useless column
            select(-starts_with("GMW"))
          tiles.crop
        
        } else {NA}
      } else {NA}
    }, mc.cores = numberOfCores) %>% .[!is.na(.)] %>%
    bind_rows()
    
    if(nrow(me.table) > 0) {
      ndvi.files.w = list.files(
        "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/West",
        full.names = T, 
        pattern = paste0("*NDVI_doy", y,".*tif"))
      ndvi.files.ne = list.files(
        "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/NorthEast",
        full.names = T, 
        pattern = paste0("*NDVI_doy", y,".*tif"))
      ndvi.files.se = list.files(
        "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/SouthEast",
        full.names = T, 
        pattern = paste0("*NDVI_doy", y,".*tif"))
      
      NDVI.w = rast(ndvi.files.w)
      NDVI.ne = rast(ndvi.files.ne)
      NDVI.se = rast(ndvi.files.se)
      
      ndvi.data = lapply(list(NDVI.w, NDVI.ne, NDVI.se), function(x) {
        exact_extract(x, me.table, stack_apply = T, fun = 'mean',
                      append_cols = c('buffer', 'lon', 'lat'), progress = F) %>%
          drop_na() %>%
          pivot_longer(cols = starts_with('mean'), values_to = 'mean', names_to = 'time') %>%
          mutate(time = str_sub(time, start = 40, end = 46)) %>%
          mutate(month = as.Date(time, '%Y%j'))
      }) %>%
        bind_rows()
      ndvi.data
      
      write.table(ndvi.data, file = paste0("./data/ndvi.data.", y, ".csv"),
                  row.names = F, sep = ",")
    }

    
  }
  
  end.time = Sys.time()
  time.taken = round(end.time - start.time,2)
  print(paste("Year", y, ": Done!", sep = " "))
  print(time.taken)
})