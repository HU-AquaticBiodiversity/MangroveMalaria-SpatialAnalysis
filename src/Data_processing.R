library(malariaAtlas)
library(dplyr)
library(sf)
library(sp)
library(rnaturalearth)
library(tidyverse)
library(osmdata)
library(zen4R)
library(rdhs)
library(pbapply)
library(terra)
library(tidyterra)
library(viridis)

##----------------------------------------------------------
## Creat low-res shapefile to locate mangrove scenes in USGS
##----------------------------------------------------------
# load names of mangrove tiles (part of overall raster)
tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_mng_union/"),
                   pattern = ".tif", full.names = T)
# load information of mangrove raster in year y
tile.list = pblapply(tiles, rast)

# aggregate raster to simplify shapefile
tiles.aggr = pblapply(tile.list, function(x) aggregate(x, fact = 4500, fun = 'mean', na.rm = T))

# split rasters into separate regions to meet size requirement for NASA download
tiles.mosaic.west = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(-20, 20, -37, 30))
tiles.mosaic.northeast = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(20, 55, 5, 30))
tiles.mosaic.southeast = crop(mosaic(sprc(tiles.aggr), fun = 'first'), ext(20, 55, -37, 5))
mangrove.shp.west = st_as_sf(as.polygons(tiles.mosaic.west))
mangrove.shp.northeast = st_as_sf(as.polygons(tiles.mosaic.northeast))
mangrove.shp.southeast = st_as_sf(as.polygons(tiles.mosaic.southeast))
write_sf(mangrove.shp.west, dsn = './data/mangrove_vcts_aggr/mangrove_vec_west.shp', append = F)
write_sf(mangrove.shp.northeast, dsn = './data/mangrove_vcts_aggr/mangrove_vec.northeast.shp', append = F)
write_sf(mangrove.shp.southeast, dsn = './data/mangrove_vcts_aggr/mangrove_vec.southeast.shp', append = F)
# code to be used to plot rasters for visualisation #
ggplot() +  
  geom_sf(data = mangrove.shp.southeast) +
  scale_fill_viridis()



##--------------------------------------##
##--------------------------------------##
## DOWNLOAD MALARIA DATA and save
##--------------------------------------##
##--------------------------------------##
 
# overview of all relevant countries
country.list = read.csv("./data/country_table.csv")
# make expanded list of years to include data from the year after (year + 1)
years.expanded = c(1996:1997, 2007:2011, 2015:2021)

# download malaria prevalences from Malaria Atlas Project
malaria.PR.raw = getPR(continent = "Africa", species = "Pf") %>%
  # filter by countries relevant to this study #
  filter(country_id %in% country.list$ISO3,
         year_start %in% years.expanded | year_end %in% years.expanded)

## add coordinates from DHS website: 
#   - data access needs to be requested through DHS website
malaria.PR.raw.dhs = fillDHSCoordinates(malaria.PR.raw %>% filter(is.na(latitude)),
                                        email = "cruzmamo@googlemail.com",
                                  project = "Malaria prevalence in mangrove forests",
                                  timeout = 60)

malaria.PR = malaria.PR.raw %>% 
  filter(!is.na(latitude)) %>%
  filter(!is.na(pr)) %>%
  as.data.frame()

malaria.PR.dhs = malaria.PR.raw.dhs %>% 
  filter(!is.na(latitude)) %>%
  filter(!is.na(pr)) %>%
  as.data.frame() %>%
  mutate(site_id = dhs_id)

nrow(malaria.PR)
nrow(malaria.PR.dhs)
# ==> almost triple the amount of PR values

#----------------------------------------------------------------------------------
## Only keep PR data in 50 km radius for different GIS layers ##
#----------------------------------------------------------------------------------
# define dataset to be used #
#  this code uses the publicly available malariaAtlas data (excluding
#  the DHS dataset). Use malaria.PR.dhs instead of malaria.PR to include
#  DHS data
PR = malaria.PR.dhs

# call coastline shape data #
coastline = ne_coastline(scale = 10, returnclass = "sf")
# limit countries to those with available prevalence data
countries.available = country.list %>% filter(ISO3 %in% PR$country_id)

# filter dataset to only include point up to 50 km off the coast
# + render shape files to infer raster data in 50-km radius down the line
coastal.PR = pblapply(1:nrow(countries.available), function(c) {
  disease.data = PR %>%
    filter(country_id == countries.available$ISO3[c]) %>%
    as.data.frame() %>%
    ## transform disease data to correct CRS ## 
    st_as_sf(coords = c('longitude', 'latitude'), 
             crs = st_crs(4326), dim = "XY") %>%
    st_transform(countries.available$CRS[c]) 

  ## Filter location in 50-km range of coastline ##
  
  ## apply CRS to coastline shape file ##
  coastline.crs = coastline %>%
    st_transform(countries.available$CRS[c])
  ## filter multilinestrings that are invalid, e.g. only a single point ##
  coastline.filtered = coastline.crs[ifelse(is.na(st_is_valid(coastline.crs)), F, T),]
  coastline.alldist = st_distance(coastline.filtered, disease.data)
  
  ## coastline data is multilinestring. Minimal distance is extracted for each point ##
  coastline.dist =  sapply(1:nrow(disease.data), function(k) {          
    ## calculate closest distance to any of the line strings ##
    min(coastline.alldist[,k], na.rm = T)/1000
  })
  
  disease.filtered = disease.data %>%
    mutate(coastline.dist = coastline.dist) %>%
    filter(coastline.dist <= 50) 
})

# check for which countries coastal PR is available
data.available = sapply(1:nrow(countries.available), function (c) {
  nrow(coastal.PR[[c]]) != 0
})

# merge all country datasets and export to a csv file
coastal.PR.final = bind_rows(
  lapply(coastal.PR, function(x) {
    x %>% 
      st_transform(crs = 4326) %>%
      dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                    lat = sf::st_coordinates(.)[,2]) %>%
      st_drop_geometry()
      })
  )

write.table(coastal.PR.final, 
            # name coastal.PR.dhs.final.csv for DHS
            file = "./data/coastal.PR.dhs.final.csv",
            sep = ",", row.names = F)

# create 1-50 km buffer around all locations to extract raster values
# in surrounding areas for downstream analyses
coastal.PR.polygons = bind_rows(pblapply(1:50, function(y) {
  bind_rows(
    lapply(coastal.PR[data.available], function(x) {
      # extract lon/lat data
      points = x %>%
        dplyr::select(geometry) %>%
        st_transform(crs = 4326) %>%
        dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                      lat = sf::st_coordinates(.)[,2]) %>%
        st_drop_geometry()
      # buffer point by 50 kms, transform back to WGS84, add lon/lat, and pick
      #   unique values
      x %>%
        st_buffer(y * 1000) %>%
        dplyr::select(geometry) %>%
        st_transform(crs = 4326) %>%
        bind_cols(points) %>%
        unique() %>%
        mutate(buffer = y)
    })
  )
  })
  )

# save polygons
st_write(coastal.PR.polygons,
         # name polygons_pr_malaria_dhs.shp for DHS
         "./data/polygons_pr_malaria_dhs.shp",
         append = F)
nrow(coastal.PR.polygons %>% filter(buffer == 50))

# plot polygons for visualisation
p = ggplot(coastal.PR.polygons %>% filter(buffer == 50) ) +
  geom_sf() +
  #geom_sf(data =coastline) +
  theme_minimal()
p
