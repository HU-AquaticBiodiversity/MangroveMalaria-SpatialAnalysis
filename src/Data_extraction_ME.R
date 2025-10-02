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


##-----------------------------------------------
## Plot example figure: calculation of mangrove variables
##-----------------------------------------------
## load polygons files and extract those for maximum radius (50 km)
polygon.ex = coastal.PR.polygons %>%
  # select specific site, here: Saloum Delta in Senegal
  filter(lat > 14 & lat < 14.1) %>%
  mutate(hJust = -(buffer-5)/9)

# load names of mangrove tiles (part of overall raster)
tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", 2008),
                   pattern = ".tif", full.names = T)
# load information of mangrove raster in year y
tile.list = lapply(tiles, rast)

# determine which mangrove tiles overlap with the 50-km polygon
tile.overlaps = pblapply(1:length(tile.list), function(x) {
  nrow(st_crop(polygon.ex[50,], tile.list[[x]])) != 0
}) %>%  unlist()

# select tile(s) within 50-km polygon
tile.sel = mosaic(sprc(tile.list[tile.overlaps]))

# import background map
ne.land = st_read("./res/ne_10m_land/ne_10m_land.shp")

# plot mangrove extent
p.mangrove = ggplot() +
  #back groun map
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  #plot mangrove land cover
  geom_spatraster(data = tile.sel, aes(fill = GMW_N14W017_2008_v3)) +
  # make mangroves green
  scale_fill_gradient(low = "grey", high = "green4", na.value = "transparent") +
  # plot radii, but only select a few
  geom_sf(data = polygon.ex %>%
            filter(buffer %in% c(1,10,20,30,40,50)), alpha = 0) +
  # add labels to radii
  geom_sf_label(data = polygon.ex %>%
                 filter(buffer %in% c(1,10,20,30,40,50)),
                # spread out labels, so they do not overlap
               aes(label = paste(buffer, "km"), hjust = hJust, vjust = -1),
               size = 2.5) +
  # cut out map to specific area
  coord_sf(xlim = c(-17.2, -16.1), ylim = c(13.5, 14.6), expand = FALSE) +
  # change axis labels
  xlab("longitude")+ ylab("latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill='skyblue', colour='red'),
        legend.position = 'none')

# find 2008 NDVI data
ndvi.files.w = list.files(
  "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/West",
  full.names = T, 
  pattern = paste0("*NDVI_doy", 2008,".*tif"))
# load information of mangrove raster in 2008
NDVI.w = lapply(ndvi.files.w, rast)

# crop mangrove area inside the 50-km range
mpoly.ex.50 = as.polygons(crop(tile.sel, polygon.ex[50,], mask = T))
# crop NDVI data inside the 50-km mangrove area
NDVI.crop = crop(NDVI.w[[15]], mpoly.ex.50, mask = T)

# plot NDVI
p.ndvi = ggplot() +
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  geom_spatraster(data = NDVI.crop, 
                  aes(fill = MOD13Q1.061__250m_16_days_NDVI_doy2008225_aid0001)) +
  scale_fill_gradient(low = "yellow", high = "green", na.value = "transparent") +
  geom_sf(data = polygon.ex %>%
            filter(buffer %in% c(1,10,20,30,40,50)), alpha = 0) +
  geom_sf_label(data = polygon.ex %>%
                  filter(buffer %in% c(1,10,20,30,40,50)),
                aes(label = paste(buffer, "km"), hjust = hJust, vjust = -1),
                size = 2.5) +
  coord_sf(xlim = c(-17.2, -16.1), ylim = c(13.5, 14.6), expand = FALSE) +
  xlab("longitude")+ ylab("latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill='skyblue', colour='red')) +
  labs(fill = "NDVI")

# put two plots in single panel
p.joint = plot_grid(p.mangrove, p.ndvi, labels = c('A', 'B'), 
                    label_size = 20, rel_widths =c(46.3,53.7))

# export joint plot
ps.options(family = "Arial")
ggsave(filename="mangroveVars_v1.svg", device = "svg", width = 270, height = 135, units = "mm")
ggsave(filename="mangroveVars_v1.png", device = "png", width = 270, height = 135, units = "mm")
ggsave(filename="mangroveVars_v1.pdf", device = cairo_pdf, width = 270, height = 135, units = "mm")
cairo_ps(filename="mangroveVars_v1.eps", width = 10.63, height = 5.315, #in inches
         pointsize = 10, fallback_resolution = 2400)
p.joint
dev.off()