##load necessary libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(sf)
library(sp)
library(lwgeom)
library(tidyverse)
library(geosphere)
library(osmdata)
library(units)
library(pbapply)
library(data.table)
library(terra)
library(tidyterra)
library(vroom)
library(exactextractr)
library(caret)
library(rnaturalearth)
library(purrr)
library(naniar)
library(parallel)
library(cowplot)
library(doParallel)

# load prevalence data
country.list = read.csv("./data/country_table.csv")
coastal.PR = read.csv("./data/coastal.PR.final.csv", sep = ",", header = T) %>%
  mutate(site_id = as.factor(site_id)) %>%
  bind_rows(read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = T))

#load coastline data
coastline = ne_coastline(scale = 10, returnclass = "sf")

# limit countries to those with available prevalence data
countries.available = country.list %>% filter(ISO3 %in% coastal.PR$country_id)

##----------------------------------------------------------------------------------
##----------------------------------------------------------------------------------
## Data extraction on local machine
##----------------------------------------------------------------------------------
##----------------------------------------------------------------------------------
## load malaria polygons with 10-km radius ##
all_polygons.10 = st_read("./data/polygons_pr_malaria.shp") %>%
  bind_rows(st_read("./data/polygons_pr_malaria_dhs.shp")) %>%
  filter(buffer %in% c(5,10,20))

nrow(all_polygons.10)
#--> 4260 different locations

##------------------------------------------------------------##
##  extract values for 'constant' rasters,
##  -- i.e. rasters, which are only available for a single year
##------------------------------------------------------------##

 # load rasters for saltwater species and merge #
prob.list = sprc(rast('./res/mosquito_distributions/2017_Anopheles_melas.Mean_Decompressed.geotiff'),
                 rast('./res/mosquito_distributions/2017_Anopheles_merus.Mean_Decompressed.geotiff'))
prob.salt = terra::mosaic(prob.list, fun = "sum")

 # load rasters for 'constant' health management data #
t_city = rast('./res/Health_management/201501_Global_Travel_Time_to_Cities_2015/201501_Global_Travel_Time_to_Cities_2015.tif')
t_healthcare_motor = rast('./res/Health_management/202001_Global_Motorized_Travel_Time_to_Healthcare_2019/202001_Global_Motorized_Travel_Time_to_Healthcare_2019.tif')
t_healthcare_walk = rast('./res/Health_management/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.tif')
 
 # code to be used to plot rasters for visualisation #
ggplot() +  
  geom_spatraster(data = t_city) +
  scale_fill_viridis()

 # extracting median and sd values for all rasters #
const.raster.data = pblapply(
  list(prob.salt, t_city, t_healthcare_motor, t_healthcare_walk),
  FUN = function(y) {
    as.data.frame(
      exact_extract(y, all_polygons.10,
              fun= c('median', 'stdev'))
      )
    }) %>% bind_cols() %>%
  mutate(lon = all_polygons.10$lon,
         lat = all_polygons.10$lat,
         buffer = all_polygons.10$buffer) %>% 
  filter(buffer == 10) %>% dplyr::select(-buffer)

 # change column names into something more meaningful #
names(const.raster.data) = c('mosq.median', 'mosq.sd', 't_city.median', 't_city.sd', 
                               't_healthcare_motor.median', 't_healthcare_motor.sd',
                               't_healthcare_walk.median', 't_healthcare_walk.sd',
                               'lon', 'lat')


gg_miss_upset(const.raster.data)
write.table(const.raster.data, "./data/const.raster.data.csv", sep = ",", row.names = F)

##-----------------------------------------------------------##
##  extract values for other rasters,
##  -- i.e. rasters, which are available for a multiple years
##-----------------------------------------------------------##
 # load all raster files for health management #
ITN_use = rast(list.files('./res/health_management/2024_GBD2023_Africa_ITN_2000/',
                         pattern = "2024_GBD2023", full.names = T))
ITN_access = rast(list.files('./res/health_management/ITN_2000_access_mean/',
                                 pattern = "ITN", full.names = T))
IRS = rast(list.files('./res/health_management/2024_GBD2023_Africa_IRS_2000/',
                      pattern = "2024_GBD2023", full.names = T))
AMT = rast(list.files('./res/health_management/2024_GBD2023_Global_Antimalarial_EFT_2000/',
                      pattern = "2024_GBD2023", full.names = T))
POP = rast(list.files('./res/population_density/', full.names = T))

# code to be used to plot rasters for visualisation #
ggplot() +  
  geom_spatraster(data = POP[[1]]) +
  scale_fill_viridis()

 # extracting median and sd values for all rasters #
annual.raster.data = lapply(list(ITN_use, ITN_access, IRS, AMT, POP), function(x) {
  exact_extract(x, all_polygons.10,
                fun= c('median', 'stdev'))
  }) %>% bind_cols() %>%
  mutate(lon = all_polygons.10$lon,
         lat = all_polygons.10$lat,
         buffer = all_polygons.10$buffer) %>%
  # rearrange table to extract metadata in names of files #
  pivot_longer(!lon & !lat & !buffer, names_to = 'dummy', values_to = 'value') %>%
  mutate(
    # extract year of measurement #
    year_start = as.numeric(ifelse(str_ends(dummy,"\\d"), 
                        str_sub(dummy, -4),
                        str_extract(dummy, "(\\d)+"))),
    # extract statistics: mean or sd #
    statistic = str_extract(dummy, "[^.]+"),
    # extract variable name
    name = str_remove_all(str_split(dummy, fixed("."), simplify = T)[,2],
                      "(\\d)+"),
    # alter names of population density layers
    name = ifelse(name == "Af_more_qr_clust_" | name == "Ex_AF_Cot_NoNeg_",
                  "pop_dens", name)
    ) %>%
  mutate(name = str_remove(str_remove(str_replace(name, '__','_'), "^\\_"), "\\_$")) %>%
  select(-dummy) %>%
  pivot_wider(names_from = c(name, statistic), values_from = value)

write.table(annual.raster.data, "./data/annual.raster.data.csv", sep = ",", row.names = F)

##-----------------------------------------------------------##
##  extract values for monthly rasters,
##  -- i.e. rasters, which are available for multiple years AND months
## 1. Cropland
##-----------------------------------------------------------##

# create vector for names of layers (year)
layer.names = paste0(
  "lc_",
  str_extract(list.files('./res/landcover/East/', full.names = T), pattern = "\\d{4}")
  )

## import raster layers: these layers were downloaded separately 
#  for West and East Africa
lc_east = rast(list.files('./res/landcover/East/', full.names = T), 
               subds = 'lccs_class')
lc_west = rast(list.files('./res/landcover/West/', full.names = T), 
               subds = 'lccs_class')
# re-classify rasters: cropland = 1, other = NA
lc_east_bin = ifel(lc_east %in% c(10,20,30,40), 1, NA)
lc_west_bin = ifel(lc_west %in% c(10,20,30,40), 1, NA)
#split polygons into West and East to match rasters
poly_west = all_polygons.10[all_polygons.10$lon < 28,]
poly_east = all_polygons.10[all_polygons.10$lon > 28,]

# extract surface area of croplands in km²
lc_west = lapply(1:nlyr(lc_west_bin), function(x) {
  exact_extract(cellSize(lc_west_bin[[x]], unit = 'km', mask = T),
                   poly_west, fun = 'sum', progress = T)
})

lc_east = lapply(1:nlyr(lc_east_bin), function(x) {
  exact_extract(cellSize(lc_east_bin[[x]], unit = 'km', mask = T),
                poly_east, fun = 'sum', progress = T)
})

#create dataframe with cropland areas
lc_west.tab = do.call(cbind, lc_west)
lc_east.tab = do.call(cbind, lc_east)
colnames(lc_west.tab) = layer.names
colnames(lc_east.tab) = layer.names
lc_all = st_drop_geometry(rbind(cbind(poly_west, lc_west.tab),
                       cbind(poly_east, lc_east.tab)))

#export cropland dataframe
write.csv(lc_all, "./data/cropland.csv", row.names = F)

##-----------------------------------------------------------##
## 2. Weather data
##-----------------------------------------------------------##

lyr.names = str_sub(list.files("./res/climate/"), end = -10)

weather.raster = terra::project(
  # change projection to WGS84
  rast(list.files("./res/climate/", full.names = T)),
  '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
  # to speed up calculation
  use_gdal = T, threads = 7
  )

names(weather.raster) = lyr.names

weather.data = exact_extract(weather.raster, all_polygons.10, fun= 'mean') %>% 
  bind_cols(all_polygons.10 %>% st_drop_geometry() %>% dplyr::select(lon, lat, buffer))

#export cropland dataframe
write.csv(weather.data, "./data/weather_data.csv", row.names = F)

##-----------------##
## HDI ##
##-----------------##
 # load shapefiles #
## shp.subnat = st_read("./res/GDL_Shapefiles_V6.3/GDL Shapefiles V6.3 large.shp") %>%
##  dplyr::select(-iso_code)

nHDI.table = read.csv("./res/GDL_Shapefiles_V6.3/GDL-Subnational-HDI-data.csv") %>%
  pivot_longer(X1990:X2021, names_to = 'year', values_to = 'nHDI') %>%
  mutate(year = as.numeric(str_sub(year, start = 2))) %>%
  dplyr::select(GDLCODE, year, nHDI)

# export nHDI data to file
write.table(nHDI.table, "./data/nHDI.data.csv", sep = ",", row.names = F)

##----------------------------------------- ##
##----------------------------------------- ##
## merging all data together ##
##----------------------------------------- ##
##----------------------------------------- ##

## load shapefiles for GDL data ##
## NOTE: files are partially erronous --> use 'st_make_valid'
shp.subnat = st_read("./res/GDL_Shapefiles_V6.3/GDL Shapefiles V6.3 large.shp") %>%
  dplyr::select(-iso_code) %>%
  st_make_valid()

## adapt cropland dataframe
lc_all = read.csv("./data/cropland.csv") %>%
  pivot_longer(cols = starts_with("lc_"),
               names_to = "year",
               values_to = "lc_crop",
               values_drop_na = TRUE) %>%
  mutate(year = as.numeric(str_extract(year, pattern = "\\d{4}"))) %>%
  # arrange 5 and 20-km values as additional columns
  mutate(buffer = ifelse(buffer == 10, "lc_crop", paste0("lc_crop", "_", buffer))) %>%
  pivot_wider(names_from = "buffer", values_from = "lc_crop")

## load mangrove cover data frame ##
mangrove.cover = vroom(list.files("data/MangroveCover", full.names = T),
      # use file name as ID
      id = 'year') %>%
  # extract year from file name: "\\d{4}" --> four digits
  mutate(year = as.numeric(str_extract(year, pattern = "\\d{4}")))

##  add constant raster data to prevalence data table ##
coastal.PR.all = read.csv("./data/coastal.PR.final.csv", sep = ",", header = T) %>%
  mutate(site_id = as.character(site_id)) %>%
  bind_rows(read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = T)) %>%
  dplyr::select(pr, lat, lon, examined, method, coastline.dist,
                month_start, year_start, month_end, year_end) %>%
  st_as_sf(coords = c('lon', 'lat'), remove = F,
           crs = st_crs(4326), dim = "XY") %>%
  st_join(shp.subnat) %>%
  st_drop_geometry() %>%
  
  ##  add constant raster data ##
  left_join(read.csv("./data/const.raster.data.csv")) %>%
  # add cropland data
  left_join(lc_all, by = c("year_start" = "year", 'lon', 'lat')) %>%
  # add annual raster data
  left_join(
    read.csv("./data/annual.raster.data.csv") %>%
      # select only default buffer of 10 km
      filter(buffer == 10) %>%
      dplyr::select(-buffer) %>%
      # make 5 and 20-km buffers additional columns for pop_dens
      mutate(
        pop_dens_median_5 = read.csv("./data/annual.raster.data.csv") %>%
          filter(buffer == 5) %>% dplyr::select(pop_dens_median) %>% unlist(),
        pop_dens_median_20 = read.csv("./data/annual.raster.data.csv") %>%
          filter(buffer == 20) %>% dplyr::select(pop_dens_median) %>% unlist()
        ),
            by = c('year_start', 'lon', 'lat')) %>%
  # add mangrove cover
  right_join(mangrove.cover, by = c("year_start" = "year", 'lon', 'lat'),
             relationship = "many-to-many") %>%
  filter(!is.na(pr)) %>%
  left_join(mangrove.cover %>% 
               rename(mangrove.cover.min1 = mangrove.cover) %>%
               mutate(year = year + 1) ,
             by = c("year_start" = "year", 'lon', 'lat', 'buffer'),
             relationship = "many-to-many") %>%
  # add nHDI
  left_join(read.csv("./data/nHDI.data.csv"),
            by = c("gdlcode" = "GDLCODE", "year_start" = "year")) %>%
  select(- continent)

## -------------------------------------------------- ##
## calculating mean ndvi for each row of disease data ##
## -------------------------------------------------- ##
# load ndvi data
ndvi.data = list.files("./data/NDVI.data/", full.names = T) %>%
  map_df(~read_csv(.))

# detect the number of cores
n.cores <- detectCores()
n.cores

# set up parallel computation for better speed (parLapply function)
clust <- makeCluster(n.cores - 1)
clusterExport(clust, c("coastal.PR.all", "ndvi.data"))
system.time({
  ndvi.mean = parSapply(clust, 1:nrow(coastal.PR.all), function(x) {

    # need to load R packages again in cluster
    library(dplyr)
    library(lubridate)
    a = ndvi.data %>%
      # only select relevant rows of NDVI data
      filter(buffer == coastal.PR.all[x,]$buffer &
               lat == coastal.PR.all[x,]$lat &
               lon == coastal.PR.all[x,]$lon) %>%
      # select rows between dates specified in prevalence data
      filter(between(as.Date(month), 
                     # start month
                     as.Date(paste(coastal.PR.all[x,]$year_start,
                                   coastal.PR.all[x,]$month_start, "1", sep = "-")),
                     # end month +1 to include all days of end month
                     as.Date(paste(coastal.PR.all[x,]$year_end,
                                   coastal.PR.all[x,]$month_end,
                                   "1", sep = "-")) %m+% months(1)))
    #calculate average NDVI for selected NDVI data
    a
    mean(a$mean)
  })
})
# cluster needs to be stopped after run 
stopCluster(clust)

## ----------------------------------------------------- ##
## calculating weather data for each row of disease data ##
## ----------------------------------------------------- ##
weather = function(b) {
  ## --> weather table needs to be prepared
  weather.table = read.csv("./data/weather_data.csv") %>%
    filter(buffer == b) %>%
    dplyr::select(-buffer) %>%
    # pivot to process column names
    pivot_longer(c(!lat & ! lon), names_to = "measure", values_to = "value") %>%
    # remove unnecessary text
    mutate(measure = str_replace(
      str_replace(
        str_replace(measure, "_1991\\.2020", ""),
        "mean\\.1month_", ""),
      "_Global_ea", "")) %>%
    # split information
    separate(measure, c('measure.1', 'measure.2', 'date'), sep = "_") %>%
    # merge variable info
    unite(measure, c(measure.1, measure.2)) %>%
    # place variables next to each other
    pivot_wider(names_from = measure, values_from = value) %>%
    # split month and year
    mutate(date = paste(str_sub(date, 1, 4), str_sub(date, 5, 6), "01", sep = "-"))
  
  ## 1. Run for weather WITHIN survey period ##
  # make cluster for parallel computation
  clust <- makeCluster(n.cores - 1)
  clusterExport(clust, c("coastal.PR.all", "weather.table"), envir=environment())
  
  system.time({
    weather.data = parLapply(clust, 1:nrow(coastal.PR.all), function(x) {
      # need to load R packages again in cluster
      library(dplyr)
      library(lubridate)
      
      a = weather.table %>%
        filter(lat == coastal.PR.all[x,]$lat &
                 lon == coastal.PR.all[x,]$lon) %>%
        filter(between(as.Date(date), 
                       as.Date(paste(coastal.PR.all[x,]$year_start,
                                     coastal.PR.all[x,]$month_start, "1", sep = "-")),
                       as.Date(paste(coastal.PR.all[x,]$year_end,
                                     coastal.PR.all[x,]$month_end,
                                     "1", sep = "-")) %m+% months(1))) %>%
        select(anomaly_2t:mean_tp)
      colMeans(a)
    }) %>%
      bind_rows()
  })
  stopCluster(clust)
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
  ## 1. Run for weather BEFORE survey period ##
  # make cluster for parallel computation
  clust <- makeCluster(n.cores - 1)
  clusterExport(clust, c("coastal.PR.all", "weather.table"), 
                envir=environment())
  
  system.time({
    weather.data.min6 = parLapply(clust, 1:nrow(coastal.PR.all), function(x) {
      # need to load R packages again in cluster
      library(dplyr)
      library(lubridate)
      
      a = weather.table %>%
        filter(lat == coastal.PR.all[x,]$lat &
                 lon == coastal.PR.all[x,]$lon) %>%
        filter(between(as.Date(date), 
                       as.Date(paste(coastal.PR.all[x,]$year_start,
                                     coastal.PR.all[x,]$month_start, 
                                     "1", sep = "-")) %m-% months(6),
                       as.Date(paste(coastal.PR.all[x,]$year_end,
                                     coastal.PR.all[x,]$month_end,
                                     "1", sep = "-")) )) %>%
        select(anomaly_2t:mean_tp)
      colSums(a)
    }) %>%
      bind_rows() %>%
      rename_with(~ paste(., "6m", sep = "_"))
  })
  
  stopCluster(clust)
  
  t = cbind(weather.data, weather.data.min6)
  
  if(b == 10) {
    t
  } else {
    t %>%
      rename_with(~ paste(., b, sep = "_"))
  }
  
}

weather_all = pblapply(c(5,10,20), function(x) weather(x)) %>%
  bind_cols()

## ----------------- ##
## FINAL DATA MERGER ##
## ----------------- ##
all.data = coastal.PR.all %>%
  mutate(mean.ndvi = ndvi.mean,
         mangrove.cover = mangrove.cover/(pi*(buffer**2))) %>%
  bind_cols(weather_all) %>%
  filter(mangrove.cover != 0) %>%
  select(-gdlcode)

## plot prevalence vs. mangrove land cover and mangrove NDVI
plot.c = ggplot(data = all.data %>% filter(buffer == 40), 
       aes(x = mangrove.cover*100, y = pr, size = examined)) +
  geom_point() +
  xlab("Mangrove land cover [%]") + ylab("Malaria prevalence") +
  theme_bw() + 
  theme(text = element_text(size = 15))

plot.b = ggplot(data = all.data %>% filter(buffer == 3), 
       aes(x = mean.ndvi, y = pr, size = examined)) +
  geom_point() +
  xlab("Mangrove NDVI") + ylab("Malaria prevalence") +
  theme_bw() + 
  theme(text = element_text(size = 15))

# put two plots in single panel
plot.joint = plot_grid(plot.b, plot.c, labels = c('B', 'C'), 
                    label_size = 20, rel_widths =c(1,1))

plot.joint

# export joint plot
ps.options(family = "Arial")
ggsave(filename="./Figures/mangroveVars_points_v2.svg", device = "svg", width = 270, height = 135, units = "mm")
ggsave(filename="./Figures/mangroveVars_points_v2.png", device = "png", width = 270, height = 135, units = "mm")
ggsave(filename="./Figures/mangroveVars_points_v2.pdf", device = cairo_pdf, width = 270, height = 135, units = "mm")
cairo_ps(filename="./Figures/mangroveVars_points_v2.eps", width = 10.63, height = 5.315, #in inches
         pointsize = 10, fallback_resolution = 2400)
plot.joint
dev.off()

##----------------------------------------- ##
##----------------------------------------- ##
## SUMMARISING DATA ##
##----------------------------------------- ##
##----------------------------------------- ##

##-------------------------##
## illustrate missing data ##
##-------------------------##
gg_miss_upset(all.data)
gg_miss_fct(x = all.data, fct = year_start)
vis_miss(all.data, warn_large_data = F)

# which years have the most missing data for NDVI
hist(as.numeric((filter(all.data, is.na(mean.ndvi))$year_start)), nclass = 20)
# 1996 --> no NDVI data
nrow(filter(all.data, year_start == 2008, 
            is.na(mean.ndvi)))/nrow(filter(all.data, year_start == 2008))
nrow(filter(all.data, year_start == 2019,
            is.na(mean.ndvi)))/nrow(filter(all.data, year_start == 2019))
# 2007 (12 %) and 2018 (25 %) --> unknown issue

nrow(filter(all.data, is.na(mean.ndvi)))/nrow(filter(all.data))
nrow(filter(all.data, is.na(mangrove.cover.min1)))/nrow(filter(all.data))
# missing data: mangrove.cover.min1 - 19%
#               mean.ndvi - 8%

##----------------------------------------------##
## export dataset without missing mangrove data ##
##   --> for robustness analysis                ##
##----------------------------------------------##
all.data.small = all.data %>%
  filter(!is.na(mangrove.cover.min1), !is.na(mean.ndvi))
pc.small = preProcess(all.data.small %>% select(-method),
                method = c("center", "scale"))
alldata.small.2 = predict(pc.small, all.data.small %>%
                            select(-method)) %>%
  # revert to original values (not centred or scaled) for some columns
  mutate(lat = all.data.small$lat,
         lon = all.data.small$lon,
         buffer = all.data.small$buffer, 
         pr = all.data.small$pr, 
         examined = all.data.small$examined,
         year_start = all.data.small$year_start,
         month_start = all.data.small$month_start,
         year_end = all.data.small$year_end,
         month_end = all.data.small$month_end) %>%
  cbind(., all.data.small %>% select(method))

write.table(alldata.small.2, "./data/alldata.small.csv", sep = ",", row.names = F)

##----------------------------------------- ##
##----------------------------------------- ##
## KNN imputation                           ##
##       WARNING: process might take long!! ##
##----------------------------------------- ##
##----------------------------------------- ##
pc = preProcess(all.data %>% select(-method),
                method = "knnImpute", k = 10)
alldata.impute.1 = predict(pc, all.data %>% select(-method))

alldata.impute.2 = alldata.impute.1 %>%
  # revert to original values (not centred or scaled) for some columns
  mutate(lat = all.data$lat,
         lon = all.data$lon,
         buffer = all.data$buffer, 
         pr = all.data$pr, 
         examined = all.data$examined,
         year_start = all.data$year_start,
         month_start = all.data$month_start,
         year_end = all.data$year_end,
         month_end = all.data$month_end) %>%
  cbind(., all.data %>% select(method))

## export final dataset ##
write.table(alldata.impute.2, "./data/alldata.impute.csv", sep = ",", row.names = F)
