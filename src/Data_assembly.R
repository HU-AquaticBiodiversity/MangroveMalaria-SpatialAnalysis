# ==============================================================================
# DATA PREPARATION & SPATIAL EXTRACTION PIPELINE
# Description: Merges point prevalence data with multi-scale spatial rasters 
# (weather, NDVI, human impact, mosquito distributions) using parallelised 
# extraction. Prepares the final imputed dataset for SEM modeling.
# ==============================================================================

# --- 1. Load Libraries ---
# Grouped by functionality for better readability
# Data Manipulation & I/O
library(dplyr)
library(tidyverse)
library(data.table)
library(vroom)
library(stringr)
library(purrr)
library(lubridate)

# Spatial & GIS
library(sf)
library(sp)
library(lwgeom)
library(geosphere)
library(osmdata)
library(units)
library(terra)
library(tidyterra)
library(exactextractr)
library(rnaturalearth)

# Parallel Processing
library(parallel)
library(doParallel)
library(pbapply)

# Machine Learning & Missing Data
library(caret)
library(naniar)

# Visualization
library(ggplot2)
library(viridis)
library(cowplot)


# --- 2. Load Core Datasets ---
country.list <- read.csv("./data/country_table.csv")

# Load prevalence data and bind DHS data
coastal.PR <- read.csv("./data/coastal.PR.final.csv", sep = ",", header = TRUE) %>%
  mutate(site_id = as.factor(site_id)) %>%
  bind_rows(read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = TRUE))

# Load coastline data
coastline <- ne_coastline(scale = 10, returnclass = "sf")

# Limit countries to those with available prevalence data
countries.available <- country.list %>% filter(ISO3 %in% coastal.PR$country_id)

# Load malaria polygons with 5, 10, and 20-km radii (Result: 4260 locations)
all_polygons.10 <- st_read("./data/polygons_pr_malaria.shp") %>%
  bind_rows(st_read("./data/polygons_pr_malaria_dhs.shp")) %>%
  filter(buffer %in% c(5, 10, 20))


# ==============================================================================
# A. CONSTANT RASTER EXTRACTION (Single-Year Variables)
# ==============================================================================


# 1. Mosquito Distributions (Saltwater species)
prob.list <- sprc(
  rast('./res/mosquito_distributions/2017_Anopheles_melas.Mean_Decompressed.geotiff'),
  rast('./res/mosquito_distributions/2017_Anopheles_merus.Mean_Decompressed.geotiff')
)
prob.salt <- terra::mosaic(prob.list, fun = "sum")

# 2. Health Management Data (Travel times)
t_city <- rast('./res/Health_management/201501_Global_Travel_Time_to_Cities_2015/201501_Global_Travel_Time_to_Cities_2015.tif')
t_healthcare_motor <- rast('./res/Health_management/202001_Global_Motorized_Travel_Time_to_Healthcare_2019/202001_Global_Motorized_Travel_Time_to_Healthcare_2019.tif')
t_healthcare_walk <- rast('./res/Health_management/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.tif')

# Extract median and SD values for all constant rasters (10km buffer only)
const.raster.data <- pblapply(
  list(prob.salt, t_city, t_healthcare_motor, t_healthcare_walk),
  FUN = function(y) {
    as.data.frame(exact_extract(y, all_polygons.10, fun = c('median', 'stdev')))
  }
) %>% 
  bind_cols() %>%
  mutate(lon = all_polygons.10$lon,
         lat = all_polygons.10$lat,
         buffer = all_polygons.10$buffer) %>% 
  filter(buffer == 10) %>% 
  dplyr::select(-buffer)

# Rename columns
names(const.raster.data) <- c(
  'mosq.median', 'mosq.sd', 't_city.median', 't_city.sd', 
  't_healthcare_motor.median', 't_healthcare_motor.sd',
  't_healthcare_walk.median', 't_healthcare_walk.sd',
  'lon', 'lat'
)

write.table(const.raster.data, "./data/const.raster.data.csv", sep = ",", row.names = FALSE)


# ==============================================================================
# B. ANNUAL RASTER EXTRACTION (Multi-Year Variables)
# ==============================================================================

# Load all annual raster files
ITN_use <- rast(list.files('./res/health_management/2024_GBD2023_Africa_ITN_2000/', pattern = "2024_GBD2023", full.names = TRUE))
ITN_access <- rast(list.files('./res/health_management/ITN_2000_access_mean/', pattern = "ITN", full.names = TRUE))
IRS <- rast(list.files('./res/health_management/2024_GBD2023_Africa_IRS_2000/', pattern = "2024_GBD2023", full.names = TRUE))
AMT <- rast(list.files('./res/health_management/2024_GBD2023_Global_Antimalarial_EFT_2000/', pattern = "2024_GBD2023", full.names = TRUE))
POP <- rast(list.files('./res/population_density/', full.names = TRUE))

# Extract median and SD values for annual rasters
annual.raster.data <- lapply(list(ITN_use, ITN_access, IRS, AMT, POP), function(x) {
  exact_extract(x, all_polygons.10, fun = c('median', 'stdev'))
}) %>% 
  bind_cols() %>%
  mutate(lon = all_polygons.10$lon, lat = all_polygons.10$lat, buffer = all_polygons.10$buffer) %>%
  pivot_longer(!lon & !lat & !buffer, names_to = 'dummy', values_to = 'value') %>%
  mutate(
    # Extract metadata from filenames via regex
    year_start = as.numeric(ifelse(str_ends(dummy,"\\d"), str_sub(dummy, -4), str_extract(dummy, "(\\d)+"))),
    statistic = str_extract(dummy, "[^.]+"),
    name = str_remove_all(str_split(dummy, fixed("."), simplify = TRUE)[,2], "(\\d)+"),
    name = ifelse(name %in% c("Af_more_qr_clust_", "Ex_AF_Cot_NoNeg_"), "pop_dens", name),
    name = str_remove(str_remove(str_replace(name, '__','_'), "^\\_"), "\\_$")
  ) %>%
  dplyr::select(-dummy) %>%
  pivot_wider(names_from = c(name, statistic), values_from = value)

write.table(annual.raster.data, "./data/annual.raster.data.csv", sep = ",", row.names = FALSE)


# ==============================================================================
# C. MONTHLY RASTER EXTRACTION (Cropland & Weather)
# ==============================================================================

# 1. Cropland (Split by East/West Africa to match raster extents)
layer.names <- paste0("lc_", str_extract(list.files('./res/landcover/East/', full.names = TRUE), pattern = "\\d{4}"))

lc_east <- rast(list.files('./res/landcover/East/', full.names = TRUE), subds = 'lccs_class')
lc_west <- rast(list.files('./res/landcover/West/', full.names = TRUE), subds = 'lccs_class')

# Re-classify: cropland = 1, other = NA
lc_east_bin <- ifel(lc_east %in% c(10,20,30,40), 1, NA)
lc_west_bin <- ifel(lc_west %in% c(10,20,30,40), 1, NA)

poly_west <- all_polygons.10[all_polygons.10$lon < 28,]
poly_east <- all_polygons.10[all_polygons.10$lon > 28,]

# Extract surface area of croplands in km²
lc_west_list <- lapply(1:nlyr(lc_west_bin), function(x) {
  exact_extract(cellSize(lc_west_bin[[x]], unit = 'km', mask = TRUE), poly_west, fun = 'sum', progress = TRUE)
})
lc_east_list <- lapply(1:nlyr(lc_east_bin), function(x) {
  exact_extract(cellSize(lc_east_bin[[x]], unit = 'km', mask = TRUE), poly_east, fun = 'sum', progress = TRUE)
})

lc_west.tab <- do.call(cbind, lc_west_list)
lc_east.tab <- do.call(cbind, lc_east_list)
colnames(lc_west.tab) <- layer.names
colnames(lc_east.tab) <- layer.names

lc_all <- st_drop_geometry(rbind(cbind(poly_west, lc_west.tab), cbind(poly_east, lc_east.tab)))
write.csv(lc_all, "./data/cropland.csv", row.names = FALSE)

# 2. Weather Data (Projected to WGS84 for alignment)
lyr.names <- str_sub(list.files("./res/climate/"), end = -10)
weather.raster <- terra::project(
  rast(list.files("./res/climate/", full.names = TRUE)),
  '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
  use_gdal = TRUE, threads = 7
)
names(weather.raster) <- lyr.names

weather.data <- exact_extract(weather.raster, all_polygons.10, fun = 'mean') %>% 
  bind_cols(all_polygons.10 %>% st_drop_geometry() %>% dplyr::select(lon, lat, buffer))

write.csv(weather.data, "./data/weather_data.csv", row.names = FALSE)


# ==============================================================================
# D. MASTER DATA MERGE
# ==============================================================================

# HDI Data
nHDI.table <- read.csv("./res/GDL_Shapefiles_V6.3/GDL-Subnational-HDI-data.csv") %>%
  pivot_longer(X1990:X2021, names_to = 'year', values_to = 'nHDI') %>%
  mutate(year = as.numeric(str_sub(year, start = 2))) %>%
  dplyr::select(GDLCODE, year, nHDI)
write.table(nHDI.table, "./data/nHDI.data.csv", sep = ",", row.names = FALSE)

# Shapefiles for GDL (Using st_make_valid to fix topological errors)
shp.subnat <- st_read("./res/GDL_Shapefiles_V6.3/GDL Shapefiles V6.3 large.shp") %>%
  dplyr::select(-iso_code) %>%
  st_make_valid()

# Format Cropland and Mangrove data for the merge
lc_all_formatted <- read.csv("./data/cropland.csv") %>%
  pivot_longer(cols = starts_with("lc_"), names_to = "year", values_to = "lc_crop", values_drop_na = TRUE) %>%
  mutate(
    year = as.numeric(str_extract(year, pattern = "\\d{4}")),
    buffer = ifelse(buffer == 10, "lc_crop", paste0("lc_crop_", buffer))
  ) %>%
  pivot_wider(names_from = "buffer", values_from = "lc_crop")

mangrove.cover <- vroom(list.files("data/MangroveCover", full.names = TRUE), id = 'year') %>%
  mutate(year = as.numeric(str_extract(year, pattern = "\\d{4}")))

# Execute Master Join Pipeline
coastal.PR.all <- coastal.PR %>%
  dplyr::select(pr, lat, lon, examined, method, coastline.dist, month_start, year_start, month_end, year_end) %>%
  st_as_sf(coords = c('lon', 'lat'), remove = FALSE, crs = 4326) %>%
  st_join(shp.subnat) %>%
  st_drop_geometry() %>%
  left_join(read.csv("./data/const.raster.data.csv"), by = c("lon", "lat")) %>%
  left_join(lc_all_formatted, by = c("year_start" = "year", 'lon', 'lat')) %>%
  left_join(
    read.csv("./data/annual.raster.data.csv") %>%
      filter(buffer == 10) %>%
      dplyr::select(-buffer) %>%
      mutate(
        pop_dens_median_5 = read.csv("./data/annual.raster.data.csv") %>% filter(buffer == 5) %>% pull(pop_dens_median),
        pop_dens_median_20 = read.csv("./data/annual.raster.data.csv") %>% filter(buffer == 20) %>% pull(pop_dens_median)
      ),
    by = c('year_start', 'lon', 'lat')
  ) %>%
  right_join(mangrove.cover, by = c("year_start" = "year", 'lon', 'lat'), relationship = "many-to-many") %>%
  filter(!is.na(pr)) %>%
  left_join(
    mangrove.cover %>% rename(mangrove.cover.min1 = mangrove.cover) %>% mutate(year = year + 1),
    by = c("year_start" = "year", 'lon', 'lat', 'buffer'),
    relationship = "many-to-many"
  ) %>%
  left_join(nHDI.table, by = c("gdlcode" = "GDLCODE", "year_start" = "year")) %>%
  dplyr::select(-continent)


# ==============================================================================
# E. TARGETED TEMPORAL EXTRACTIONS (NDVI & Weather)
# ==============================================================================
n.cores <- detectCores()

# 1. Parallel NDVI Extraction (Matches exact survey dates)
ndvi.data <- list.files("./data/NDVI.data/", full.names = TRUE) %>% map_df(~read_csv(.))

clust <- makeCluster(n.cores - 1)
clusterExport(clust, c("coastal.PR.all", "ndvi.data"))

system.time({
  ndvi.mean <- parSapply(clust, 1:nrow(coastal.PR.all), function(x) {
    library(dplyr)
    library(lubridate)
    
    a <- ndvi.data %>%
      filter(buffer == coastal.PR.all[x,]$buffer & lat == coastal.PR.all[x,]$lat & lon == coastal.PR.all[x,]$lon) %>%
      filter(between(
        as.Date(month), 
        as.Date(paste(coastal.PR.all[x,]$year_start, coastal.PR.all[x,]$month_start, "1", sep = "-")),
        as.Date(paste(coastal.PR.all[x,]$year_end, coastal.PR.all[x,]$month_end, "1", sep = "-")) %m+% months(1)
      ))
    mean(a$mean)
  })
})
stopCluster(clust)

# 2. Parallel Weather Extraction Function (Survey period & 6 months prior)
weather_extract <- function(b, cores = n.cores) {
  
  weather.table <- read.csv("./data/weather_data.csv") %>%
    filter(buffer == b) %>%
    dplyr::select(-buffer) %>%
    pivot_longer(cols = -c(lat, lon), names_to = "measure", values_to = "value") %>%
    mutate(
      measure = str_replace_all(measure, "_1991\\.2020|mean\\.1month_|_Global_ea", ""),
      date = paste(str_sub(str_extract(measure, "\\d{6}"), 1, 4), str_sub(str_extract(measure, "\\d{6}"), 5, 6), "01", sep = "-"),
      measure = str_remove(measure, "_\\d{6}")
    ) %>%
    pivot_wider(names_from = measure, values_from = value)
  
  # A. Run for weather WITHIN survey period
  clust <- makeCluster(cores - 1)
  clusterExport(clust, c("coastal.PR.all", "weather.table"), envir = environment())
  
  weather.data <- parLapply(clust, 1:nrow(coastal.PR.all), function(x) {
    library(dplyr)
    library(lubridate)
    
    weather.table %>%
      filter(lat == coastal.PR.all[x,]$lat & lon == coastal.PR.all[x,]$lon) %>%
      filter(between(
        as.Date(date), 
        as.Date(paste(coastal.PR.all[x,]$year_start, coastal.PR.all[x,]$month_start, "1", sep = "-")),
        as.Date(paste(coastal.PR.all[x,]$year_end, coastal.PR.all[x,]$month_end, "1", sep = "-")) %m+% months(1)
      )) %>%
      dplyr::select(anomaly_2t:mean_tp) %>%
      colMeans()
  }) %>% bind_rows()
  stopCluster(clust)
  
  # B. Run for weather 6 MONTHS PRIOR to survey period
  clust <- makeCluster(cores - 1)
  clusterExport(clust, c("coastal.PR.all", "weather.table"), envir = environment())
  
  weather.data.min6 <- parLapply(clust, 1:nrow(coastal.PR.all), function(x) {
    library(dplyr)
    library(lubridate)
    
    weather.table %>%
      filter(lat == coastal.PR.all[x,]$lat & lon == coastal.PR.all[x,]$lon) %>%
      filter(between(
        as.Date(date), 
        as.Date(paste(coastal.PR.all[x,]$year_start, coastal.PR.all[x,]$month_start, "1", sep = "-")) %m-% months(6),
        as.Date(paste(coastal.PR.all[x,]$year_end, coastal.PR.all[x,]$month_end, "1", sep = "-"))
      )) %>%
      dplyr::select(anomaly_2t:mean_tp) %>%
      colSums()
  }) %>% bind_rows() %>%
    rename_with(~ paste(., "6m", sep = "_"))
  stopCluster(clust)
  
  # Combine and rename based on buffer size
  t <- cbind(weather.data, weather.data.min6)
  if(b == 10) return(t) else return(t %>% rename_with(~ paste(., b, sep = "_")))
}

# Run extraction across 5, 10, and 20 km buffers
weather_all <- pblapply(c(5, 10, 20), function(x) weather_extract(x)) %>% bind_cols()


# ==============================================================================
# F. FINAL MERGE & MISSING DATA IMPUTATION
# ==============================================================================

# Create Full Dataset
all.data <- coastal.PR.all %>%
  mutate(mean.ndvi = ndvi.mean,
         mangrove.cover = mangrove.cover / (pi * (buffer^2))) %>%
  bind_cols(weather_all) %>%
  filter(mangrove.cover != 0) %>%
  dplyr::select(-gdlcode)

write.table(all.data, "./data/alldata.full.csv", sep = ",", row.names = FALSE)

# (Optional: Plotting Code Omitted for Brevity - Works as written)


# --- Create Complete Case (Small) Dataset ---
all.data.small <- all.data %>% filter(!is.na(mangrove.cover.min1) & !is.na(mean.ndvi))

pc.small <- preProcess(all.data.small %>% dplyr::select(-method), method = c("center", "scale"))
alldata.small.2 <- predict(pc.small, all.data.small %>% dplyr::select(-method)) %>%
  mutate(
    lat = all.data.small$lat, lon = all.data.small$lon, buffer = all.data.small$buffer, 
    pr = all.data.small$pr, examined = all.data.small$examined,
    year_start = all.data.small$year_start, month_start = all.data.small$month_start,
    year_end = all.data.small$year_end, month_end = all.data.small$month_end
  ) %>%
  bind_cols(all.data.small %>% dplyr::select(method))

write.table(alldata.small.2, "./data/alldata.small.csv", sep = ",", row.names = FALSE)

# --- Perform KNN Imputation on Full Dataset ---
all.data.10 <- all.data %>% dplyr::select(-ends_with("_5"), -ends_with("_20"), -method)
all.data.520 <- all.data %>% dplyr::select(ends_with("_5"), ends_with("_20"), -method)

set.seed(123)
# Impute missing data for 10-km dataset
pc.10 <- preProcess(all.data.10, method = "knnImpute", k = 10)
impute.10 <- predict(pc.10, all.data.10)

# Center and scale 5 and 20-km variables (No imputation)
pc.520 <- preProcess(all.data.520, method = c("center", "scale"))
alldata.proc.520 <- predict(pc.520, all.data.520)

alldata.impute <- impute.10 %>%
  mutate(
    lat = all.data$lat, lon = all.data$lon, buffer = all.data$buffer,
    pr = all.data$pr, examined = all.data$examined,
    year_start = all.data$year_start, month_start = all.data$month_start,
    year_end = all.data$year_end, month_end = all.data$month_end
  ) %>%
  bind_cols(alldata.proc.520) %>%
  bind_cols(all.data %>% dplyr::select(method))

write.table(alldata.impute, "./data/alldata.impute.csv", sep = ",", row.names = FALSE)

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
