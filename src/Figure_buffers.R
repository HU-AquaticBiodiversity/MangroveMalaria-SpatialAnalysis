# ==============================================================================
# Script Name: Spatial Variable Visualisation (3-Panel Map)
# Description: Generates a high-resolution, 3-panel figure illustrating 
#              mangrove land cover, mangrove NDVI (vegetation health), and 
#              population density within specific spatial buffer zones (1-50 km) 
#              around a selected site in the Saloum Delta, Senegal.
# ==============================================================================

## =============================================================================
## LOAD NECESSARY LIBRARIES
## =============================================================================
# Spatial data handling (Vector & Raster)
library(sf)             # For vector data handling (polygons, points)
library(sp)             # Legacy spatial classes
library(lwgeom)         # Advanced geometry operations
library(geosphere)      # Spherical trigonometry for distance/bearing
library(terra)          # Modern, fast package for raster data handling
library(exactextractr)  # Fast extraction of raster values by polygons

# Data wrangling & API access
library(dplyr)          # Data manipulation
library(tidyverse)      # Core tidyverse packages
library(data.table)     # Fast data table manipulation
library(osmdata)        # OpenStreetMap data retrieval
library(rnaturalearth)  # For base maps (coastlines, landmasses)
library(pbapply)        # For progress bars on apply loops
library(units)          # Measurement units handling
library(ape)            # Phylogenetic and evolutionary analysis

# Visualisation & Plotting
library(ggplot2)        # Core plotting engine
library(tidyterra)      # Allows ggplot to plot terra SpatRasters directly
library(viridis)        # Colorblind-friendly palettes (scale_fill_viridis)
library(cowplot)        # For combining multiple plots into a grid
library(ggrepel)        # For repelling overlapping text labels

## =============================================================================
## 1. DATA PREPARATION (GLOBAL)
## =============================================================================

# Load malaria prevalence buffer polygons and filter to the maximum 50 km radius.
# This serves as the outer boundary for our spatial extractions.
coastal.PR.polygons = st_read("./data/polygons_pr_malaria.shp")
coastal.PR.polygons.50 = coastal.PR.polygons %>%
  filter(buffer == 50)

# Load the raw prevalence coordinate data table
coastal.PR = read.csv("./data/coastal.PR.final.csv", sep = ",", header = TRUE)

## =============================================================================
## 2. SET UP EXAMPLE EXTENT (SALOUM DELTA, SENEGAL)
## =============================================================================

# Isolate the specific survey location in Senegal by filtering latitude
polygon.ex = coastal.PR.polygons %>%
  filter(lat > 14 & lat < 14.1) %>%
  # Dynamically calculate the precise geographic coordinates for label placement
  mutate(
    # X coordinate: The far-right edge of the buffer (East) for placing labels
    label_x = sapply(st_geometry(.), function(p) st_bbox(p)["xmax"]),
    # Y coordinate: The vertical center of the buffer for vertical alignment
    label_y = sapply(st_geometry(.), function(p) mean(c(st_bbox(p)["ymin"], st_bbox(p)["ymax"])))
  )

## =============================================================================
## 3. PANEL A: MANGROVE EXTENT MAP (p.mangrove)
## =============================================================================


# Locate and load the Global Mangrove Watch (GMW) raster tiles for the year 2008
tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", 2008),
                   pattern = ".tif", full.names = TRUE)
tile.list = lapply(tiles, rast) # Convert file paths to terra SpatRaster objects

# Optimisation: Identify ONLY the raster tiles that intersect with our 50km polygon
tile.overlaps = pblapply(1:length(tile.list), function(x) {
  nrow(st_crop(polygon.ex[50,], tile.list[[x]])) != 0
}) %>% unlist()

# Mosaic (merge) only the required overlapping tiles into a single seamless raster
tile.sel = mosaic(sprc(tile.list[tile.overlaps]))

# Load Natural Earth land polygon for the background reference map
ne.land = st_read("./res/ne_10m_land/ne_10m_land.shp")

# Build Panel A: Mangrove Cover Map
p.mangrove = ggplot() +
  # 1. Base map: Fill landmasses with a distinct background color
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  # 2. Raster data: Plot the mangrove extent
  geom_spatraster(data = tile.sel, aes(fill = GMW_N14W017_2008_v3)) +
  scale_fill_gradient(low = "grey", high = "green4", na.value = "transparent") +
  # 3. Vector data: Plot invisible buffer rings (used purely for mapping label coordinates)
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)), alpha = 0) +
  # 4. Labels: Tag the buffer rings dynamically
  geom_label_repel(
    data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)),
    aes(x = label_x, y = label_y, label = paste(buffer, "km")),
    direction = "y",          # Forces labels to arrange strictly vertically
    nudge_x = 0.05,           # Pushes the stack of labels to the right into a neat column
    hjust = 0,                # Left-aligns the text in the column
    segment.size = 0.3,       # Thickness of the line pointing to the circle
    segment.color = "grey30", # Color of the pointing line
    size = 2.5,
    fill = alpha("white", 0.7), 
    label.size = 0,
    seed = 42                 # Ensures the layout looks identical on every run
  ) +
  # 5. Cropping: Zoom tightly into the Saloum Delta bounding box
  coord_sf(xlim = c(-17.2, -16.1), ylim = c(13.5, 14.6), expand = FALSE) +
  # 6. Aesthetics & Theming
  xlab("longitude") + ylab("latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'skyblue', colour = 'red'),
        legend.position = 'none',
        text = element_text(family = "sans")) +
  ggtitle("Mangrove land cover")

## =============================================================================
## 4. PANEL B: NDVI (VEGETATION HEALTH) MAP (p.ndvi)
## =============================================================================


# Locate and load the 2008 NDVI raster data for West Africa
ndvi.files.w = list.files(
  "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/West",
  full.names = TRUE, 
  pattern = paste0("*NDVI_doy", 2008,".*tif"))

NDVI.w = lapply(ndvi.files.w, rast)

# Create a spatial mask: Convert the cropped mangrove raster into a vector polygon
mpoly.ex.50 = as.polygons(crop(tile.sel, polygon.ex[50,], mask = TRUE))

# Mask the NDVI data so we ONLY visualise vegetation health *inside* the mangrove boundaries
NDVI.crop = crop(NDVI.w[[15]], mpoly.ex.50, mask = TRUE)

# Build Panel B: Mangrove NDVI Map
p.ndvi = ggplot() +
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  # Plot the strictly masked NDVI raster
  geom_spatraster(data = NDVI.crop, 
                  aes(fill = MOD13Q1.061__250m_16_days_NDVI_doy2008225_aid0001)) +
  scale_fill_gradient(low = "yellow", high = "green", na.value = "transparent") +
  # Add invisible buffer rings for labeling
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)), alpha = 0) +
  geom_label_repel(
    data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)),
    aes(x = label_x, y = label_y, label = paste(buffer, "km")),
    direction = "y",       
    nudge_x = 0.05,        
    hjust = 0,             
    segment.size = 0.3,    
    segment.color = "grey30", 
    size = 2.5,
    fill = alpha("white", 0.7), 
    label.size = 0,
    seed = 42              
  ) +
  coord_sf(xlim = c(-17.2, -16.1), ylim = c(13.5, 14.6), expand = FALSE) +
  xlab("longitude") + ylab("latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'skyblue', colour = 'red'),
        text = element_text(family = "sans")) +
  labs(fill = "NDVI") +
  ggtitle("Mangrove NDVI")

## =============================================================================
## 5. PANEL C: POPULATION DENSITY MAP (p.pop)
## =============================================================================


# Load global population raster and crop roughly to our specific coordinate bounds
POP.sel = rast("./res/population_density/GlobPOP_Count_30arc_2008_I32.tiff")
POP.sel.crop = crop(POP.sel, ext(-17.2, -16.1, 13.5, 14.6))

# Build Panel C: Population Density Map
p.pop = ggplot() +
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  # Plot the population density raster
  geom_spatraster(data = POP.sel.crop) +
  # Log-transform the color scale. This is crucial so highly dense urban clusters 
  # don't visually wash out the sparser rural populations.
  scale_fill_viridis(na.value = "transparent",
                     option = "plasma",
                     trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # Add a visible solid border specifically for the 10km buffer
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(10)), 
          alpha = 0, colour = "black") +
  # Add visible dashed borders for the 5km and 20km buffers to highlight key zones
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(5,20)), 
          alpha = 0, linetype = "dashed", colour = "black") +
  # Label the 5, 10, and 20km buffers
  geom_label_repel(
    data = polygon.ex %>% filter(buffer %in% c(5,10,20)),
    aes(x = label_x, y = label_y, label = paste(buffer, "km")),
    direction = "y",       
    nudge_x = 0.05,        
    hjust = 0,             
    segment.size = 0.3,    
    segment.color = "grey30", 
    size = 2.5,
    fill = alpha("white", 0.7), 
    label.size = 0,
    seed = 42              
  ) +
  coord_sf(xlim = c(-17.2, -16.1), ylim = c(13.5, 14.6), expand = FALSE) +
  xlab("longitude") + ylab("latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'skyblue', colour = 'red'),
        text = element_text(family = "sans")) +
  # Expression used to render scientific notation nicely in the legend title
  labs(fill = expression("Population\ndensity [km"^{-2}* "]")) +
  ggtitle("Fixed variables (e.g. population density)")

## =============================================================================
## 6. ASSEMBLE AND EXPORT FINAL FIGURE
## =============================================================================

# Combine the three individual plots into a cohesive figure using cowplot.
# Top row: Mangrove (A) and NDVI (B) plotted side-by-side.
# Bottom row: Population map (C) spanning the full width beneath them.
p.joint = plot_grid(
  plot_grid(p.mangrove, p.ndvi, labels = c('A', 'B'), 
            label_size = 20, rel_widths = c(46.3, 53.7)),
  p.pop, 
  labels = c('', 'C'), # Empty string ensures 'C' aligns correctly with the bottom plot
  label_size = 20, rel_heights = c(1, 1),
  nrow = 2 # Forces p.pop onto the second row
)

# Export the combined figure in multiple high-resolution formats for publication
ggsave(filename="./Figures/mangroveVars_v2.svg", plot = p.joint, device = "svg", width = 270, height = 240, units = "mm")
ggsave(filename="./Figures/mangroveVars_v2.png", plot = p.joint, device = "png", width = 270, height = 240, units = "mm")
ggsave(filename="./Figures/mangroveVars_v2.pdf", plot = p.joint, device = cairo_pdf, width = 270, height = 240, units = "mm")

# Special handling for EPS export (often required by journals)
cairo_ps(filename="./Figures/mangroveVars_v2.eps", width = 10.63, height = 9.45, 
         pointsize = 10, fallback_resolution = 2400)
print(p.joint)
dev.off()