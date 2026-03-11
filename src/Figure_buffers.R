## =============================================================================
## LOAD NECESSARY LIBRARIES
## =============================================================================
library(dplyr)
library(sf)             # For vector data handling (polygons, points)
library(sp)
library(rnaturalearth)  # For base maps (coastlines, landmasses)
library(ggplot2)
library(lwgeom)
library(tidyverse)
library(geosphere)
library(osmdata)
library(units)
library(data.table)
library(exactextractr)
library(terra)          # Modern package for raster data handling
library(pbapply)        # For progress bars on apply functions
library(ape)
library(tidyterra)      # Allows ggplot to plot terra SpatRasters directly
library(viridis)        # ADDED: Needed for scale_fill_viridis() in p.pop
library(cowplot)        # ADDED: Needed for plot_grid() at the end
library(ggrepel)

## =============================================================================
## 1. DATA PREPARATION (GLOBAL)
## =============================================================================

# Load malaria prevalence polygons and filter to the maximum 50km radius
coastal.PR.polygons = st_read("./data/polygons_pr_malaria.shp")
coastal.PR.polygons.50 = coastal.PR.polygons %>%
  filter(buffer == 50)

# Load the raw prevalence data table
coastal.PR = read.csv("./data/coastal.PR.final.csv", sep = ",", header = TRUE)


## =============================================================================
## 2. SET UP EXAMPLE EXTENT (SALOUM DELTA, SENEGAL)
## =============================================================================

polygon.ex = coastal.PR.polygons %>%
  filter(lat > 14 & lat < 14.1) %>%
  # Calculate the exact geographic coordinates for the right edge of each circle
  mutate(
    # X coordinate: The far-right edge of the buffer (East)
    label_x = sapply(st_geometry(.), function(p) st_bbox(p)["xmax"]),
    # Y coordinate: The vertical center of the buffer
    label_y = sapply(st_geometry(.), function(p) mean(c(st_bbox(p)["ymin"], st_bbox(p)["ymax"])))
  )


## =============================================================================
## 3. PANEL A: MANGROVE EXTENT MAP (p.mangrove)
## =============================================================================

# Locate and load the Global Mangrove Watch (GMW) raster tiles for 2008
tiles = list.files(paste0("./res/mangrove_raster/gmw_v3_", 2008),
                   pattern = ".tif", full.names = TRUE)
tile.list = lapply(tiles, rast) # Convert file paths to terra SpatRaster objects

# Identify which raster tiles intersect with our 50km polygon to avoid loading extra data
tile.overlaps = pblapply(1:length(tile.list), function(x) {
  nrow(st_crop(polygon.ex[50,], tile.list[[x]])) != 0
}) %>% unlist()

# Mosaic (merge) only the overlapping tiles into a single raster
tile.sel = mosaic(sprc(tile.list[tile.overlaps]))

# Load natural earth land polygon for the background map
ne.land = st_read("./res/ne_10m_land/ne_10m_land.shp")

# Build Panel A: Mangrove Cover
p.mangrove = ggplot() +
  # 1. Base map: Landmasses
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  # 2. Raster data: Mangrove extent
  geom_spatraster(data = tile.sel, aes(fill = GMW_N14W017_2008_v3)) +
  scale_fill_gradient(low = "grey", high = "green4", na.value = "transparent") +
  # 3. Vector data: Buffer rings (invisible borders, used only for mapping coordinates)
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)), alpha = 0) +
  # 4. Labels: Tag the buffer rings
  geom_label_repel(
    data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)),
    aes(x = label_x, y = label_y, label = paste(buffer, "km")),
    direction = "y",       # FIX: Forces labels to arrange strictly vertically
    nudge_x = 0.05,        # Pushes the stack of labels to the right into a neat column
    hjust = 0,             # Left-aligns the text in the column
    segment.size = 0.3,    # Thickness of the line pointing to the circle
    segment.color = "grey30", # Color of the pointing line
    size = 2.5,
    fill = alpha("white", 0.7), 
    label.size = 0,
    seed = 42              # Ensures the layout looks the exact same every time you run it
  ) +
  # 5. Cropping: Zoom tightly into the Saloum Delta
  coord_sf(xlim = c(-17.2, -16.1), ylim = c(13.5, 14.6), expand = FALSE) +
  # 6. Aesthetics
  xlab("longitude") + ylab("latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'skyblue', colour = 'red'),
        legend.position = 'none',
        text = element_text(family = "sans")) +
  ggtitle("Mangrove land cover")


## =============================================================================
## 4. PANEL B: NDVI (VEGETATION HEALTH) MAP (p.ndvi)
## =============================================================================

# Locate and load the 2008 NDVI raster data
ndvi.files.w = list.files(
  "/vsc-hard-mounts/leuven-data/347/vsc34705/Mangrove_PL/VI_raster/West",
  full.names = TRUE, 
  pattern = paste0("*NDVI_doy", 2008,".*tif"))

NDVI.w = lapply(ndvi.files.w, rast)

# Create a spatial mask: Convert the 50km cropped mangrove area into a polygon
mpoly.ex.50 = as.polygons(crop(tile.sel, polygon.ex[50,], mask = TRUE))

# Mask the NDVI data so we ONLY see vegetation health *inside* the mangroves
NDVI.crop = crop(NDVI.w[[15]], mpoly.ex.50, mask = TRUE)

# Build Panel B: NDVI
p.ndvi = ggplot() +
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  # Plot the masked NDVI raster
  geom_spatraster(data = NDVI.crop, 
                  aes(fill = MOD13Q1.061__250m_16_days_NDVI_doy2008225_aid0001)) +
  scale_fill_gradient(low = "yellow", high = "green", na.value = "transparent") +
  # Add invisible buffer rings for labeling
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)), alpha = 0) +
  geom_label_repel(
    data = polygon.ex %>% filter(buffer %in% c(1,10,20,30,40,50)),
    aes(x = label_x, y = label_y, label = paste(buffer, "km")),
    direction = "y",       # FIX: Forces labels to arrange strictly vertically
    nudge_x = 0.05,        # Pushes the stack of labels to the right into a neat column
    hjust = 0,             # Left-aligns the text in the column
    segment.size = 0.3,    # Thickness of the line pointing to the circle
    segment.color = "grey30", # Color of the pointing line
    size = 2.5,
    fill = alpha("white", 0.7), 
    label.size = 0,
    seed = 42              # Ensures the layout looks the exact same every time you run it
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

# Load global population raster and crop roughly to our coordinate bounds
POP.sel = rast("./res/population_density/GlobPOP_Count_30arc_2008_I32.tiff")
POP.sel.crop = crop(POP.sel, ext(-17.2, -16.1, 13.5, 14.6))

# Build Panel C: Population Density
p.pop = ggplot() +
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  # Plot the population raster
  geom_spatraster(data = POP.sel.crop) +
  # Log transform the color scale so highly dense areas don't wash out rural areas
  scale_fill_viridis(na.value = "transparent",
                     option = "plasma",
                     trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # Add visible solid border for 10km buffer
  geom_sf(data = polygon.ex %>% filter(buffer %in% c(10)), 
          alpha = 0, colour = "black") +
  # Add visible dashed borders for 5km and 20km buffers
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
  labs(fill = expression("Population\ndensity [km"^{-2}* "]")) +
  ggtitle("Fixed variables (e.g. population density)")


## =============================================================================
## 6. ASSEMBLE AND EXPORT FINAL FIGURE
## =============================================================================

# Combine the three plots using cowplot::plot_grid
# Top row: Mangrove (A) and NDVI (B). Bottom row: Population (spanning full width)
p.joint = plot_grid(
  plot_grid(p.mangrove, p.ndvi, labels = c('A', 'B'), 
            label_size = 20, rel_widths = c(46.3, 53.7)),
  p.pop, 
  labels = c('', 'C'), # Added a label for the bottom plot to match A and B
  label_size = 20, rel_heights = c(1, 1),
  nrow = 2 # Forces p.pop to be on the second row
)

# Export the figure in multiple high-resolution formats
ggsave(filename="./Figures/mangroveVars_v2.svg", plot = p.joint, device = "svg", width = 270, height = 240, units = "mm")
ggsave(filename="./Figures/mangroveVars_v2.png", plot = p.joint, device = "png", width = 270, height = 240, units = "mm")
ggsave(filename="./Figures/mangroveVars_v2.pdf", plot = p.joint, device = cairo_pdf, width = 270, height = 240, units = "mm")

cairo_ps(filename="./Figures/mangroveVars_v2.eps", width = 10.63, height = 9.45, 
         pointsize = 10, fallback_resolution = 2400)
print(p.joint)
dev.off()