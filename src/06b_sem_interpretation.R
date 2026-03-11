# ==============================================================================
# Script Name: Post-Modelling Visualisation and Analysis
# Description: Generates heatmaps, structural equation model (SEM) path diagrams, 
#              effect size bar charts, and spatial maps for the Mangrove-Malaria 
#              analysis. Evaluates model robustness across spatial resolutions.
# ==============================================================================

##----------------------------------------------------------##
## 1. Load Required Libraries
##----------------------------------------------------------##
# Data wrangling
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)

# Statistics and modelling
library(piecewiseSEM)
library(DHARMa)
library(nlme)
library(performance)
library(MASS)

# Parallel processing & progress bars
library(pbapply)
library(doParallel)

# Plotting & Visualisation
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(GGally)
library(ComplexHeatmap)
library(circlize)

# Spatial data processing
library(sf)
library(ggsflabel)

# Graph and SVG generation (DiagrammeR)
library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)

##----------------------------------------------------------##
## 2. Load Custom Functions
##----------------------------------------------------------##
# Load helper functions for model fitting and comparison
source("./src/sem_functions.R")

##----------------------------------------------------------##
## 3. Robustness Checks (Heatmap Generation)
##----------------------------------------------------------##

# Create a lookup table for publication-ready labels to use in plots
pretty.labels = data.frame(
  nodes = c("pr.new", "mean_ndvi", "mangrove_cover", "mangrove.cover.min1",
            "coastline.dist", "lc_crop", "pop_dens_median",
            "anomaly_2t", "anomaly_tp",
            "mean_2t_6m", "mean_tp_6m",
            "mean_2t", "mean_tp", 
            "anomaly_2t_6m", "anomaly_tp_6m"),
  new.labels.short = c(
    "MP", "MN", 
    "MC","MC(-1)",
    "CD", "AL", "PD",
    "T_Anom", "P_Anom", "T_6m", "P_6m", "T", "P", "T_Anom_6m", "P_Anom_6m"
  )
)

# Find and catalogue all model files in the data directory
file.list = list.files("./data/", pattern = "sem_results.+\\.rds")
file.table = as.data.frame(file.list) %>%
  mutate(fl2 = str_sub(file.list, start = 13, end = -5)) %>%
  separate(fl2, c("round", "type"), sep = "_") %>%
  mutate(type = ifelse(is.na(type), "", type))

# Load all model files, extract estimates, and reorganise them for comparison
opt.models = lapply(1:6, function(mtype){
  
  # Identify the latest optimisation round for each model type
  last.round = (filter(file.table, type == unique(file.table$type)[mtype]) %>%
                  slice_max(round))$file.list
  
  # Load the respective optimised model
  best.model = readRDS(paste0("./data/", last.round))
  
  # Extract standardised estimates across all 50 radii (buffers)
  pblapply(1:50, function(buffer) {
    if(class(best.model[[buffer]]) != "try-error") {
      best.model[[buffer]][[4]] %>%
        dplyr::select(Response, Predictor, P.Value, Std.Estimate) %>%
        mutate(buffer = buffer, mtype = unique(file.table$type)[mtype])
    }
  }) %>% bind_rows() %>%
    # Strip suffixes from 5km and 20km model predictors for clean joining
    mutate(
      Predictor = ifelse(endsWith(Predictor, "_20") | endsWith(Predictor,"_5"),
                         str_replace(Predictor, "\\_(5|20)", ""), Predictor)) %>%
    # Filter out quadratic weather terms and extraneous health variables
    filter(!grepl("sqr$", Predictor) &
             !(Predictor %in% c("t_healthcare_motor.median", 
                                "ITN_access_mean_median"))) %>%
    # Join with publication-ready labels
    left_join(pretty.labels , by = c("Predictor" = "nodes"))  %>%
    left_join(pretty.labels , by = c("Response" = "nodes")) %>%
    mutate(Predictor = new.labels.short.x,
           Response = new.labels.short.y) %>%
    dplyr::select(-starts_with("new"))
})

# Isolate the main (best) model
best.model = opt.models[[1]] %>% dplyr::select(-mtype)

# Create difference matrices comparing alternative models to the best model
robust.mat = lapply(2:6, function(eff){
  
  eff.size = opt.models[[eff]] %>%
    dplyr::select(-mtype)
  
  best.model %>%
    left_join(eff.size, by = c("Response", "Predictor", "buffer")) %>%
    mutate(Estimate.diff = Std.Estimate.x - Std.Estimate.y) %>%
    dplyr::select(-P.Value.x, -P.Value.y, -Std.Estimate.x, -Std.Estimate.y) %>%
    pivot_wider(names_from = buffer, values_from = Estimate.diff) %>%
    unite("relationship", Response:Predictor, sep = "-") %>%
    column_to_rownames("relationship")
})

# Build the robustness heatmap list
ht_list = NULL  ## Heatmap(...) + NULL yields a HeatmapList object
title.vec = c("20 km", "5 km", "reduced", "new variables", "weather-sqr")

for(s in c(3,2,1,4,5)) {
  ht_list = ht_list + Heatmap(as.matrix(robust.mat[[s]][,-1]), 
                              cluster_rows = F, cluster_columns = F,
                              col=colorRamp2(c(-1, 0, 1), c("red", "white", "blue")),
                              column_title = title.vec[s],
                              show_heatmap_legend = ifelse(s == 1, T, F),
                              row_names_side = "left",
                              column_names_gp = gpar(fontsize = 10),
                              column_names_rot = 45,
                              row_title = "relationships",
                              heatmap_legend_param = list(
                                legend_height = unit(4, "cm"),
                                title = "Deviation"
                              ))
}

# Print heatmap
ht_list

# Save figure with heatmaps
svg(file="./Figures/heatmap_robustness_v1.svg", width=11.7, height=8.3)
draw(ht_list)
dev.off()

# Export best model as RDS file for the interactive Shiny App
sem_results_3 = readRDS("./data/sem_results_3.rds")
saveRDS(sem_results_3, 
        file = "./src/Mangrove-Malaria_ShinyApp/sem_results_3.rds")


##----------------------------------------------------------------
## 4. Support of Optimisation Steps (Violin Plots)
##----------------------------------------------------------------

# Load all models for Fisher's C extraction
models = pblapply(1:nrow(file.table), function(x) {
  readRDS(paste0("./data/", file.table[x, "file.list"]))
})

# Extract test results (Fisher's C) for each optimisation round and radius
model_support = lapply(1:length(models),
                       function(m){
                         # check for errors
                         lapply(1:50, function(x) {
                           if(class(models[[m]][[x]]) != "try-error") {
                             models[[m]][[x]][[3]] %>%
                               mutate(buffer = x)
                           }
                         })  %>%
                           bind_rows() %>%
                           mutate(model = str_sub(file.table$file.list, end = 13)[m],
                                  round = file.table$round[m],
                                  type = file.table$type[m])
                       }) %>%
  bind_rows() %>%
  mutate(type = factor(type, levels = c("", "small", "5k", "20k", "newvars", "sqr")))


# Create violin plot of Fisher's C statistic across models
opt.plot = ggplot(data = model_support, aes(x = round, y = Fisher.C, fill = type)) +
  geom_violin(trim = F, position = "dodge") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2", name = "Model type",
                    labels = c("main", "reduced", "5-km", "20-km",   
                               "new variables",  "weather-qdr")) +
  facet_grid( ~ round, scale = "free") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text.x = element_blank()
  ) +
  xlab('Optimisation step') + ylab("Fisher's C statistic")
opt.plot

# Export violin plots in multiple formats
ps.options(family = "Arial")
ggsave(filename="./Figures/ModelOptim_v3.svg", device = "svg", 
       width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/ModelOptim_v3.png", device = "png", 
       width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/ModelOptim_v3.pdf", device = cairo_pdf, 
       width = 210, height = 105, units = "mm")
cairo_ps(filename="./Figures/ModelOptim_v3.eps", 
         width = 8.3, height = 4.15, # in inches
         pointsize = 10, fallback_resolution = 2400)
opt.plot
dev.off()

##-------------------------------------------------------
## 5. Plotting Path Diagrams (DiagrammeR)
##-------------------------------------------------------


# Plotting function: creates both standard plots and plots with enlarged text
#                    (suffix "_large")
#                    - input_model: models to be plotted
#                    - ID: identifier to be added to output file name
#                    - best_rs = spatial resolution r of models to be selected
plotting.fct = function(input_model, ID, best_rs) {
  # Filter for the best fitting models based on Shipley's d-test (P.Value > 0.05)
  f_test_best = model_interpretation(input_model)[[1]] %>%
    filter(P.Value > .05)
  print(f_test_best, n = nrow(f_test_best))
  best_model = input_model
  
  # Calculate robustness across radii of all paths 
  paths_robust = pblapply(f_test_best$buffer, function(x) {
    if(class(best_model[[x]]) != "try-error") {
      # Remove columns with significant levels and select effects with p > 0.05
      best_model[[x]][[4]][,1:8] %>%
        filter(P.Value <= 0.05) %>%
        # Merge and select response-predictor combinations without support
        unite('paths', c(Response, Predictor), sep = '-') %>%
        dplyr::select(paths) %>% unlist
    }
  })
  
  # Create table for support of paths, linking back to internal numeric node IDs
  path.support = data.frame(paths = unlist(paths_robust)) %>%
    count(paths) %>%
    separate(paths, c('Response', 'Predictor'), sep = "-") %>%
    left_join(best_model[[4]][[1]] %>% get_node_df() %>% 
                dplyr::select(id, nodes),
              by = c('Response' = 'nodes')) %>%
    left_join(best_model[[4]][[1]] %>% get_node_df() %>% 
                dplyr::select(id, nodes),
              by = c('Predictor' = 'nodes')) %>%
    rename(to = id.x, from = id.y, penwidth = n) %>%
    dplyr::select(-Response, -Predictor)
  
  # Create labels for temperature and precipitation cluster nodes
  clusters = data.frame(
    label = c("anomaly\n(survey period)", "mean\n(-6 months)",
              "mean\n(survey period)", "anomaly\n(-6 months)"),
    x = 11,
    y = c(8.5, 6.5, 4.5, 2.5),
    fontsize = 14,
    color = "white",
    fillcolor = "white",
    width = 1
  )
  
  # Define static position (x, y coordinates) and labels of variables
  attr.table = data.frame(
    y = c(10, 7, 4, 1, 1, 4, 7,  9,  8,  7,  6,  5, 4,  3,  2),
    x = c( 5, 5, 4, 6, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10),
    nodes = c("pr.new", "mean_ndvi", "mangrove_cover", "mangrove.cover.min1",
              "coastline.dist", "lc_crop", "pop_dens_median",
              "anomaly_2t", "anomaly_tp",
              "mean_2t_6m", "mean_tp_6m",
              "mean_2t", "mean_tp", 
              "anomaly_2t_6m", "anomaly_tp_6m"),
    new.labels = c(
      "Malaria\nprevalence", "Mangrove\nNDVI", 
      "Mangrove cover\n(current year)","Mangrove cover\n(previous year)",
      "Coastline\ndistance", "Agricultural\nland cover", "Population\ndensity",
      "T", "P", "T", "P", "T", "P", "T", "P"
    ),
    new.labels.short = c(
      "MP", "MN", 
      "MC","MC\n(-1)",
      "CD", "AL", "PD",
      "T", "P", "T", "P", "T", "P", "T", "P"
    ),
    shape = "circle")
  
  # Add health variables to the attribute table if ID is 'newvars'
  if(ID == "newvars") {
    attr.table = rbind(attr.table,
                       data.frame(
                         y = c(11, 11), x = c(4, 6),
                         nodes = c("t_healthcare_motor.median", "ITN_access_mean_median"),
                         new.labels = c("Motorised\ntravel", "ITN\naccess"),
                         new.labels.short = c("MT", "ITN"),
                         shape = "circle"
                       ))
  }
  
  # Function to generate the individual SEM SVG plots
  graph_plot = function(r, graph.type) {
    graph.plot = best_model[[r]][[1]] %>%
      
      # Strip suffixes from 5k and 20k models so they match the master attr.table
      mutate_node_attrs(nodes = gsub("_5|_20", "", nodes)) %>%
      
      join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
      mutate_node_attrs(label = get(ifelse(graph.type == "large",
                                           "new.labels",
                                           "new.labels.short")),
                        width = if(graph.type == 'large') {
                          ifelse(x == 10, 0.5, 1.5)
                        } else {ifelse(x == 10, 1, 1.8)},
                        fontsize = if(graph.type == "large") {14} else{25}
      ) %>%
      drop_node_attrs(new.labels) %>%
      join_edge_attrs(path.support) %>%
      mutate_edge_attrs(
        # Set arrow colours: red for negative associations, blue for positive
        color = ifelse(
          style == "solid",ifelse(as.numeric(label) < 0, "red", "blue"), "grey"
        ),
        arrowsize = 1,
        # Define thickness of paths based on robustness across spatial radii
        penwidth = penwidth/5,
        taillabel = label,
        label = NA,
        fontsize = if(graph.type == "large") {
          ifelse(style == "solid", 14, 0)
        } else{ifelse(style == "solid", 25, 0)}) %>%
      add_global_graph_attrs(attr = "splines",
                             value = "spline", 
                             attr_type = "graph") %>%
      add_nodes_from_table(clusters, label_col = label, set_type = "cluster")
    
    
    # Export SEM figure as SVG
    export_graph(graph.plot, paste0("./Figures/psem_", r, "_", ID,
                                    "_", graph.type, "_v3.svg"),
                 file_type = "svg", title = paste0("r = ", r, " km"),
                 width = 1200, height = 1000) 
  }
  
  # Apply plotting function to selected spatial radii
  pblapply(best_rs, function(x) graph_plot(x, "large"))
  pblapply(best_rs, function(x) graph_plot(x, "small"))
}

# Run the plotting wrapper for the main and alternative models
plotting.fct(sem_results_3, "", best_rs = c(3,28,40))

sem_results_4_small = readRDS("./data/sem_results_4_small.rds")
sem_results_5_5k = readRDS("./data/sem_results_5_5k.rds")
sem_results_6_20k = readRDS("./data/sem_results_6_20k.rds")

plotting.fct(sem_results_4_small, "small", best_rs = c(3,28,40))
plotting.fct(sem_results_5_5k, "5km", best_rs = c(3,28,40))
plotting.fct(sem_results_6_20k, "20km", best_rs = c(3,28,40))


##-------------------------------------------------------
## 6. Plotting Initially Hypothesised Model Structure
##-------------------------------------------------------
sem_results_1 = readRDS("./data/sem_results_1.rds")

# Define position and labels of variables for the unoptimised baseline model
attr.table = data.frame(
  y = c(10, 7, 4, 1, 1, 4, 7,  9,  8,  7,  6,  5, 4,  3,  2),
  x = c( 5, 5, 4, 6, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10),
  nodes = c("pr.new", "mean_ndvi", "mangrove_cover", "mangrove.cover.min1",
            "coastline.dist", "lc_crop", "pop_dens_median",
            "anomaly_2t", "anomaly_tp",
            "mean_2t_6m", "mean_tp_6m",
            "mean_2t", "mean_tp", 
            "anomaly_2t_6m", "anomaly_tp_6m"),
  new.labels = c(
    "Malaria\nprevalence", "Mangrove\nNDVI", 
    "Mangrove cover\n(current year)","Mangrove cover\n(previous year)",
    "Coastline\ndistance", "Agricultural\nland cover", "Population\ndensity",
    "T", "P", "T", "P", "T", "P", "T", "P"
  ),
  new.labels.short = c(
    "MP", "MN", 
    "MC","MC\n(-1)",
    "CD", "AL", "PD",
    "T", "P", "T", "P", "T", "P", "T", "P"
  ),
  shape = "circle")

clusters = data.frame(
  label = c("anomaly\n(survey period)", "mean\n(-6 months)",
            "mean\n(survey period)", "anomaly\n(-6 months)"),
  x = 11,
  y = c(8.5, 6.5, 4.5, 2.5),
  fontsize = 14,
  color = "white",
  fillcolor = "white",
  width = 1
)

graph.plot.init = sem_results_1[[16]][[1]] %>%
  join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
  mutate_node_attrs(label = new.labels,
                    width = ifelse(x == 10, 0.5, 1.5),
                    fontsize = 14) %>%
  drop_node_attrs(new.labels) %>%
  mutate_edge_attrs(
    color = "black",
    arrowsize = 1,
    taillabel = NA,
    label = NA,
    style = "solid") %>%
  add_global_graph_attrs(attr = "splines",
                         value = "spline", 
                         attr_type = "graph")  %>%
  add_nodes_from_table(clusters, label_col = label, set_type = "cluster")

export_graph(graph.plot.init, paste0("./Figures/initial_hypothesis_v2.svg"),
             file_type = "svg",
             width = 1200, height = 1000) 

##---------------------------------------------------------
## 7. Plotting Standardised Effect Sizes
##---------------------------------------------------------

# Extract and aggregate significant effect sizes by variable groups
eff.size = pblapply(f_test_best$buffer, function(x) {
  if(class(best_model[[x]]) != "try-error") {
    best_model[[x]][[4]][,1:8] %>%
      mutate(buffer = x)
  }
}) %>% bind_rows() %>%
  unite('path', c(Response, Predictor), sep = "-") %>%
  # Only include significant effect sizes (p <= 0.05)
  mutate(Estimate = ifelse(P.Value <= 0.05, abs(Std.Estimate), NA)) %>%
  dplyr::select(path, Estimate, buffer) %>%
  pivot_wider(names_from = path, values_from = Estimate) %>%
  mutate(
    # Aggregate effect sizes by broader variable categories
    weather = rowSums(dplyr::select(., contains("tp") | contains("2t")), na.rm = T),
    ndvi = rowSums(dplyr::select(., contains("-mean_ndvi")), na.rm = T),
    mangrove_landcover = rowSums(dplyr::select(., contains("-mangrove_cover")), na.rm = T),
    agriculture = rowSums(dplyr::select(., contains("crop")), na.rm = T),
    population = rowSums(dplyr::select(., contains("pop")), na.rm = T),
    coastline = rowSums(dplyr::select(., contains("coast")), na.rm = T)
  ) %>%
  # Reorganise data into long format for ggplotting
  dplyr::select(weather, ndvi, mangrove_landcover, agriculture, 
                population, coastline, buffer) %>%
  pivot_longer(!buffer, names_to = 'path', values_to = 'Estimate')

# Plot aggregated effect sizes as stacked bar chart
effect.plot = ggplot(eff.size, aes(x = buffer, y = Estimate, fill = path)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2", 
                    labels = c("Agricultural land cover",
                               "Coastline distance",
                               "Mangrove land cover",
                               "Mangrove NDVI",
                               "Population density",
                               "Weather variables (sum)")) +
  ylab('Sum of standardised estimates') +
  xlab('Spatial scale for calculating mangrove variables [km]') +
  labs(fill = "Causal\nrelationships")
effect.plot

# Export effect size plots
ps.options(family = "Arial")
ggsave(filename="./Figures/EffectSizes_v2.svg", device = "svg", width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/EffectSizes_v2.png", device = "png", width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/EffectSizes_v2.pdf", device = cairo_pdf, width = 210, height = 105, units = "mm")
cairo_ps(filename="./Figures/EffectSizes_v2.eps", width = 8.3, height = 4.15, # in inches
         pointsize = 10, fallback_resolution = 2400)
effect.plot
dev.off()

##---------------------------------------------------
## 8. Additional Info & Summary Statistics
##---------------------------------------------------

# How many unique geographical coordinates were analysed overall?
read.csv("./data/alldata.impute.csv") %>%
  #distinct(lat, lon) %>% 
  nrow()

# How many unique coordinates at 1 km vs 50 km radii?
read.csv("./data/alldata.impute.csv") %>%
  filter(buffer == 1) %>%
  nrow()
read.csv("./data/alldata.impute.csv") %>%
  filter(buffer == 50) %>%
  nrow()

# How many unique coordinates from original pre-filtered dataset (50 km off coastline)?
read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = T)  %>%
  distinct(lat, lon) %>% nrow()
read.csv("./data/coastal.PR.final.csv", sep = ",", header = T)  %>%
  distinct(lat, lon) %>% nrow()

test = read.csv("./data/alldata.impute.csv") %>%
  filter(buffer == 50) %>% 
  group_by(lat, lon) %>% filter(n()>1) %>% summarize(n=n())

##----------------------------------------------------------
## 9. Plot Number of Observations per Spatial Resolution (r)
##----------------------------------------------------------
radius_plot = ggplot(data = read.csv("./data/alldata.impute.csv") %>%
                       count(buffer),
                     aes(x = buffer, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  xlab(expression(paste("radius ", italic("r"), " [km]"))) +
  ylab("n(mangrove data)")

##--------------------------------------------------------------
## 10. Plot Map of Unique Geographical Locations (Malaria Data)
##--------------------------------------------------------------
unique_coord = read.csv("./data/alldata.impute.csv") %>%
  distinct(lat, lon)

# Import background land shapefile (Natural Earth)
ne.land = st_read("./res/ne_10m_land/ne_10m_land.shp")

coord_plot = ggplot() + 
  geom_sf(data = ne.land, fill = "lemonchiffon3") +
  geom_point(data = unique_coord, 
             aes(x = lon, y = lat),
             colour = "red",
             size = 1) + 
  coord_sf(xlim = c(-19,52), ylim = c(-35, 38)) + 
  theme_bw() +
  theme(panel.background = element_rect(fill='skyblue', colour='red'),
        legend.position = 'none',
        text=element_text(size=11)) +
  xlab("Longitude") + ylab("Latitude")

# Combine map and bar chart into a single composite panel
data_overview = plot_grid(coord_plot, radius_plot, labels = c('A', 'B'), 
                          label_size = 12, rel_widths =c(5,5))

# Export geographic data overview maps
ps.options(family = "Arial")
ggsave(filename="./Figures/datacount_v2.svg", device = "svg", 
       width = 200, height = 100, units = "mm")
ggsave(filename="./Figures/datacount_v2.png", device = "png", 
       width = 200, height = 100, units = "mm")
ggsave(filename="./Figures/datacount_v2.pdf", device = cairo_pdf, 
       width = 200, height = 100, units = "mm")
cairo_ps(filename="./Figures/datacount_v1.eps", width = 8, height = 4, # in inches
         pointsize = 10, fallback_resolution = 2400)
data_overview
dev.off()

##----------------------------------------------------------------------------
## 11. Plot and Test Multicollinearity of Weather Variables
##----------------------------------------------------------------------------
# Note: Ensure the `all.data` object is loaded in your environment prior to running
X = all.data %>% dplyr::select(anomaly_2t:mean_tp_6m)

ps.options(family = "Arial")
ggsave(filename="Autocorrelation_v1.svg", device = "svg", 
       width = 210, height = 210, units = "mm")
ggsave(filename="Autocorrelation_v1.png", device = "png", 
       width = 210, height = 210, units = "mm")
ggsave(filename="Autocorrelation_v1.pdf", device = cairo_pdf, 
       width = 210, height = 210, units = "mm")
cairo_ps(filename="Autocorrelation_v1.eps", 
         width = 8.3, height = 8.3, # in inches
         pointsize = 10, fallback_resolution = 2400)
ggpairs(X)
dev.off()