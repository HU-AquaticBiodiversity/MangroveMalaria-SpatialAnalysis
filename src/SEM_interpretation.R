library(dplyr)
library(piecewiseSEM)
library(DHARMa)
library(stringr)
library(nlme)
library(performance)
library(MASS)
library(pbapply)
library(doParallel)
library(DiagrammeR)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(rsvg)
library(DiagrammeRsvg)
library(magrittr)
library(GGally)
library(ggsflabel)
library(sf)
library(cowplot)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)

# load functions for model comparison
source("./src/functions.model_fitting.R")

##----------------------------------------------------------##
## Robustness checks
##----------------------------------------------------------##

# create data frame including shortened labels
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

# find model files
file.list = list.files("./data/", pattern = "sem_results.+\\.rds")
file.table = as.data.frame(file.list) %>%
  mutate(fl2 = str_sub(file.list, start = 13, end = -5)) %>%
  separate(fl2, c("round", "type"), sep = "_") %>%
  mutate(type = ifelse(is.na(type), "", type))

# load model files and re-organise
opt.models = lapply(1:6, function(mtype){

  # which is the last round for each model type?
  last.round = (filter(file.table, type == unique(file.table$type)[mtype]) %>%
     slice_max(round))$file.list
  
  # load respective optimised model
  best.model = readRDS(paste0("./data/", last.round))
  
  # extract estimates
  pblapply(1:50, function(buffer) {
    if(class(best.model[[buffer]]) != "try-error") {
      best.model[[buffer]][[4]] %>%
        dplyr::select(Response, Predictor, P.Value, Std.Estimate) %>%
        mutate(buffer = buffer, mtype = unique(file.table$type)[mtype])
    }
  }) %>% bind_rows() %>%
    mutate(
      Predictor = ifelse(endsWith(Predictor, "_20") | endsWith(Predictor,"_5"),
                         str_replace(Predictor, "\\_(5|20)", ""), Predictor)) %>%
    filter(!grepl("sqr$", Predictor) &
             !(Predictor %in% c("t_healthcare_motor.median", 
                                "ITN_access_mean_median"))) %>%
    left_join(pretty.labels , by = c("Predictor" = "nodes"))  %>%
    left_join(pretty.labels , by = c("Response" = "nodes")) %>%
    mutate(Predictor = new.labels.short.x,
           Response = new.labels.short.y) %>%
    dplyr::select(-starts_with("new"))
})

# define best model
best.model = opt.models[[1]] %>% dplyr::select(-mtype)

# create matrices
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

# build heatmaps
ht_list = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
title.vec = c("20 km", "5 km", "reduced", "new variables", "weather-sqr")
for(s in c(3,2,1,4,5)) {
  ht_list = ht_list + Heatmap(as.matrix(robust.mat[[s]][,-1]), 
                              cluster_rows = F, cluster_columns = F,
                              col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                              column_title = title.vec[s],
                              show_heatmap_legend = ifelse(s == 1, T, F),
                              row_names_side = "left",
                              column_names_gp = gpar(fontsize = 10),
                              column_names_rot = 45,
                              row_title = "relationships",
                              #column_title = "test",
                              heatmap_legend_param = list(
                                legend_height = unit(4, "cm"),
                                title = "Deviation"
                                ))
}

ht_list

# save figure with heatmaps
svg(file="./Figures/heatmap_robustness_v1.svg", width=11.7, height=8.3)
draw(ht_list)
dev.off()

# export best model as rds file for Shiny App
sem_results_3 = readRDS("./data/sem_results_3.rds")
saveRDS(sem_results_3, 
        file = "./src/Mangrove-Malaria_ShinyApp/sem_results_3.rds")


##----------------------------------------------------------------
## Support of optimisation steps
##----------------------------------------------------------------

models = pblapply(1:nrow(file.table), function(x) {
  readRDS(paste0("./data/", file.table[x, "file.list"]))
})

rm(m)

# extract test results
model_support = lapply(1:length(models),
                       function(m){
                         # check for errors
                         lapply(1:50, function(x) {
                           if(class(models[[m]][[x]]) != "try-error") {
                             models[[m]][[x]][[3]]
                           }
                         })  %>%
                           bind_rows() %>%
                           mutate(model = str_sub(file.table$file.list, end = 13)[m],
                                  round = file.table$round[m],
                                  type = file.table$type[m])
                       }) %>%
  bind_rows() %>%
  mutate(type = factor(type, levels = c("", "small", "5k", "20k", "newvars", "sqr")))


# create plot with test results as violin plots
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

# export plots
ps.options(family = "Arial")
ggsave(filename="./Figures/ModelOptim_v3.svg", device = "svg", 
       width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/ModelOptim_v3.png", device = "png", 
       width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/ModelOptim_v3.pdf", device = cairo_pdf, 
       width = 210, height = 105, units = "mm")
cairo_ps(filename="./Figures/ModelOptim_v3.eps", 
         width = 8.3, height = 4.15, #in inches
         pointsize = 10, fallback_resolution = 2400)
opt.plot
dev.off()

##-------------------------------------------------------
## Plotting path diagrams
##------------------------------------------------------

# plotting function: creates both standard plots and plots with enlarged text
#                    (suffix "_large')
#                    - input_model: models to be plotted
#                    - identifier to be added to output file name
#                    - best_rs = spatial resolution r of models to be selected
plotting.fct = function(input_model, ID, best_rs) {
  # --> best model is in sem_results_3
  f_test_best = model_interpretation(input_model)[[1]] %>%
    filter(P.Value > .05)
  print(f_test_best, n = nrow(f_test_best))
  best_model = input_model

  # calculate robustness across radii of all paths #
  paths_robust = pblapply(f_test_best$buffer, function(x) {
    # check for errors
    if(class(best_model[[x]]) != "try-error") {
      # remove columns with significant levels
      best_model[[x]][[4]][,1:8] %>%
        # select effects with p > 0.05
        filter(P.Value <= 0.05) %>%
        # merge and select response-predictor combinations without support
        unite('paths', c(Response, Predictor), sep = '-') %>%
        dplyr::select(paths) %>% unlist
    }
  })

  # create table for support of paths #
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
  
  # create labels for temperature and precipitation variables
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

  # define position and labels of variables
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

  # plotting SEMs
  graph_plot = function(r, graph.type) {
    graph.plot = best_model[[r]][[1]] %>%
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
        color = ifelse(
          style == "solid",ifelse(as.numeric(label) < 0, "red", "blue"), "grey"
          ),
        arrowsize = 1,
        # define width of paths
        penwidth = penwidth/5,
        # position of path labels
        taillabel = label,
        label = NA,
        fontsize = if(graph.type == "large") {
          ifelse(style == "solid", 14, 0)
        } else{ifelse(style == "solid", 25, 0)}) %>%
      add_global_graph_attrs(attr = "splines",
                             value = "spline", 
                             attr_type = "graph") %>%
      add_nodes_from_table(clusters, label_col = label, set_type = "cluster")

    
    # export figure
    export_graph(graph.plot, paste0("./Figures/psem_", r, "_", ID,
                                    "_", graph.type, "_v3.svg"),
                 file_type = "svg", title = paste0("r = ", r, " km"),
                 width = 1200, height = 1000) 
  }
  
  pblapply(best_rs, function(x) graph_plot(x, "large"))
  pblapply(best_rs, function(x) graph_plot(x, "small"))
}

plotting.fct(sem_results_3, "", best_rs = c(3,28,40))

sem_results_1_newvars = readRDS("./data/sem_results_1_newvars.rds")

plotting.fct(sem_results_1_newvars, "newvars", best_rs = c(3,28,40))
sem_results_1_newvars[[1]]

plotting.fct(sem_results_2_small, "small_test", best_rs = c(21,36,2))

# plotting initially hypothesised model structure
graph.plot.init = sem_results_1[[16]][[1]] %>%
  join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
  mutate_node_attrs(label = new.labels,
                    width = ifelse(x == 10, 0.5, 1.5)) %>%
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

export_graph(graph.plot.init, paste0("initial_hypothesis_v1.png"),
             file_type = "png",
             width = 1200, height = 1000) 

##---------------------------------------------------------
## Plotting effect sizes
##---------------------------------------------------------

eff.size = pblapply(f_test_best$buffer, function(x) {
  if(class(best_model[[x]]) != "try-error") {
    best_model[[x]][[4]][,1:8] %>%
      mutate(buffer = x)
  }
}) %>% bind_rows() %>%
  unite('path', c(Response, Predictor), sep = "-") %>%
  # only include significant effect sizes
  mutate(Estimate = ifelse(P.Value <= 0.05, abs(Std.Estimate), NA)) %>%
  dplyr::select(path, Estimate, buffer) %>%
  pivot_wider(names_from = path, values_from = Estimate) %>%
  mutate(
    # aggregate effect sizes by variables groups
    weather = rowSums(dplyr::select(., contains("tp") | contains("2t")), na.rm = T),
    ndvi = rowSums(dplyr::select(., contains("-mean_ndvi")), na.rm = T),
    mangrove_landcover = rowSums(dplyr::select(., contains("-mangrove_cover")), na.rm = T),
    agriculture = rowSums(dplyr::select(., contains("crop")), na.rm = T),
    population = rowSums(dplyr::select(., contains("pop")), na.rm = T),
    coastline = rowSums(dplyr::select(., contains("coast")), na.rm = T)
  ) %>%
  # re-organise data for plotting
  dplyr::select(weather, ndvi, mangrove_landcover, agriculture, 
                population, coastline, buffer) %>%
  pivot_longer(!buffer, names_to = 'path', values_to = 'Estimate')

# plotting effect sizes
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

# export plot
ps.options(family = "Arial")
ggsave(filename="./Figures/EffectSizes_v2.svg", device = "svg", width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/EffectSizes_v2.png", device = "png", width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/EffectSizes_v2.pdf", device = cairo_pdf, width = 210, height = 105, units = "mm")
cairo_ps(filename="./Figures/EffectSizes_v2.eps", width = 8.3, height = 4.15, #in inches
         pointsize = 10, fallback_resolution = 2400)
effect.plot
dev.off()

##---------------------------------------------------
# Additional info
#---------------------------------------------------
# How many unique coordinates? That were analysed
read.csv("./data/alldata.impute.csv") %>%
  #distinct(lat, lon) %>% 
  nrow()

# How many unique coordinates? At 1 and 50 km
read.csv("./data/alldata.impute.csv") %>%
  filter(buffer == 1) %>%
  nrow()
read.csv("./data/alldata.impute.csv") %>%
  filter(buffer == 50) %>%
  nrow()

# How many unique coordinates from original dataset (50 km off coastline)
read.csv("./data/coastal.PR.dhs.final.csv", sep = ",", header = T)  %>%
  distinct(lat, lon) %>% nrow()
read.csv("./data/coastal.PR.final.csv", sep = ",", header = T)  %>%
  distinct(lat, lon) %>% nrow()

test = read.csv("./data/alldata.impute.csv") %>%
  filter(buffer == 50) %>% 
  group_by(lat, lon) %>% filter(n()>1) %>% summarize(n=n())

##----------------------------------------------------------
## plot number of observation per spatial resolution r
##----------------------------------------------------------
radius_plot = ggplot(data = read.csv("./data/alldata.impute.csv") %>%
                       count(buffer),
                     aes(x = buffer, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  xlab("radius r [km]") + ylab("n(mangrove data)")

##--------------------------------------------------------------
## plot map of all unique geographic locations of malaria data
##--------------------------------------------------------------
unique_coord = all.data %>%
  distinct(lat, lon)

# import background map
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

ps.options(family = "Arial")
ggsave(filename="datacount_v2.svg", device = "svg", 
       width = 200, height = 100, units = "mm")
ggsave(filename="datacount_v2.png", device = "png", 
       width = 200, height = 100, units = "mm")
ggsave(filename="datacount_v2.pdf", device = cairo_pdf, 
       width = 200, height = 100, units = "mm")
cairo_ps(filename="datacount_v1.eps", width = 8, height = 4, #in inches
         pointsize = 10, fallback_resolution = 2400)
# put two plots in single panel
data_overview= plot_grid(coord_plot, radius_plot, labels = c('A', 'B'), 
                         label_size = 12, rel_widths =c(5,5))
data_overview
dev.off()

##----------------------------------------------------------------------------
## plot and test multicollinearity of weather variables
##----------------------------------------------------------------------------
X = all.data %>% dplyr::select(anomaly_2t:mean_tp_6m)

ps.options(family = "Arial")
ggsave(filename="Autocorrelation_v1.svg", device = "svg", 
       width = 210, height = 210, units = "mm")
ggsave(filename="Autocorrelation_v1.png", device = "png", 
       width = 210, height = 210, units = "mm")
ggsave(filename="Autocorrelation_v1.pdf", device = cairo_pdf, 
       width = 210, height = 210, units = "mm")
cairo_ps(filename="Autocorrelation_v1.eps", 
         width = 8.3, height = 8.3, #in inches
         pointsize = 10, fallback_resolution = 2400)
ggpairs(X)
dev.off()