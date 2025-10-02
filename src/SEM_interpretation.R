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

##-------------------------------------------
## Functions for checking models
##-------------------------------------------
model_interpretation = function(sem_model) {

  # list results of f test and filter out significant ones
  f_test = pblapply(1:50, function(x) {
    if(class(sem_model[[x]]) != "try-error") {
      unlist(c(sem_model[[x]][[3]], buffer = x))
    }
  }) %>% bind_rows() %>%
    filter(P.Value > .05)

  ## find all paths that need to be removed ##
  ## --> values indicate the share of models that do NOT support the relationship
  ##     1.00 means no model supports it
  paths_to_drop = pblapply(f_test$buffer, function(x) {
    # check for errors
    if(class(sem_model[[x]]) != "try-error") {
      # remove columns with significant levels
      sem_model[[x]][[4]][,1:8] %>%
        # select effects with p > 0.05
        filter(P.Value > 0.05) %>%
        # merge and select response-predictor combinations without support
        unite('paths', c(Response, Predictor), sep = '-') %>%
        dplyr::select(paths) %>% unlist
    }
  })
  
  ## find all paths that need to be added ##
  paths_to_add = pblapply(f_test$buffer, function(x) {
    # check for errors
    if(class(sem_model[[x]]) != "try-error") {
      # remove columns with signficant levels
      sem_model[[x]][[2]][,1:5] %>%
        # select effects with p > 0.05
        filter(P.Value < 0.05) %>%
        dplyr::select(Independ.Claim) %>% unlist
    }
  })
  
  #output results
  output = list(f_test, 
                # print list of all paths unsupported across support models
                table(unlist(paths_to_drop))/nrow(f_test),
                # print list of paths to be added
                table(unlist(paths_to_add))/nrow(f_test))
  names(output) = c("F_Test", "paths_to_drop", "paths_to_add")
  output
}

model_optim_info = function(sem_results, name) {
  ptd = data.frame(model_interpretation(sem_results)$paths_to_drop) %>%
    separate(col = Var1, c("dependent", "independent"), sep = "-") %>%
    filter(Freq == 1)
  
  write.table(ptd, file = paste0("./data/", name,"_ptd.csv"), sep = ",")
  
  pta = data.frame(model_interpretation(sem_results)$paths_to_add) %>%
    separate(col = Var1, c("dependent", "x", "independent", "y", "z"), sep = " ") %>%
    dplyr::select(-x, -y, -z) %>%
    filter(Freq >= .1)
  
  write.table(pta, paste0(file = "./data/", name,"_pta.csv"), sep = ",")
}

##----------------------------------------------------------##
## Model comparison
##----------------------------------------------------------##

# load results
load("./data/sem_results_1.Rdata")
load("./data/sem_results_small_1.Rdata")
load("./data/sem_results_1_sqr.Rdata")
load("./data/sem_results_2.Rdata")
load("./data/sem_results_small_2.Rdata")
load("./data/sem_results_2_sqr.Rdata")
load("./data/sem_results_3.Rdata")
load("./data/sem_results_small_3.Rdata")
load("./data/sem_results_3_sqr.Rdata")
load("./data/sem_results_4_sqr.Rdata")
load("./data/sem_results_5_sqr.Rdata")
load("./data/sem_results_6_sqr.Rdata")

# export best model as rds file for Shiny App
saveRDS(sem_results_3, 
        file = "./src/Mangrove-Malaria_ShinyApp/sem_results_3.rds")

# show model diagnostics
model_interpretation(sem_results_1)
model_interpretation(sem_results_small_1)
model_interpretation(sem_results_1_sqr)

model_interpretation(sem_results_2)
model_interpretation(sem_results_small_2)
model_interpretation(sem_results_2_sqr)

model_interpretation(sem_results_3)
model_interpretation(sem_results_small_3)
model_interpretation(sem_results_3_sqr)

model_interpretation(sem_results_4_sqr)

model_interpretation(sem_results_5_sqr)

model_interpretation(sem_results_6_sqr)

# export paths to drop and add to optimise models
model_optim_info(sem_results_1, "sem_results_1")
model_optim_info(sem_results_small_1, "sem_results_small_1")
model_optim_info(sem_results_1_sqr, "sem_results_1_sqr")

model_optim_info(sem_results_2, "sem_results_2")
model_optim_info(sem_results_small_2, "sem_results_small_2")
model_optim_info(sem_results_2_sqr, "sem_results_2_sqr")

model_optim_info(sem_results_3, "sem_results_3")
model_optim_info(sem_results_small_3, "sem_results_small_3")
model_optim_info(sem_results_3_sqr, "sem_results_3_sqr")

model_optim_info(sem_results_4_sqr, "sem_results_4_sqr")

model_optim_info(sem_results_5_sqr, "sem_results_5_sqr")

model_optim_info(sem_results_6_sqr, "sem_results_6_sqr")

# show formulas of optimised models
load("./data/start.formulas.3.Rdata")
load("./data/start.formulas.3.small.Rdata")
load("./data/start.formulas.6.sqr.Rdata")
start.formulas.3
start.formulas.3.small
start.formulas.6.sqr

##-------------------------------------------------------
## Plotting path diagrams
##------------------------------------------------------

# plotting function: creates both standard plots and plots with enlarged text
#                    (suffix "_large')
#                    - input_model: models to be plotted
#                    - identifier to be added to output file name
#                    - best_rs = spatial resolution r of models to be selected
plotting.fct = function(input_model, ID, best_rs){
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
    fontsize = 14,
    shape = "circle")
  
  # plotting SEMs
  pblapply(best_rs, function(r) {
    graph.plot = best_model[[r]][[1]] %>%
      join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
      mutate_node_attrs(label = new.labels,
                        width = ifelse(x == 10, 0.5, 1.5)) %>%
      drop_node_attrs(new.labels) %>%
      join_edge_attrs(path.support) %>%
      mutate_edge_attrs(color = ifelse(
        style == "solid",ifelse(as.numeric(label) < 0, "red", "blue"), "grey"
      ),
      arrowsize = 1,
      # define width of paths
      penwidth = penwidth/5,
      # position of path labels
      taillabel = label,
      label = NA,
      fontsize = ifelse(style == "solid", 14, 0)) %>%
      add_global_graph_attrs(attr = "splines",
                             value = "spline", 
                             attr_type = "graph") %>%
      add_nodes_from_table(clusters, label_col = label, set_type = "cluster")
    
    # export figure
    export_graph(graph.plot, paste0("./Figures/psem_", r, "_", ID, "_v2.svg"),
                 file_type = "svg", title = paste0("r = ", r, " km"),
                 width = 1200, height = 1000) 
  })
  
  # make larger version of paths above
  pblapply(best_rs, function(r) {
    graph.plot = best_model[[r]][[1]] %>%
      join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
      mutate_node_attrs(label = new.labels,
                        width = ifelse(x == 10, 1, 1.8)) %>%
      drop_node_attrs(new.labels) %>%
      join_edge_attrs(path.support) %>%
      mutate_edge_attrs(color = ifelse(
        style == "solid",ifelse(as.numeric(label) < 0, "red", "blue"), "grey"
      ),
      arrowsize = 1,
      penwidth = penwidth/5,
      taillabel = label,
      label = NA,
      fontsize = ifelse(style == "solid", 20, 0)) %>%
      add_global_graph_attrs(attr = "splines",
                             value = "spline", 
                             attr_type = "graph") %>%
      add_nodes_from_table(clusters, label_col = label, set_type = "cluster")
    
    export_graph(graph.plot, paste0("./Figures/psem_large_", r, "_", ID,"_v2.svg"),
                 file_type = "svg", title = paste0("r = ", r, " km"),
                 width = 1200, height = 1000) 
  })
}

plotting.fct(sem_results_, "", best_rs = c(3,28,40))
sem_results_small_3[[40]]
plotting.fct(sem_results_small_3, "small", best_rs = c(3,28,40))

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

##----------------------------------------------------------------
## plot model support of optimisation steps
##----------------------------------------------------------------
# create two vectors/list to extract results of Fisher's exact test across
#  all models
mod.vcts = data.frame(mod = c('m1','m2', 'm3', 'm1','m2','m3','m4', 'm5', 'm6'),
                      sqr = c('A','A','A','B','B', 'B','B','B','B')
)
models = list(sem_results_1, sem_results_2, sem_results_3,
              sem_results_1_sqr, sem_results_2_sqr, sem_results_3_sqr,
              sem_results_4_sqr, sem_results_5_sqr, sem_results_6_sqr)



# extract test results
model_support = lapply(1:nrow(mod.vcts),
                       function(m){
                         # check for errors
                         lapply(1:50, function(x) {
                           if(class(models[[m]][[x]]) != "try-error") {
                             models[[m]][[x]][[3]]
                           }
                         })  %>%
                           bind_rows() %>%
                           mutate(model = mod.vcts$mod[m],
                                  sqr = mod.vcts$sqr[m])
                       }) %>%
  bind_rows()

# create plot with test results as violin plots
opt.plot = ggplot(data = model_support, aes(x = model, y = Fisher.C, fill = model)) +
  geom_violin(trim = F) +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(~ sqr, scale = "free") +
  theme(legend.position = "none",
        strip.text.x = element_text(
          size = 12, face = "bold",
          hjust = 0
        )) +
  xlab('Model') + ylab("Fisher's C statistic")
opt.plot

# export plots
ps.options(family = "Arial")
ggsave(filename="./Figures/ModelOptim_v2.svg", device = "svg", 
       width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/ModelOptim_v2.png", device = "png", 
       width = 210, height = 105, units = "mm")
ggsave(filename="./Figures/ModelOptim_v2.pdf", device = cairo_pdf, 
       width = 210, height = 105, units = "mm")
cairo_ps(filename="./Figures/ModelOptim_v2.eps", 
         width = 8.3, height = 4.15, #in inches
         pointsize = 10, fallback_resolution = 2400)
opt.plot
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