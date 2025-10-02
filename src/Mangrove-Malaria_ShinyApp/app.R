library(shiny)

## install missing packages ##
list.of.packages <- c("DiagrammeR", "dplyr", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(DiagrammeR)
library(dplyr)
library(tidyr)

##-------------------------------------------
## LOAD results of models
##------------------------------------------
sem_results_3 = readRDS("sem_results_3.rds")

##-------------------------------------------
## Function for checking models
##-------------------------------------------
model_interpretation = function(sem_model) {
  # list results of f test and filter out significant ones
  f_test = lapply(1:50, function(x) {
    if(class(sem_model[[x]]) != "try-error") {
      unlist(c(sem_model[[x]][[3]], buffer = x))
    }
  }) %>% bind_rows() %>%
    filter(P.Value > .05)
  
  ## find all paths that need to be removed ##
  ## --> values indicate the share of models that do NOT support the relationship
  ##     1.00 means no model supports it
  paths_to_drop = lapply(f_test$buffer, function(x) {
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
  paths_to_add = lapply(f_test$buffer, function(x) {
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
  list(f_test, 
       # print list of all paths unsupported across support models
       table(unlist(paths_to_drop))/nrow(f_test),
       # print list of paths to be added
       table(unlist(paths_to_add))/nrow(f_test))
}

## --> best model is in sem_results_3
f_test_best = model_interpretation(sem_results_3)[[1]] %>%
  filter(P.Value > .05)
best_model = sem_results_3

## calculate robustness across radii of all paths ##
paths_robust = lapply(f_test_best$buffer, function(x) {
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

#create group labels for weather variables
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

# create list of attributed for path plot in DiagrammeR
attr.table = data.frame(
  # x & y positions of variables
  y = c(10, 7, 4, 1, 1, 4, 7,  9,  8,  7,  6,  5, 4,  3,  2),
  x = c( 5, 5, 4, 6, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10),
  # names of variables in dataset
  nodes = c("pr.new", "mean_ndvi", "mangrove_cover", "mangrove.cover.min1",
            "coastline.dist", "lc_crop", "pop_dens_median",
            "anomaly_2t", "anomaly_tp",
            "mean_2t_6m", "mean_tp_6m",
            "mean_2t", "mean_tp", 
            "anomaly_2t_6m", "anomaly_tp_6m"),
  # names of variables on plots
  new.labels = c(
    "Malaria\nprevalence", "Mangrove\nNDVI", 
    "Mangrove cover\n(current year)","Mangrove cover\n(previous year)",
    "Coastline\ndistance", "Agricultural\nland cover", "Population\ndensity",
    "T", "P", "T", "P", "T", "P", "T", "P"
  ),
  fontsize = 14,
  shape = "circle")

##------------------------------------------------------
## Set up Shiny app
##------------------------------------------------------

ui <- fluidPage(
  titlePanel("Path diagrams of the mangrove-malaria relationship"),
  # set up input 1-50 km (default at 22 km)
  sliderInput("r", "Spatial resolution r in km at which mangrove land cover and mangrove NDVI are calculated", value = 22, min = 1, max = 50),
  # size of plot
  grVizOutput('diagram', width = "100%", height = "760px") 
)

server <- function(input, output) {
  # render diagrame
  output$diagram <- renderDiagrammeR(render_graph(
    best_model[[input$r]][[1]] %>%
      join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
      mutate_node_attrs(label = new.labels,
                        width = ifelse(x == 10, 0.5, 1.6)) %>%
      drop_node_attrs(new.labels) %>%
      join_edge_attrs(path.support) %>%
      # path colour by negative (blue), positive (red), or unsupported (grey)
      mutate_edge_attrs(color = ifelse(
        style == "solid",ifelse(as.numeric(label) < 0, "red", "blue"), "grey"
      ),
      arrowsize = 1,
      # width of paths = number of models supporting the path, divided by 5
      penwidth = penwidth/5,
      # effect sizes at beginning of paths
      taillabel = label,
      # turn of path labels at centre of paths
      label = NA,
      fontsize = ifelse(style == "solid", 20, 0)) %>%
      # make paths as splines to avoid crossing nodes
      add_global_graph_attrs(attr = "splines",
                             value = "spline", 
                             attr_type = "graph") %>%
      add_nodes_from_table(clusters, label_col = label, set_type = "cluster")
  ))
}

##----------------
## RUN APP
##----------------
shinyApp(ui = ui, server = server)
