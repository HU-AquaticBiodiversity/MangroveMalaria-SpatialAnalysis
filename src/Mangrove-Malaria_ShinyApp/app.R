# ==============================================================================
# Script Name: Mangrove-Malaria SEM Interactive Dashboard (Shiny App)
# Description: An interactive web application to explore Structural Equation 
#              Models (SEMs) of the mangrove-malaria relationship across 
#              different spatial resolutions (radii from 1 to 50 km).
# ==============================================================================

library(shiny)

##-------------------------------------------
## 1. Install & Load Required Packages
##-------------------------------------------
# Automatically check for and install missing packages before loading
list.of.packages <- c("DiagrammeR", "dplyr", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(DiagrammeR)
library(dplyr)
library(tidyr)

##-------------------------------------------
## 2. Load Pre-computed Model Results
##-------------------------------------------
# Load the main optimised models (containing outputs for all 50 spatial radii)
sem_results_3 = readRDS("sem_results_3.rds")

##-------------------------------------------
## 3. Helper Function for Checking Models
##-------------------------------------------
# Evaluates model fit and extracts path significance to calculate robustness
model_interpretation = function(sem_model) {
  
  # List results of Shipley's d-test (model fit) and filter for well-supported models
  # A p-value > 0.05 means the model structure is NOT rejected by the data
  f_test = lapply(1:50, function(x) {
    if(class(sem_model[[x]]) != "try-error") {
      unlist(c(sem_model[[x]][[3]], buffer = x))
    }
  }) %>% bind_rows() %>%
    filter(P.Value > .05)
  
  # Find all paths that are NOT significant (p > 0.05) across models
  # These are tracked to see which relationships drop out at certain spatial scales
  paths_to_drop = lapply(f_test$buffer, function(x) {
    if(class(sem_model[[x]]) != "try-error") {
      sem_model[[x]][[4]][,1:8] %>%
        filter(P.Value > 0.05) %>%
        unite('paths', c(Response, Predictor), sep = '-') %>%
        dplyr::select(paths) %>% unlist
    }
  })
  
  # Find paths that the directed separation (d-test) suggests are missing
  # and should theoretically be added to improve the model fit
  paths_to_add = lapply(f_test$buffer, function(x) {
    if(class(sem_model[[x]]) != "try-error") {
      sem_model[[x]][[2]][,1:5] %>%
        filter(P.Value < 0.05) %>%
        dplyr::select(Independ.Claim) %>% unlist
    }
  })
  
  # Output a list containing the test statistics and ratios of missing/unsupported paths
  list(f_test, 
       # Print proportion of models that do NOT support a specific existing path
       table(unlist(paths_to_drop))/nrow(f_test),
       # Print proportion of models suggesting a new path should be added
       table(unlist(paths_to_add))/nrow(f_test))
}

##-------------------------------------------
## 4. Extract Robustness Metrics for Plotting
##-------------------------------------------
# Identify the best supported models from the dataset
f_test_best = model_interpretation(sem_results_3)[[1]] %>%
  filter(P.Value > .05)
best_model = sem_results_3

# Calculate robustness across all spatial radii for all paths
# (i.e., how many times is a specific causal link statistically significant?)
paths_robust = lapply(f_test_best$buffer, function(x) {
  if(class(best_model[[x]]) != "try-error") {
    best_model[[x]][[4]][,1:8] %>%
      filter(P.Value <= 0.05) %>%
      unite('paths', c(Response, Predictor), sep = '-') %>%
      dplyr::select(paths) %>% unlist
  }
})

# Create a lookup table linking the robustness count to the internal graph node IDs.
# The 'n' variable will later dictate the thickness (penwidth) of the arrows.
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

##-------------------------------------------
## 5. Define Graph Aesthetics (Nodes & Clusters)
##-------------------------------------------
# Create group labels and bounding boxes for the weather variables
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

# Create a master attribute table defining the static x/y layout and text labels
attr.table = data.frame(
  # x & y positions of variables on the canvas
  y = c(10, 7, 4, 1, 1, 4, 7,  9,  8,  7,  6,  5, 4,  3,  2),
  x = c( 5, 5, 4, 6, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10),
  # Exact names of variables as they appear in the model output
  nodes = c("pr.new", "mean_ndvi", "mangrove_cover", "mangrove.cover.min1",
            "coastline.dist", "lc_crop", "pop_dens_median",
            "anomaly_2t", "anomaly_tp",
            "mean_2t_6m", "mean_tp_6m",
            "mean_2t", "mean_tp", 
            "anomaly_2t_6m", "anomaly_tp_6m"),
  # Publication-ready labels for the plots
  new.labels = c(
    "Malaria\nprevalence", "Mangrove\nNDVI", 
    "Mangrove cover\n(current year)","Mangrove cover\n(previous year)",
    "Coastline\ndistance", "Agricultural\nland cover", "Population\ndensity",
    "T", "P", "T", "P", "T", "P", "T", "P"
  ),
  fontsize = 14,
  shape = "circle")

##------------------------------------------------------
## 6. Set up Shiny App UI (User Interface)
##------------------------------------------------------


ui <- fluidPage(
  titlePanel("Path diagrams of the mangrove-malaria relationship"),
  
  # User input slider for spatial resolution r (1-50 km). Default set to 22 km.
  sliderInput("r", 
              "Spatial resolution r in km at which mangrove land cover and mangrove NDVI are calculated", 
              value = 22, min = 1, max = 50),
  
  # Allocate space for the DiagrammeR output
  grVizOutput('diagram', width = "100%", height = "760px") 
)

##------------------------------------------------------
## 7. Set up Shiny App Server (Backend Logic)
##------------------------------------------------------
server <- function(input, output) {
  
  # Render the reactive diagram based on the slider input (input$r)
  output$diagram <- renderDiagrammeR(render_graph(
    
    # Select the model corresponding to the user's chosen radius
    best_model[[input$r]][[1]] %>%
      
      # Attach the static layout coordinates and labels
      join_node_attrs(df = attr.table, by_df = "nodes", by_graph ="nodes") %>%
      mutate_node_attrs(label = new.labels,
                        width = ifelse(x == 10, 0.5, 1.6)) %>%
      drop_node_attrs(new.labels) %>%
      
      # Attach the calculated robustness metrics to the edges
      join_edge_attrs(path.support) %>%
      
      # Format edge colours: red for negative associations, blue for positive, grey for unsupported
      mutate_edge_attrs(color = ifelse(
        style == "solid",ifelse(as.numeric(label) < 0, "red", "blue"), "grey"
      ),
      arrowsize = 1,
      # Arrow width scaled by robustness (number of models supporting the path / 5)
      penwidth = penwidth/5,
      # Move effect sizes to the beginning of the paths (tails)
      taillabel = label,
      # Remove default centre path labels to reduce visual clutter
      label = NA,
      fontsize = ifelse(style == "solid", 20, 0)) %>%
      
      # Format paths as splines to prevent arrows from passing straight through nodes
      add_global_graph_attrs(attr = "splines",
                             value = "spline", 
                             attr_type = "graph") %>%
      
      # Overlay the weather variable bounding boxes
      add_nodes_from_table(clusters, label_col = label, set_type = "cluster")
  ))
}


##------------------------------------------------------
## 8. Run App
##------------------------------------------------------
shinyApp(ui = ui, server = server)