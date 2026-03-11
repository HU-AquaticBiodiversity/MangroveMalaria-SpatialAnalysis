# ==============================================================================
# ==============================================================================
# SPATIAL SEM BOOTSTRAPPING & BACK-TRANSFORMATION SCRIPT
#
# Description: 
# PART 1: Generates a cluster bootstrap for a piecewise SEM across spatial scales,
# accounting for spatial autocorrelation, and saves the raw iterations.
# PART 2: Filters non-converged runs, calculates CIs, and back-transforms direct, 
# indirect, and total effects into comparative units and Odds Ratios.
# ==============================================================================
# ==============================================================================

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
library(piecewiseSEM) 
library(dplyr)        
library(tidyverse)    
library(doParallel)   
library(MASS)         
library(nlme)         
library(performance)  
library(parallel)     

# ==============================================================================
# PART 1: BOOTSTRAPPING PIPELINE (DATA GENERATION)
# ==============================================================================

# --- 1A. Global Setup ---
source("./src/functions.model_fitting.R")

# Detect and set the number of CPU cores
numberOfCores <- 36

# Load the previously identified 'best' model formulas
f <- readRDS("./data/start.formulas_4.rds")

# --- 1B. Main Bootstrap Wrapper Function ---
boot.fct <- function(r) {
  
  cat("\n==================================================\n")
  cat("STARTING RADIUS:", r, "km\n")
  cat("==================================================\n")

  # 1. Data Loading & Preparation
  input_set <- read.csv("./data/alldata.impute.csv") %>% 
    data.prep() %>% 
    filter(buffer == r)
  
  # 2. Initial Master Models
  # MODEL A: Prevalence (Binomial glmmPQL with spatial autocorrelation)
  model.pr.1 <- glmmPQL(
    fixed = as.formula(paste(deparse(f[[1]]), collapse ="")),
    random = ~ 1 | group,
    correlation = corExp(1, form = ~ lat + lon | group, nugget = TRUE),
    data = input_set,
    weights = input_set$examined, 
    family = binomial(link = "logit"),
    control = lmeControl(msMaxIter = 1000, msMaxEval = 1000), 
    verbose = FALSE
  )
  
  # MODEL B: Mangrove NDVI (Gaussian)
  model.vi.1 <- update(
    model.pr.1,
    fixed = as.formula(paste(deparse(f[[2]]), collapse ="")),
    family = gaussian
  )
  
  # MODEL C: Mangrove Cover (Gaussian)
  model.mc.1 <- update(
    model.vi.1,
    fixed = as.formula(paste(deparse(f[[3]]), collapse =""))
  )
  
  # 3. Parallel Cluster Configuration
  n_boots = 1000
  group_var = "group"
  
  cl <- makeCluster(numberOfCores, outfile = "") 
  clusterSetRNGStream(cl, iseed = 123)
  clusterExport(cl, varlist = c("n_boots", "group_var", "model.pr.1", 
                                "model.mc.1", "model.vi.1", "input_set", "r"),
                envir = environment())
  
  # 4. Parallel Bootstrapping Loop
  system.time({
    all_boots <- parLapply(cl, 1:n_boots, function(i) {
      
      suppressPackageStartupMessages({
        library(piecewiseSEM)
        library(dplyr)
        library(nlme)
        library(MASS)
      })
      
      # --- Step A: Resample Groups ---
      groups <- unique(input_set[[group_var]])
      sampled_groups <- sample(groups, size = length(groups), replace = TRUE)
      
      boot_data <- do.call(rbind, lapply(sampled_groups, function(g) {
        input_set[input_set[[group_var]] == g, ]
      }))
      
      # --- Step B: Jitter Coordinates ---
      boot_data$lat <- jitter(boot_data$lat, amount = 0.00001)
      boot_data$lon <- jitter(boot_data$lon, amount = 0.00001)
      
      # --- Step C: Environment Assignment ---
      assign("boot_data", boot_data, envir = .GlobalEnv)
      
      # --- Step D: Update Models ---
      mod_A_new <- try(update(model.pr.1, data = boot_data, weights = boot_data$examined), silent = TRUE)
      mod_B_new <- try(update(model.vi.1, data = boot_data, weights = boot_data$examined), silent = TRUE)
      mod_C_new <- try(update(model.mc.1, data = boot_data, weights = boot_data$examined), silent = TRUE)
      
      if (inherits(mod_A_new, "try-error") || 
          inherits(mod_B_new, "try-error") || 
          inherits(mod_C_new, "try-error")) {
        return(NULL)
      }
      
      # --- Step E: Extract SEM Coefficients ---
      sem_new <- suppressWarnings(psem(mod_A_new, mod_B_new, mod_C_new, data = boot_data))
      stats <- try(coefs(sem_new), silent = TRUE)
      
      if (inherits(stats, "try-error")) return(NULL)
      
      res <- stats[, c("Response", "Predictor", "Estimate")]
      res$Iteration <- i
      
      if (i %% 50 == 0) cat("Radius", r, "km - Completed Iteration:", i, "\n")
      
      return(res)
    }) %>% bind_rows() 
  })
  
  stopCluster(cl)
  
  # 5. Save Raw Results 
  raw_out_file <- paste0("./data/raw_bootstraps_", r, "km.rds")
  saveRDS(all_boots, file = raw_out_file)
  cat("Successfully saved RAW iterations to:", raw_out_file, "\n")
}

# --- 1C. Execute Bootstrap ---
lapply(c(3, 28, 40), function(r) boot.fct(r))


# ==============================================================================
# PART 2: BACK-TRANSFORMATION & COMPOUND PATHS (DATA PROCESSING)
# ==============================================================================

# --- 2A. Setup & Prep for Back-transformation ---
sem_results_3 <- readRDS("./data/sem_results_3.rds") 
all.data.raw <- read.csv("./data/alldata.full.csv") %>% data.prep()
path_map <- read.csv("./data/all_paths_input.csv", stringsAsFactors = FALSE)

all_vars <- unique(c(sem_results_3[[40]][[4]]$Predictor, sem_results_3[[40]][[4]]$Response))
sd.table <- lapply(all_vars, function(x) {
  value <- sd(unlist(dplyr::select(all.data.raw, all_of(x))), na.rm = TRUE)
  data.frame(var.name = x, sd.value = value)
}) %>% bind_rows()

comp.unit.table <- sd.table %>%
  mutate(comp.unit = ifelse(str_detect(var.name, "2t") | var.name == "coastline.dist", 1,
                            ifelse(str_detect(var.name, "tp"), 0.01,
                                   ifelse(var.name == "pop_dens_median", 100, 0.1)))) %>%
  dplyr::select(var.name, comp.unit)

calc_path_per_iter <- function(path_name, n1, n2, n3, n4, boots) {
  if (is.na(n3) | n3 == "NA" | n3 == "") {
    boots %>% filter(Predictor == n1 & Response == n2) %>%
      dplyr::select(Iteration, Path_Est = Estimate) %>%
      mutate(Predictor = n1, Final_Response = n2, Path_Name = path_name)
  } else if (is.na(n4) | n4 == "NA" | n4 == "") {
    l1 <- boots %>% filter(Predictor == n1 & Response == n2) %>% dplyr::select(Iteration, E1 = Estimate)
    l2 <- boots %>% filter(Predictor == n2 & Response == n3) %>% dplyr::select(Iteration, E2 = Estimate)
    l1 %>% inner_join(l2, by = "Iteration") %>%
      mutate(Path_Est = E1 * E2, Predictor = n1, Final_Response = n3, Path_Name = path_name) %>%
      dplyr::select(Iteration, Predictor, Final_Response, Path_Name, Path_Est)
  } else {
    l1 <- boots %>% filter(Predictor == n1 & Response == n2) %>% dplyr::select(Iteration, E1 = Estimate)
    l2 <- boots %>% filter(Predictor == n2 & Response == n3) %>% dplyr::select(Iteration, E2 = Estimate)
    l3 <- boots %>% filter(Predictor == n3 & Response == n4) %>% dplyr::select(Iteration, E3 = Estimate)
    l1 %>% inner_join(l2, by = "Iteration") %>% inner_join(l3, by = "Iteration") %>%
      mutate(Path_Est = E1 * E2 * E3, Predictor = n1, Final_Response = n4, Path_Name = path_name) %>%
      dplyr::select(Iteration, Predictor, Final_Response, Path_Name, Path_Est)
  }
}

calc_orig_path <- function(n1, n2, n3, n4, orig_coefs) {
  get_e <- function(p, r) {
    val <- orig_coefs$Estimate[orig_coefs$Predictor == p & orig_coefs$Response == r]
    if(length(val) == 0) return(NA) else return(val)
  }
  if (is.na(n3) | n3 == "NA" | n3 == "") return(get_e(n1, n2))
  if (is.na(n4) | n4 == "NA" | n4 == "") return(get_e(n1, n2) * get_e(n2, n3))
  return(get_e(n1, n2) * get_e(n2, n3) * get_e(n3, n4))
}

# --- 2B. Master Lapply for Processing Output ---
final_outputs <- lapply(c(3, 28, 40), function(r) {
  
  cat("\n===========================================\n")
  cat("Processing back-transformations for:", r, "km...\n")
  
  orig_coefs <- sem_results_3[[r]][[4]][, c("Predictor", "Response", "Estimate")]
  
  path_map_orig <- path_map %>%
    rowwise() %>%
    mutate(Orig_Scaled_Est = calc_orig_path(Node_1, Node_2, Node_3, Node_4, orig_coefs)) %>%
    ungroup()
  
  # -------------------------------------------------------------------------
  # LOAD AND FILTER RAW BOOTSTRAPS
  # -------------------------------------------------------------------------
  cat("  -> Loading raw bootstraps and filtering non-converged runs...\n")
  raw_boots <- readRDS(paste0("./data/raw_bootstraps_", r, "km.rds"))
  
  # Identify and drop any iteration containing extreme non-converged log-odds
  bad_iterations <- raw_boots %>%
    filter(Response == "pr.new" & abs(Estimate) > 20) %>%
    pull(Iteration) %>%
    unique()
  
  if (length(bad_iterations) > 0) {
    cat("     Dropped", length(bad_iterations), "bad iterations.\n")
    raw_boots <- raw_boots %>% filter(!Iteration %in% bad_iterations)
  }
  
  # -------------------------------------------------------------------------
  # SECTION I: Main Estimates (Direct Paths)
  # -------------------------------------------------------------------------
  # Calculate 95% CIs for direct effects from the filtered bootstraps
  summary_direct_results <- raw_boots %>%
    group_by(Response, Predictor) %>%
    summarise(
      Mean_Estimate = mean(Estimate, na.rm = TRUE),
      Median_Estimate = median(Estimate, na.rm = TRUE),
      Lower_CI_5 = quantile(Estimate, 0.025, na.rm = TRUE),
      Upper_CI_95 = quantile(Estimate, 0.975, na.rm = TRUE),
      SD = sd(Estimate, na.rm = TRUE),
      n_successful = n(),
      .groups = 'drop'
    )
  
  res.table <- sem_results_3[[r]][[4]][, 1:8] %>%
    left_join(summary_direct_results, by = c("Predictor", "Response")) %>% 
    left_join(sd.table, by = c("Predictor" = "var.name")) %>% rename(sd.value_pred = sd.value) %>%
    left_join(sd.table, by = c("Response" = "var.name")) %>% rename(sd.value_resp = sd.value) %>%
    left_join(comp.unit.table, by = c("Predictor" = "var.name")) %>%
    mutate(
      Orig_Unscaled_Est = ifelse(Response == "pr.new",
                                 Estimate / sd.value_pred * comp.unit,
                                 Estimate * (sd.value_resp / sd.value_pred) * comp.unit),
      Lower_Unscaled_Est = ifelse(Response == "pr.new",
                                  Lower_CI_5 / sd.value_pred * comp.unit,
                                  Lower_CI_5 * (sd.value_resp / sd.value_pred) * comp.unit),
      Upper_Unscaled_Est = ifelse(Response == "pr.new",
                                  Upper_CI_95 / sd.value_pred * comp.unit,
                                  Upper_CI_95 * (sd.value_resp / sd.value_pred) * comp.unit),
      
      OR_Est = ifelse(Response == "pr.new", exp(Orig_Unscaled_Est), NA),
      OR_Lower = ifelse(Response == "pr.new", exp(Lower_Unscaled_Est), NA),
      OR_Upper = ifelse(Response == "pr.new", exp(Upper_Unscaled_Est), NA)
    ) %>%
    dplyr::select(-sd.value_pred, -sd.value_resp, -comp.unit)
  
  write.csv(res.table, paste0("./data/EstimatesMain_", r, "km.csv"), row.names = FALSE)
  
  # -------------------------------------------------------------------------
  # SECTION II: Specific Paths & Total Effects
  # -------------------------------------------------------------------------
  cat("  -> Processing Indirect and Total Effects...\n")
  
  results_list <- lapply(1:nrow(path_map), function(i) {
    calc_path_per_iter(path_map$Path_Name[i], path_map$Node_1[i], path_map$Node_2[i], 
                       path_map$Node_3[i], path_map$Node_4[i], raw_boots)
  })
  all_iterations_data <- bind_rows(results_list)
  
  # --- II-A. Specific Paths (Indirect Routes) ---
  summary_specific_paths_final <- all_iterations_data %>%
    group_by(Path_Name, Predictor, Final_Response) %>%
    summarise(
      Lower_CI_Scaled = quantile(Path_Est, 0.025, na.rm = TRUE),
      Upper_CI_Scaled = quantile(Path_Est, 0.975, na.rm = TRUE),
      Valid_Bootstraps = n(), .groups = 'drop'
    ) %>%
    left_join(path_map_orig %>% dplyr::select(Path_Name, Orig_Scaled_Est), by = "Path_Name") %>%
    left_join(sd.table, by = c("Predictor" = "var.name")) %>% rename(sd.value_pred = sd.value) %>%
    left_join(sd.table, by = c("Final_Response" = "var.name")) %>% rename(sd.value_resp = sd.value) %>%
    left_join(comp.unit.table, by = c("Predictor" = "var.name")) %>%
    mutate(
      Orig_Unscaled_Est = ifelse(Final_Response == "pr.new", 
                                 Orig_Scaled_Est / sd.value_pred * comp.unit, 
                                 Orig_Scaled_Est * (sd.value_resp / sd.value_pred) * comp.unit),
      Lower_Unscaled_Est = ifelse(Final_Response == "pr.new", 
                                  Lower_CI_Scaled / sd.value_pred * comp.unit, 
                                  Lower_CI_Scaled * (sd.value_resp / sd.value_pred) * comp.unit),
      Upper_Unscaled_Est = ifelse(Final_Response == "pr.new", 
                                  Upper_CI_Scaled / sd.value_pred * comp.unit, 
                                  Upper_CI_Scaled * (sd.value_resp / sd.value_pred) * comp.unit),
      
      OR_Est = ifelse(Final_Response == "pr.new", exp(Orig_Unscaled_Est), NA),
      OR_Lower = ifelse(Final_Response == "pr.new", exp(Lower_Unscaled_Est), NA),
      OR_Upper = ifelse(Final_Response == "pr.new", exp(Upper_Unscaled_Est), NA)
    ) %>%
    dplyr::select(Path_Name, Predictor, Final_Response, Valid_Bootstraps, 
                  Orig_Scaled_Est, Lower_CI_Scaled, Upper_CI_Scaled,
                  Orig_Unscaled_Est, Lower_Unscaled_Est, Upper_Unscaled_Est,
                  OR_Est, OR_Lower, OR_Upper)
  
  write.csv(summary_specific_paths_final, paste0("./data/Specific_Paths_", r, "km.csv"), row.names = FALSE)
  
  # --- II-B. Total Effects (Sum of all specific routes) ---
  orig_total_effects <- summary_specific_paths_final %>%
    group_by(Predictor, Final_Response) %>%
    summarise(Total_Orig_Scaled_Est = sum(Orig_Scaled_Est, na.rm = TRUE), .groups = 'drop')
  
  summary_total_effects_final <- all_iterations_data %>%
    group_by(Iteration, Predictor, Final_Response) %>%
    summarise(Total_Est_Per_Iter = sum(Path_Est, na.rm = TRUE), .groups = 'drop') %>%
    group_by(Predictor, Final_Response) %>%
    summarise(
      Lower_CI_Scaled = quantile(Total_Est_Per_Iter, 0.025, na.rm = TRUE),
      Upper_CI_Scaled = quantile(Total_Est_Per_Iter, 0.975, na.rm = TRUE),
      Valid_Bootstraps = n(), .groups = 'drop'
    ) %>%
    left_join(orig_total_effects, by = c("Predictor", "Final_Response")) %>%
    left_join(sd.table, by = c("Predictor" = "var.name")) %>% rename(sd.value_pred = sd.value) %>%
    left_join(sd.table, by = c("Final_Response" = "var.name")) %>% rename(sd.value_resp = sd.value) %>%
    left_join(comp.unit.table, by = c("Predictor" = "var.name")) %>%
    mutate(
      Effect_Type = "Total_Effect",
      Total_Unscaled_Est = ifelse(Final_Response == "pr.new", 
                                  Total_Orig_Scaled_Est / sd.value_pred * comp.unit, 
                                  Total_Orig_Scaled_Est * (sd.value_resp / sd.value_pred) * comp.unit),
      Total_Lower_Unscaled = ifelse(Final_Response == "pr.new", 
                                    Lower_CI_Scaled / sd.value_pred * comp.unit, 
                                    Lower_CI_Scaled * (sd.value_resp / sd.value_pred) * comp.unit),
      Total_Upper_Unscaled = ifelse(Final_Response == "pr.new", 
                                    Upper_CI_Scaled / sd.value_pred * comp.unit, 
                                    Upper_CI_Scaled * (sd.value_resp / sd.value_pred) * comp.unit),
      
      Total_OR_Est = ifelse(Final_Response == "pr.new", exp(Total_Unscaled_Est), NA),
      Total_OR_Lower = ifelse(Final_Response == "pr.new", exp(Total_Lower_Unscaled), NA),
      Total_OR_Upper = ifelse(Final_Response == "pr.new", exp(Total_Upper_Unscaled), NA)
    ) %>%
    dplyr::select(Effect_Type, Predictor, Final_Response, Valid_Bootstraps,
                  Total_Orig_Scaled_Est, Lower_CI_Scaled, Upper_CI_Scaled,
                  Total_Unscaled_Est, Total_Lower_Unscaled, Total_Upper_Unscaled,
                  Total_OR_Est, Total_OR_Lower, Total_OR_Upper)
  
  write.csv(summary_total_effects_final, paste0("./data/Total_Effects_", r, "km.csv"), row.names = FALSE)
  
  cat("  -> Output complete!\n")
  
  return(list(
    radius = r,
    main_estimates = res.table,
    specific_paths = summary_specific_paths_final,
    total_effects = summary_total_effects_final
  ))
})

names(final_outputs) <- paste0("Radius_", c(3, 28, 40), "km")