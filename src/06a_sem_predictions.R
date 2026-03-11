# ==============================================================================
# SEM PREDICTION AND PLOTTING SCRIPT (V4 - Integrated & '=' Assignment)
# Plots raw data against glmmPQL trendlines and bootstrapped SEM effects
# ==============================================================================

library(ggplot2)
library(dplyr)
library(nlme)
library(cowplot)

# ------------------------------------------------------------------------------
# 1. LOAD DATA & MODELS
# ------------------------------------------------------------------------------
source("./src/sem_functions.R")

all.data.raw = read.csv("./data/alldata.full.csv") %>% data.prep()
all.data.impute = read.csv("./data/alldata.impute.csv") %>% data.prep()

pr.fct = function(r) {
  cat("\n==================================================\n")
  cat("STARTING RADIUS:", r, "km\n")
  cat("==================================================\n")
  
  # 1. Data Loading & Preparation
  input_set = all.data.impute %>% 
    filter(buffer == r)
  
  # 2. Initial Master Models
  # MODEL A: Prevalence (Binomial glmmPQL with spatial autocorrelation)
  glmmPQL(
    fixed = as.formula(paste(deparse(f[[1]]), collapse ="")),
    random = ~ 1 | group,
    correlation = corExp(1, form = ~ lat + lon | group, nugget = TRUE),
    data = input_set,
    weights = input_set$examined, 
    family = binomial(link = "logit"),
    control = lmeControl(msMaxIter = 1000, msMaxEval = 1000), 
    verbose = FALSE
  )
}

model.pr.3 = pr.fct(3)
model.pr.40 = pr.fct(40)

# Define the scaled datasets used for the models
input_set_3 = all.data.impute %>% filter(buffer == 3)
input_set_40 = all.data.impute %>% filter(buffer == 40)

# Load the saved SEM results for the 40km scale to extract betas programmatically
# Make sure the file paths match your working directory
res_main_40 = read.csv("./data/EstimatesMain_40km.csv")
res_ind_40 = read.csv("./data/Specific_Paths_40km.csv")
res_tot_40 = read.csv("./data/Total_Effects_40km.csv")

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

# Function A: Standard Prediction with 95% CIs (for Plot B)
get_glmmPQL_preds = function(model, var_name, scaled_data, raw_data, scale_factor = 1) {
  beta = fixef(model)
  V = vcov(model)
  
  x_seq = seq(min(scaled_data[[var_name]], na.rm = TRUE), 
              max(scaled_data[[var_name]], na.rm = TRUE), length.out = 100)
  
  X_df = data.frame(matrix(0, nrow = 100, ncol = length(beta)))
  colnames(X_df) = names(beta)
  X_df[[var_name]] = x_seq
  X_df$'(Intercept)' = 1 
  
  X = as.matrix(X_df)
  eta = X %*% beta
  se_eta = sqrt(diag(X %*% V %*% t(X)))
  
  inv_logit = function(x) exp(x) / (1 + exp(x))
  
  raw_mean = mean(raw_data[[var_name]], na.rm = TRUE)
  raw_sd = sd(raw_data[[var_name]], na.rm = TRUE)
  
  data.frame(
    plot_x = (x_seq * raw_sd + raw_mean) * scale_factor,
    fit_pr = inv_logit(eta),
    lwr_pr = inv_logit(eta - 1.96 * se_eta),
    upr_pr = inv_logit(eta + 1.96 * se_eta)
  )
}

# Function B: Compound Prediction with Bootstrapped CIs (for Plot C)
get_compound_preds = function(model, var_name, scaled_data, raw_data, 
                              beta_direct, beta_indirect, beta_total,
                              ci_direct, ci_indirect, ci_total, 
                              scale_factor = 1) {
  
  beta_0 = fixef(model)["(Intercept)"]
  
  x_seq = seq(min(scaled_data[[var_name]], na.rm = TRUE), 
              max(scaled_data[[var_name]], na.rm = TRUE), length.out = 100)
  
  # Calculate center lines (log-odds)
  eta_direct = beta_0 + (beta_direct * x_seq)
  eta_indirect = beta_0 + (beta_indirect * x_seq)
  eta_total = beta_0 + (beta_total * x_seq)
  
  # Calculate bounds (pmin/pmax properly handles negative slopes crossing the origin)
  eta_dir_lwr = pmin(beta_0 + ci_direct[1]*x_seq, beta_0 + ci_direct[2]*x_seq)
  eta_dir_upr = pmax(beta_0 + ci_direct[1]*x_seq, beta_0 + ci_direct[2]*x_seq)
  
  eta_ind_lwr = pmin(beta_0 + ci_indirect[1]*x_seq, beta_0 + ci_indirect[2]*x_seq)
  eta_ind_upr = pmax(beta_0 + ci_indirect[1]*x_seq, beta_0 + ci_indirect[2]*x_seq)
  
  eta_tot_lwr = pmin(beta_0 + ci_total[1]*x_seq, beta_0 + ci_total[2]*x_seq)
  eta_tot_upr = pmax(beta_0 + ci_total[1]*x_seq, beta_0 + ci_total[2]*x_seq)
  
  inv_logit = function(x) exp(x) / (1 + exp(x))
  
  raw_mean = mean(raw_data[[var_name]], na.rm = TRUE)
  raw_sd = sd(raw_data[[var_name]], na.rm = TRUE)
  plot_x = (x_seq * raw_sd + raw_mean) * scale_factor
  
  data.frame(
    plot_x = rep(plot_x, 3),
    fit_pr = c(inv_logit(eta_direct), inv_logit(eta_indirect), inv_logit(eta_total)),
    lwr_pr = c(inv_logit(eta_dir_lwr), inv_logit(eta_ind_lwr), inv_logit(eta_tot_lwr)),
    upr_pr = c(inv_logit(eta_dir_upr), inv_logit(eta_ind_upr), inv_logit(eta_tot_upr)),
    effect_type = factor(rep(c("Direct", "Indirect", "Total"), each = 100),
                         levels = c("Total", "Direct", "Indirect"))
  )
}

# ------------------------------------------------------------------------------
# 3. EXTRACT SEM ESTIMATES & GENERATE PREDICTIONS
# ------------------------------------------------------------------------------

# --- Plot B: NDVI Prediction (3km) ---
pred_ndvi = get_glmmPQL_preds(
  model = model.pr.3, 
  var_name = "mean_ndvi", 
  scaled_data = input_set_3, 
  raw_data = all.data.raw %>% filter(buffer == 3),
  scale_factor = 1
)

# --- Plot C: Mangrove Cover Compound Predictions (40km) ---

# 1. Direct Effect
my_beta_direct = res_main_40 %>% filter(Predictor == "mangrove_cover" & Response == "pr.new") %>% pull(Estimate)
ci_direct = c(res_main_40 %>% filter(Predictor == "mangrove_cover" & Response == "pr.new") %>% pull(Lower_CI_5),
              res_main_40 %>% filter(Predictor == "mangrove_cover" & Response == "pr.new") %>% pull(Upper_CI_95))

# 2. Total Effect
my_beta_total = res_tot_40 %>% filter(Predictor == "mangrove_cover" & Final_Response == "pr.new") %>% pull(Total_Orig_Scaled_Est)
ci_total = c(res_tot_40 %>% filter(Predictor == "mangrove_cover" & Final_Response == "pr.new") %>% pull(Lower_CI_Scaled),
             res_tot_40 %>% filter(Predictor == "mangrove_cover" & Final_Response == "pr.new") %>% pull(Upper_CI_Scaled))

# 3. Indirect Effect (NOTE: Update "Mangrove_NDVI_Mediation" to your actual exact path name!)
my_beta_indirect = res_ind_40 %>% filter(Path_Name == "Cover_via_NDVI") %>% pull(Orig_Scaled_Est)
ci_indirect = c(res_ind_40 %>% filter(Path_Name == "Cover_via_NDVI") %>% pull(Lower_CI_Scaled),
                res_ind_40 %>% filter(Path_Name == "Cover_via_NDVI") %>% pull(Upper_CI_Scaled))

pred_compound_cover = get_compound_preds(
  model = model.pr.40, 
  var_name = "mangrove_cover", 
  scaled_data = input_set_40, 
  raw_data = all.data.raw %>% filter(buffer == 40),
  beta_direct = my_beta_direct, ci_direct = ci_direct,
  beta_indirect = my_beta_indirect, ci_indirect = ci_indirect,
  beta_total = my_beta_total, ci_total = ci_total,
  scale_factor = 100
)

# ------------------------------------------------------------------------------
# 4. PLOTTING
# ------------------------------------------------------------------------------

# Plot B: NDVI with Standard CI Ribbon
plot.b = ggplot(data = all.data.raw %>% filter(buffer == 3), 
                aes(x = mean_ndvi, y = pr)) +
  geom_point(aes(size = examined), alpha = 0.5) +
  geom_ribbon(data = pred_ndvi, aes(x = plot_x, y = fit_pr, ymin = lwr_pr, ymax = upr_pr), 
              fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = pred_ndvi, aes(x = plot_x, y = fit_pr), 
            color = "blue", size = 1.2, inherit.aes = FALSE) +
  xlab("Mangrove NDVI") + ylab("Malaria prevalence") +
  theme_bw() + 
  theme(text = element_text(size = 15), legend.position = "none")

# Plot C: Mangrove Cover with Decomposed SEM Lines AND Ribbons
plot.c = ggplot(data = all.data.raw %>% filter(buffer == 40), 
                aes(x = mangrove_cover*100, y = pr)) +
  geom_point(aes(size = examined), alpha = 0.3, color = "grey40") +
  
  # CIs (alpha set low at 0.1 to avoid making the plot too dense)
  geom_ribbon(data = pred_compound_cover, 
              aes(x = plot_x, ymin = lwr_pr, ymax = upr_pr, fill = effect_type), 
              alpha = 0.1, inherit.aes = FALSE) +
  
  # Trendlines
  geom_line(data = pred_compound_cover, 
            aes(x = plot_x, y = fit_pr, color = effect_type, linetype = effect_type), 
            size = 1.2, inherit.aes = FALSE) +
  
  # Color and styling mapping
  scale_color_manual(values = c("Total" = "black", "Direct" = "blue", "Indirect" = "red"), name = "Effect Pathway") +
  scale_fill_manual(values = c("Total" = "black", "Direct" = "blue", "Indirect" = "red"), name = "Effect Pathway") +
  scale_linetype_manual(values = c("Total" = "solid", "Direct" = "dashed", "Indirect" = "dotted"), name = "Effect Pathway") +
  
  xlab("Mangrove land cover [%]") + ylab("Malaria prevalence") +
  theme_bw() + 
  theme(text = element_text(size = 15),
        legend.position = c(0.75, 0.8), 
        legend.background = element_rect(fill = alpha("white", 0.8), color = "grey50"))

# Combine Plots
plot.joint = plot_grid(plot.b, plot.c, labels = c('B', 'C'), 
                       label_size = 20, rel_widths = c(1, 1))

print(plot.joint)

# ------------------------------------------------------------------------------
# 5. EXPORT PLOTS
# ------------------------------------------------------------------------------

ps.options(family = "Arial")

ggsave(filename="./Figures/mangroveVars_points_v3.svg", plot = plot.joint, device = "svg", width = 270, height = 135, units = "mm")
ggsave(filename="./Figures/mangroveVars_points_v3.png", plot = plot.joint, device = "png", width = 270, height = 135, units = "mm", dpi = 300)
ggsave(filename="./Figures/mangroveVars_points_v3.pdf", plot = plot.joint, device = cairo_pdf, width = 270, height = 135, units = "mm")

cairo_ps(filename="./Figures/mangroveVars_points_v3.eps", width = 10.63, height = 5.315, 
         pointsize = 10, fallback_resolution = 2400)
print(plot.joint)
dev.off()