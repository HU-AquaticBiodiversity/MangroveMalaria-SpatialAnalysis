library(dplyr)
library(parallel)
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

##----------------------------------------------------------##
## Preparation for models
##----------------------------------------------------------##
# detect the number of cores
numberOfCores <- 25

data.prep = function(x) {
  rename(x, mean_ndvi = mean.ndvi, mangrove_cover = mangrove.cover) %>%
    mutate(infected = round(pr * examined, 0),
           pr.new = infected/examined,
           date = as.Date(paste(year_start, month_start), '%Y%j'),
           anomaly_2t_sqr = scale(anomaly_2t^2), 
           anomaly_tp_sqr = scale(anomaly_tp^2),
           mean_2t_sqr = scale(mean_2t^2), 
           mean_tp_sqr = scale(mean_tp^2), 
           anomaly_2t_6m_sqr = scale(anomaly_2t_6m^2),
           anomaly_tp_6m_sqr = scale(anomaly_tp_6m^2),
           mean_2t_6m_sqr = scale(mean_2t_6m^2), 
           mean_tp_6m_sqr = scale(mean_tp_6m^2), 
           anomaly_2t_sqr = scale(anomaly_2t^2),
           # create group variable for random effect and spatial
           #   autocorrelation term
           group = factor(paste0(year_start, month_start, method)))
  }

# load final dataset and prep
all.data = read.csv("./data/alldata.impute.csv") %>% data.prep()

# load small dataset and prep
all.data.small = read.csv("./data/alldata.small.csv") %>% data.prep()

# define starting formulas of intially hypothesis model structure
start.formulas = list(
  as.formula("pr.new ~ mean_ndvi + mangrove_cover + coastline.dist +
  lc_crop + pop_dens_median +
  anomaly_2t + anomaly_tp +
  mean_2t + mean_tp +
  anomaly_2t_6m + anomaly_tp_6m +
  mean_2t_6m + mean_tp_6m"),
  as.formula("mean_ndvi ~ mangrove_cover + pop_dens_median +
  anomaly_2t + anomaly_tp +
  mean_2t + mean_tp +
  anomaly_2t_6m + anomaly_tp_6m +
  mean_2t_6m + mean_tp_6m"),
  as.formula("mangrove_cover ~ mangrove.cover.min1 + lc_crop + 
  pop_dens_median + coastline.dist +
  anomaly_2t + anomaly_tp +
  mean_2t + mean_tp +
  anomaly_2t_6m + anomaly_tp_6m +
  mean_2t_6m + mean_tp_6m")
)

weather.formula = "+ anomaly_2t_sqr + anomaly_tp_sqr + mean_2t_sqr + 
mean_tp_sqr + anomaly_2t_6m_sqr + anomaly_tp_6m_sqr + mean_2t_6m_sqr +
mean_tp_6m_sqr"

##---------------------------------------------------##
## Function for fitting models for r between 1-50 km
##
##  - input data: which data frame is used as input (complete vs. small dataset)
##  - input formula: starting formulas that are being updated. Either
##    'starting_formulas' defined above or output of previous optimisation step.
##  - ptd, pta: paths to drop (ptd) and paths to add (pta), need to be defined
##    prior by "SEM interpretation" script on model output from previous step.
##  - weather_sqr: should quadratic terms be added for the weather variables.
##    This part only needs to be T for the first two steps, when
##    input_formulas = start.formulas
##---------------------------------------------------##
model_fit = function(input_data, input_formula,
                     ptd = NULL, pta = NULL, 
                     weather_sqr = F) {
  
  # vector of names of dependent variables
  dep = c("pr.new", "mean_ndvi", "mangrove_cover")
  
  # create formulas for models
  f = lapply(1:3, function(s) {
    
    # check whether ptd is NULL, if so make it NA
    if(is.null(ptd)) {ptd.f = NA} else {
      # check whether ptd is empty for selected r
      if(nrow(filter(ptd, dependent == dep[s])) == 0) {
        # if so make it NA
        ptd.f = NA} else {
          # otherwise, create string of effects to remove
          ptd.f = ptd %>% filter(dependent == dep[s]) %>%
            dplyr::select("independent") %>% unlist()
          ptd.f = paste0("- ", paste(ptd.f, collapse = " - "))
        }
      }
      
    # do the same with pta  
    if(is.null(pta)) {pta.f = NA} else {
      if(nrow(filter(pta, dependent == dep[s])) == 0) {
        pta.f = NA} else {
          pta.f = pta %>% filter(dependent == dep[s]) %>%
            dplyr::select("independent") %>% unlist()
          pta.f = paste0("+ ", paste(pta.f, collapse = " + "))
        }
      }
    
    # check whether weather variables have to be added as quadratic effects
    if(weather_sqr == F){
      # if not: drop/add paths to starting formula
      update(input_formula[[s]], 
             paste0(na.omit(c("~.", ptd.f, pta.f)), collapse = ""))
    } else {
      # otherwise ...
      update(
        # .... first add weather_sqr variables
        update(input_formula[[s]], 
               paste0(na.omit(c("~.", weather.formula)), collapse = "")),
        # ... then drop and add paths
        paste0(na.omit(c("~.", ptd.f, pta.f)), collapse = "")
      )
    }
  })
  
  ## run model for r between 1 and 50 km in 1-km intervals
  res = mclapply(1:50, function(r){

    # filter input data by r
    input_set <- input_data %>% filter(buffer == r)

    ## MODEL A: prevalence
    ##          --> needs to be defined with binomial distribution
    model.pr.1 = glmmPQL(
      # Formulas cannot be read as variables with lapply and glmmPQL
      # Therefore, they need to be deparsed and back-transformed into
      # formulas directly inside of glmmPQL function to work.
      as.formula(paste(deparse(f[[1]]), collapse ="")),
      # define spatial and temporal autocorrelation term
      random = ~ 1 | group,
      correlation = corExp(1, form = ~ lat + lon | group, nugget = T),
      data = input_set,
      # set weights for binomial distribution
      weights = input_set$examined,
      # iterations were extended to assure conversion
      control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
      # important: prevalence data follow binomial structure
      family = binomial(link = "logit"),
      # to turn off automatic output
      verbose = F)
    
    # MODEL B: mangrove NDVI
    model.vi.1 = update(model.pr.1,
                        as.formula(paste(deparse(f[[2]]), collapse ="")),
                        # important: prevalence data follow binomial structure
                        family = gaussian
                        )

    # MODEL C: mangrove cover
    model.mc.1 <- update(model.vi.1,
                         as.formula(paste(deparse(f[[3]]), collapse =""))
                         )
    
    ### SEM ###
    mangrove_sem = psem(
      model.pr.1,
      model.vi.1,
      model.mc.1,
      data = input_set
    )
    
    # save SEM output
    a = list(plot(mangrove_sem, return = T), 
             dSep(mangrove_sem, .progressBar = F), 
             fisherC(mangrove_sem),
             coefs(mangrove_sem),
             AIC_psem(mangrove_sem))
    
    # print that round of calculations has been terminated
    print(paste0(r, "/50 Done!"))
    
    a
  }, mc.preschedule = FALSE, mc.cores = numberOfCores)
  
  # output results for model formulas and SEM
  list(f, res)
}

## NOTE: These models should be run in a linux environment. Otherwise, I
##       recommend using 'parLapply'.
## NOTE: Don't forget to change 'mc.cores' to the number of available cores in
##       your environment!

##----------------------------------------------------------##
## STEP 1
##----------------------------------------------------------##
system.time({
  sem_results_1 = model_fit(all.data, start.formulas)[[2]]
})

## save for export and viewing on local environment ##
save(sem_results_1, file = "./data/sem_results_1.Rdata")

system.time({
  sem_results_small_1 = model_fit(all.data.small, start.formulas)[[2]]
}) 

## save for export and viewing on local environment ##
save(sem_results_small_1, file = "./data/sem_results_small_1.Rdata")

system.time({
  sem_results_1_sqr = model_fit(all.data, start.formulas, weather_sqr = T)[[2]]
}) 

## save for export and viewing on local environment ##
save(sem_results_1_sqr, file = "./data/sem_results_1_sqr.Rdata")

##----------------------------------------------------------##
## STEP 2
##----------------------------------------------------------##
## full data ##
system.time({
  output = model_fit(all.data,
                            start.formulas,
                            ptd = read.csv("./data/sem_results_1_ptd.csv") %>%
                              filter(independent != "mangrove_cover"),
                            pta = read.csv("./data/sem_results_1_pta.csv"))
  sem_results_2 = output[[2]]
  start.formulas.3 = output[[1]]
})

# save for export and viewing on local environment #
save(sem_results_2, file = "./data/sem_results_2.Rdata")
save(start.formulas.3, file = "./data/start.formulas.3.Rdata")

## small data set
system.time({
  output = model_fit(all.data.small,
                     start.formulas,
                     ptd = read.csv("./data/sem_results_small_1_ptd.csv") %>%
                       filter(independent != "mangrove_cover"),
                     pta = read.csv("./data/sem_results_small_1_pta.csv"))
  sem_results_small_2 = output[[2]]
  start.formulas.3.small = output[[1]]
}) 

# save for export and viewing on local environment #
save(sem_results_small_2, file = "./data/sem_results_small_2.Rdata")
save(start.formulas.3.small, file = "./data/start.formulas.3.small.Rdata")

## full data with quadratic weather variables
system.time({
  output = model_fit(all.data,
                                start.formulas, 
                                ptd = read.csv("./data/sem_results_1_sqr_ptd.csv"),
                                pta = read.csv("./data/sem_results_1_sqr_pta.csv"),
                                weather_sqr = T)
  sem_results_2_sqr = output[[2]]
  start.formulas.3.sqr = output[[1]]
}) 

## save for export and viewing on local environment ##
save(sem_results_2_sqr, file = "./data/sem_results_2_sqr.Rdata")
save(start.formulas.3.sqr, file = "./data/start.formulas.3.sqr.Rdata")

##----------------------------------------------------------##
## STEP 3
##----------------------------------------------------------##

# load formulas of previous model
load("./data/start.formulas.3.Rdata")
load("./data/start.formulas.3.small.Rdata")
load("./data/start.formulas.3.sqr.Rdata")

system.time({
  output = model_fit(all.data,
                     start.formulas.3,
                     ptd = read.csv("./data/sem_results_2_ptd.csv"),
                     pta = read.csv("./data/sem_results_2_pta.csv"))
  sem_results_3 = output[[2]]
  start.formulas.4 = output[[1]]
})

## save for export and viewing on local environment ##
save(sem_results_3, file = "./data/sem_results_3.Rdata")
save(start.formulas.4, file = "./data/start.formulas.4.Rdata")

system.time({
  output = model_fit(all.data.small,
                     start.formulas.3.small,
                     ptd = read.csv("./data/sem_results_small_2_ptd.csv"),
                     pta = read.csv("./data/sem_results_small_2_pta.csv"))
  sem_results_small_3 = output[[2]]
  start.formulas.4.small = output[[1]]
})

## save for export and viewing on local environment ##
save(sem_results_small_3, file = "./data/sem_results_small_3.Rdata")
save(start.formulas.4.small, file = "./data/start.formulas.4.small.Rdata")

system.time({
  output = model_fit(input_data = all.data,
                     input_formula = start.formulas.3.sqr,
                     ptd = read.csv("./data/sem_results_2_sqr_ptd.csv"),
                     pta = read.csv("./data/sem_results_2_sqr_pta.csv"))
  sem_results_3_sqr = output[[2]]
  start.formulas.4.sqr = output[[1]]
})

## save for export and viewing on local environment ##
save(sem_results_3_sqr, file = "./data/sem_results_3_sqr.Rdata")
save(start.formulas.4.sqr, file = "./data/start.formulas.4.sqr.Rdata")

##----------------------------------------------------------##
## STEP 4
##----------------------------------------------------------##

# load formulas of previous model
load("./data/start.formulas.4.sqr.Rdata")

system.time({
  output = model_fit(input_data = all.data,
                     input_formula = start.formulas.4.sqr,
                     ptd = read.csv("./data/sem_results_3_sqr_ptd.csv"),
                     pta = read.csv("./data/sem_results_3_sqr_pta.csv"))
  sem_results_4_sqr = output[[2]]
  start.formulas.5.sqr = output[[1]]
})

## save for export and viewing on local environment ##
save(sem_results_4_sqr, file = "./data/sem_results_4_sqr.Rdata")
save(start.formulas.5.sqr, file = "./data/start.formulas.5.sqr.Rdata")

##----------------------------------------------------------##
## STEP 5
##----------------------------------------------------------##

# load formulas of previous model
load("./data/start.formulas.5.sqr.Rdata")

system.time({
  output = model_fit(input_data = all.data,
                     input_formula = start.formulas.5.sqr,
                     ptd = read.csv("./data/sem_results_4_sqr_ptd.csv"),
                     pta = read.csv("./data/sem_results_4_sqr_pta.csv"))
  sem_results_5_sqr = output[[2]]
  start.formulas.6.sqr = output[[1]]
})

## save for export and viewing on local environment ##
save(sem_results_5_sqr, file = "./data/sem_results_5_sqr.Rdata")
save(start.formulas.6.sqr, file = "./data/start.formulas.6.sqr.Rdata")

##----------------------------------------------------------##
## STEP 6
##----------------------------------------------------------##

# load formulas of previous model
load("./data/start.formulas.6.sqr.Rdata")

system.time({
  output = model_fit(input_data = all.data,
                     input_formula = start.formulas.6.sqr,
                     ptd = read.csv("./data/sem_results_5_sqr_ptd.csv"),
                     pta = read.csv("./data/sem_results_5_sqr_pta.csv"))
  sem_results_6_sqr = output[[2]]
  start.formulas.7.sqr = output[[1]]
})

## save for export and viewing on local environment ##
save(sem_results_6_sqr, file = "./data/sem_results_6_sqr.Rdata")
save(start.formulas.7.sqr, file = "./data/start.formulas.7.sqr.Rdata")