library(dplyr)
library(tidyr)

# load imputed dataset
me_impute_filter = read.csv("./data/alldata.impute.csv") %>%
  # change some of the variables
  rename(mean_ndvi = mean.ndvi, mangrove_cover = mangrove.cover) %>%
  mutate(infected = round(pr * examined, 0),
         pr.new = infected/examined,
         date = as.Date(paste(year_start, month_start), '%Y%j'))

input = me_impute_filter %>% mutate(uninfected = examined - infected) %>%
  dplyr::select(-pr, -pr.new, -examined) %>%
  pivot_longer(cols = c("infected", "uninfected"),
               names_to = "infected", values_to = "weight") %>%
  mutate(infected = ifelse(infected == "infected", 1, 0)) %>%
  filter(weight != 0) %>% 
  dplyr::select(#-contains("stdev"), -contains("sd"), -month_end, -year_end, -date,
                infected, mean_ndvi, mangrove_cover, mangrove.cover.min1,
                coastline.dist, lc_crop, 
                pop_dens_median, anomaly_2t, anomaly_tp, mean_2t, mean_tp, 
                anomaly_2t_6m, anomaly_tp_6m, mean_2t_6m, mean_tp_6m,
                lat, lon, year_start, month_start, method, weight)

write.table(input, "./data/ML_input_file.csv", row.names = F, sep = ",")
