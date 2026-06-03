library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)
library(tidyverse)
source("src/R/functions/function_create_pop_ltla_agegroups.R")
source("src/R/functions/function_create_pop_density_ltla.R")
source("src/R/functions/function_create_imd_ltla.R")
source("src/R/functions/function_utility.R")

# Monthly regression between prob. of being in the same cluster and covariates
do_epianalysis <- function(mmyy, maxIters){

  # mmyy <- "10-2020"
  # maxIters <- 2000
  mmyy_input <- mmyy
  
  # OUTCOME
  # Get prob. of being in the same cluster ready
  prob <- readRDS(paste0("data/processed/prob_clustertrend_assignment_", mmyy, ".rds")) %>%
    mutate(pair_lad = paste(pmin(location_fine1, location_fine2), pmax(location_fine1, location_fine2), sep = "_"))
  
  # RUN ONLY FOR MONTHS WHERE WE HAVE ATLEAST THREE CLUSTERS???
  
  # COVARIATES
  # Get LTLA age group specific population size ready
  # Keep proportion of individuals aged 0-64 in each LTLA
  pop_agegroups <- create_pop_ltla_agegroups("data/raw/covariates/ukpopestimatesmid2020on2021geography.xls") %>%
    mutate(age_group_new = ifelse(age_group %in% c("65-74", "75+"), "65+", "0-64")) %>%
    group_by(location_fine, age_group_new) %>%
    summarise(pop = sum(pop, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(location_fine) %>%
    mutate(pop_ltla = sum(pop),
           prop = pop/pop_ltla) %>%
    ungroup() %>%
    filter(age_group_new == "65+") %>%
    select(location_fine, prop)
  
  write.csv(pop_agegroups, file = "data/processed/covariates/pop_agegroups64.csv")
  
  # Get LTLA population density
  pop_density <- create_pop_density_ltla("data/raw/covariates/ukpopestimatesmid2020on2021geography.xls")
  write.csv(pop_density, file = "data/processed/covariates/pop_density.csv")
  
  # Get LTLA IMD ready
  imd <- create_imd_ltla("data/raw/covariates/File_10_-_IoD2019_Local_Authority_District_Summaries__lower-tier__.xlsx")
  write.csv(imd, file = "data/processed/covariates/imd.csv")
  
  # Get mobility within and between LTLAs ready
  # DIRECTLY USING PROCESSED MOBILITY DATA
  mob_within <- readRDS("data/processed/covariates/ltla_monthly_mobility_processed_within.rds") %>% 
    filter(mmyy == mmyy_input) %>%
    mutate(location_fine = sub("\\_.*", "", pair_lad)) %>%
    select(location_fine, mob_per_pop)
  
  # If there is no mobility data available, exit
  if(nrow(mob_within) == 0){
    output_no <- data.frame("covariate" = NA, "estimate" = NA, 
                            "se" = NA, "month_year" = mmyy_input)
    return(output_no)
  } 
  # If the number of clusters <=3, exit
  clustertrend <- readRDS(paste0("data/processed/all_clustertrend_assignment_", mmyy, ".rds")) %>%
    select(last_col()) %>%
    pull()
  n_clusters <- max(clustertrend)
  
  if(n_clusters <= 3){
    output_no <- data.frame("covariate" = NA, "estimate" = NA, 
                            "se" = NA, "month_year" = mmyy_input)
    return(output_no)
  }else{
    mob_between <- readRDS("data/processed/covariates/ltla_monthly_mobility_processed_between.rds") %>%
      filter(mmyy == mmyy_input)%>%
      select(pair_lad, mob_per_pop)
    
    # Get distance between LTLAs
    # DIRECTLY USING PROCESSED DISTANCE DATA from function_create_distance_between_ltlas_and_neighbours.R
    dist <- readRDS("data/processed/distance_between_ltlas.rds") 
    
    # Get if LTLAs are neighbours
    # DIRECTLY USING PROCESSED NEIGHBOUR DATA from function_create_distance_between_ltlas_and_neighbours.R
    is_neighbour <- readRDS("data/processed/is_neighbours_ltlas.rds")
    
    # Merge all data for the regression
    # # mobility within LTLA 1
    # dat <- left_join(prob, mob_within, by = c("location_fine1" = "location_fine")) %>%
    #   mutate(avg_trips_within1 = avg_trips,
    #          avg_trips = standardize(avg_trips)) %>%
    #   rename(mobility_within1 = avg_trips)
    # # mobility within LTLA 2
    # dat <- left_join(dat, mob_within, by = c("location_fine2" = "location_fine")) %>%
    #   mutate(avg_trips_within2 = avg_trips,
    #          avg_trips = standardize(avg_trips)) %>%
    #   rename(mobility_within2 = avg_trips)
    # # mobility between LTLA 1 and 2
    # dat <- left_join(dat, mob_between, by = c("pair_lad" = "pair_lad")) %>%
    #   mutate(avg_trips_between = avg_trips,
    #          avg_trips = standardize(avg_trips)) %>%
    #   rename(mobility_between = avg_trips)
    # mobility within LTLA 1
    dat <- left_join(prob, mob_within, by = c("location_fine1" = "location_fine")) %>%
      rename(mob_per_pop_within1 = mob_per_pop)
    # mobility within LTLA 2
    dat <- left_join(dat, mob_within, by = c("location_fine2" = "location_fine")) %>%
      rename(mob_per_pop_within2 = mob_per_pop)
    # mobility between LTLA 1 and 2
    dat <- left_join(dat, mob_between, by = c("pair_lad" = "pair_lad")) %>%
      rename(mob_per_pop_between = mob_per_pop)
    
    # IMD of LTLA 1 
    dat <- left_join(dat, imd, by = c("location_fine1" = "location_fine")) %>%
      mutate(imd_per10_1 = imd/10) %>%
      select(-imd)
    # IMD of LTLA 2
    dat <- left_join(dat, imd, by = c("location_fine2" = "location_fine")) %>%
      mutate(imd_per10_2 = imd/10) %>%
      select(-imd)
    
    # Pop density of LTLA 1 
    dat <- left_join(dat, pop_density, by = c("location_fine1" = "location_fine")) %>%
      mutate(pop_density_per1000_1 = pop_density/1000) %>%
      select(-pop_density)
    # Pop density of LTLA 2
    dat <- left_join(dat, pop_density, by = c("location_fine2" = "location_fine")) %>%
      mutate(pop_density_per1000_2 = pop_density/1000) %>%
      select(-pop_density)
    
    # Pop proportion aged 0-64 in LTLA 1
    dat <- left_join(dat, pop_agegroups, by = c("location_fine1" = "location_fine")) %>%
      rename(prop1 = prop)
    # Pop proportion aged 0-64 in LTLA 2
    dat <- left_join(dat, pop_agegroups, by = c("location_fine2" = "location_fine")) %>%
      rename(prop2 = prop) 
    
    # Distance between LTLAs
    dat <- merge(dat, dist, by = c("location_fine1" = "location_fine1", 
                                   "location_fine2" = "location_fine2")) %>%
      mutate(distance_km = distance_m/1000,
             distance_km_per100 = distance_km/100) %>%
      select(-distance_km, -distance_m)
    # Are LTLAs neighbours
    dat <- merge(dat, is_neighbour, by = c("location_fine1" = "location_fine1", 
                                           "location_fine2" = "location_fine2"))
    
    # Remove NA's and make probabilities -> counts
    dat <- dat %>%
      na.omit() %>%
      mutate(y = prob*(maxIters/2),
             n = (maxIters/2))
    
    # Create  new variables of difference between pairs for IMD, pop density and prop
    # dat <- dat %>%
    #   mutate(imd_per10_diff = abs(imd_per10_1 - imd_per10_2),
    #          pop_density_per100_diff = abs(pop_density_per100_1 - pop_density_per100_2),
    #          prop_diff = abs(prop1 - prop2),
    #          mob_per_pop_within_diff = abs(mob_per_pop_within1 - mob_per_pop_within2))
    
    # TRY logging prop difference
    dat <- dat %>%
      mutate(imd_per10_diff = abs(imd_per10_1 - imd_per10_2),
             pop_density_per1000_diff = abs(pop_density_per1000_1 - pop_density_per1000_2),
             prop_diff = log(abs(prop1 - prop2)),
             mob_per_pop_within_diff = abs(mob_per_pop_within1 - mob_per_pop_within2))
    
    # First try models in R (lots of zeroes in the outcome)
    # Logistic
    m1 <- glm(cbind(round(y), round(n - y)) ~ mob_per_pop_within_diff + 
                mob_per_pop_between + 
                imd_per10_diff +
                pop_density_per1000_diff +
                prop_diff +
                distance_km_per100 +
                is_neighbour, data = dat, family = "binomial")
    summary(m1)
    
    # Linear
    m2 <- lm(prob ~ mob_per_pop_within_diff + 
               mob_per_pop_between + 
               imd_per10_diff +
               pop_density_per1000_diff +
               prop_diff +
               distance_km_per100 +
               is_neighbour, data = dat)  
    summary(m2)
    p_pred <- predict(m2)
    plot(dat$prob, p_pred)
    
    # Logit transformed outcome
    # m3 <- lm(log((prob)/(1-prob)) ~ mob_per_pop_within_diff +
    #            mob_per_pop_between +
    #            imd_per10_diff +
    #            pop_density_per1000_diff +
    #            prop_diff +
    #            distance_km_per100 +
    #            is_neighbour, data = dat)
    # summary(m3)
    # p_pred <- predict(m3)
    # plot(dat$prob, p_pred)
    
    # zero-one inflated beta regression
    m_zoib <- gamlss(prob ~ mob_per_pop_within_diff +
                       mob_per_pop_between +
                       imd_per10_diff +
                       pop_density_per1000_diff +
                       prop_diff +
                       distance_km_per100 +
                       is_neighbour,
                     family = BEINF,
                     nu.fo = ~ mob_per_pop_within_diff +
                       mob_per_pop_between +
                       imd_per10_diff +
                       pop_density_per1000_diff +
                       prop_diff +
                       distance_km_per100 +
                       is_neighbour,
                     tau.fo = ~ mob_per_pop_within_diff +
                       mob_per_pop_between +
                       imd_per10_diff +
                       pop_density_per1000_diff +
                       prop_diff +
                       distance_km_per100 +
                       is_neighbour,
                     data = dat)

    summary(m_zoib)
    # Pseudo-R² for the BEINF model
    pR2(m_zoib)
    pR2(m2)

    AIC(m1, m_zoib, m2)

    # Reshape data to long format
    dat_long <- dat %>%
      select(prob,
             mob_per_pop_within_diff,
             mob_per_pop_between,
             imd_per10_diff,
             pop_density_per1000_diff,
             prop_diff,
             distance_km_per100,
             is_neighbour)
    dat_long <- pivot_longer(dat_long, cols = -prob, names_to = "variable", values_to = "value")

    # Plot outcome vs each variable in facets
    ggplot(dat_long, aes(x = value, y = prob)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      facet_wrap(~variable, scales = "free_x") +
      theme_minimal() +
      labs(title = paste0("Probability of co-clustering vs other variables, ", mmyy_input),
           x = "Predictor",
           y = "Prob of co-clustering")
    
    
    output <- summary(m2)$coef[,c("Estimate", "Std. Error")] %>% 
      cbind(month_year = mmyy_input) %>%
      as.data.frame() %>%
      rownames_to_column(var = "covariate")
    
    colnames(output) <- c("covariate", "estimate", 
                          "se", "month_year")
    
    return(output)
  }
}


# Prepare the data for stan
# Original dataframe
# X <- dat %>% select(mobility_within1, mobility_within2, mobility_between,
#                     imd1, imd2, prop1, prop2)
# # Dataframe with difference ad average between 2 LTLAs' imd, prop, and mobility (1)
# X <- dat %>% select(mobility_within_av1, mobility_within_diff1, mobility_between,
#                     imd_av, imd_diff, prop_av, prop_diff)
# X <- cbind(intercept = 1, X)
# 
# data_list <- list(
#   N = nrow(dat),
#   K = ncol(X),
#   X = X,
#   y = dat$y,
#   n = dat$n
# )
# 
# Compile the Stan model
# stan_model <- stan_model("src/stan/epi_regression.stan")
# 
# # Fit the model
# fit <- sampling(stan_model, data = data_list, chains = 4, iter = 1000)
# saveRDS(fit, file = paste0("data/processed/stan_fit_",
#                            mmyy, ".rds"))
# # stan_fit_old_mmyy is original variables and stan_fit_mmyy is with av and diff variables
# print(fit)
# 
# # Compare the observed and predicted number of time LTLAs were in the same cluster
# post_dat <- extract(fit)[["y_pred"]]
# quants <- c(0.025, 0.50, 0.0975)
# post_dat <- apply(post_dat, 2 , quantile , probs = quants , na.rm = TRUE ) %>% 
#   t() %>%
#   as.data.frame()
# post_dat <- cbind(dat, post_dat)
# 
# ggplot(data = post_dat, aes(x = y, y = `50%`)) +
#   geom_point()
# 
# 
# post_dat <- extract(fit)[["p_pred"]]
# quants <- c(0.025, 0.50, 0.0975)
# post_dat <- apply(post_dat, 2 , quantile , probs = quants , na.rm = TRUE ) %>% 
#   t() %>%
#   as.data.frame()
# post_dat <- cbind(dat, post_dat)
# 
# ggplot(data = post_dat, aes(x = prob, y = `50%`)) +
#   geom_jitter() +
#   geom_smooth(method = "lm") +
#   geom_abline()




