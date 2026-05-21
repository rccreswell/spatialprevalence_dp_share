create_monthly_mobility <- function(filename){
  mob <- read.csv("data/raw/covariates/mobility_lad_weekly.csv")

  monthly_mob <- mob %>%
    mutate(week_date = as.Date(week_date, "%Y-%m-%d"),
           mmyy = paste0(month(week_date), "-", year(week_date))) %>%
    group_by(mmyy, start_lad, end_lad, type) %>%
    summarise(trips = sum(trips, na.rm = TRUE)) %>%
    ungroup() 
  
  # Add columns about population size of LTLAs
  pop <- read_excel("data/raw/covariates/ukpopestimatesmid2020on2021geography.xls", 
                    sheet = "MYE2 - Persons", skip = 7) %>%
    filter(startsWith(Code, "E0")) %>%
    rename(location_fine = Code,
           pop = `All ages`) %>%
    select(location_fine, pop)
  
  monthly_mob <- left_join(monthly_mob, pop, by = c("start_lad" = "location_fine")) %>%
    rename(pop_start = pop)
  monthly_mob <- left_join(monthly_mob, pop, by = c("end_lad" = "location_fine")) %>%
    rename(pop_end = pop)
  
  # REMOVING DIRECTIONAILITY
  # Combine start and end LADs into a single column to consider both directions
  monthly_mob1 <- monthly_mob %>%
    mutate(pair_lad = paste(pmin(start_lad, end_lad), pmax(start_lad, end_lad), sep = "_"),
           type_within = ifelse(start_lad == end_lad, 1, 0)) %>%
    mutate(pop = round((pop_start + pop_end)/2))
  
  # Summarize the data to calculate average trips between LAD pairs
  monthly_mob1 <- monthly_mob1 %>%
    group_by(mmyy, pair_lad) %>%
    summarise(avg_trips = round(mean(trips)),
              type_within = mean(type_within),
              avg_pop = round(mean(pop))) %>%
    ungroup() %>%
    mutate(mob_per_pop = avg_trips/avg_pop)
  

  monthly_mob_within <- monthly_mob1 %>%
    filter(type_within == 1)
  
  monthly_mob_between <- monthly_mob1 %>%
    filter(type_within == 0)
  
  saveRDS(monthly_mob_within, file = "data/processed/covariates/ltla_monthly_mobility_processed_within.rds")
  saveRDS(monthly_mob_between, file = "data/processed/covariates/ltla_monthly_mobility_processed_between.rds")
}