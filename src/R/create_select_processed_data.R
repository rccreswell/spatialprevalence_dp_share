# ------------------------------------------------------------------------------
# Create mmyy dataframe
library(dplyr)
library(lubridate)
dat_og <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")

weeks_monthyear <- tibble(weeks = unique(dat_og$week_date) ) %>%
  mutate(monthyear = format(as.Date(weeks), "%Y-%m"))

# Create mmyy df
mmyy_df <- NA
i=1
for(temp_monthyear in unique(weeks_monthyear$monthyear)){
  print(paste0("Month: ", temp_monthyear))
  temp_weeks <- weeks_monthyear %>% 
    filter(monthyear == temp_monthyear) %>%
    select(weeks) %>%
    pull() %>%
    as.Date()
  
  # Plotting
  if(length(temp_weeks) <= 1)
    print("only one week in this month, so moving on to next month")
  else{
    mmyy <- paste0(month(temp_weeks[1]), "-", year(temp_weeks[1]))
    
    mmyy_df[i] <- mmyy
    
    i = i + 1
  }
}

mmyy_df <- mmyy_df %>% as.data.frame() 
colnames(mmyy_df) <- c("mmyy")

write.csv(mmyy_df, file = "data/processed/mmyy_df.csv")

# ------------------------------------------------------------------------------
# Clean ONS data
library(readxl)
library(dplyr)

dat <- read_excel("data/raw/20230203covid19infectionsurveydatasetsengland.xlsx", 
                  sheet = "UK summary - positivity", skip = 5) %>%
  select(`England\r\nTime period`, 
         `England \r\nEstimated average % of the population testing positive for COVID-19`,
         `England \r\n95% Lower confidence/credible interval`,
         `England \r\n95% Upper confidence/ credible interval`) %>%
  rename(time_period = `England\r\nTime period`,
         mean_prev = `England \r\nEstimated average % of the population testing positive for COVID-19`,
         lower = `England \r\n95% Lower confidence/credible interval`,
         upper = `England \r\n95% Upper confidence/ credible interval`) %>%
  mutate(mean_prev = mean_prev/100,
         lower = lower/100,
         upper = upper/100)

# Subtracting 2 days from ONS dates, to match our week evaluation
dat <- dat %>%
  mutate(start_date = sub("\\ to.*", "", time_period),
         Date_og = as.Date(start_date, "%d %B %Y"),
         Date = Date_og - 2) %>% 
  filter(Date>= "2020-07-25", Date<= "2022-04-02") %>%
  select(-time_period, -start_date)

saveRDS(dat, "data/processed/ons_prev_england.rds")

# ------------------------------------------------------------------------------
# Create monthly mobility processed data
source("src/R/functions/function_create_monthly_mobility.R")
create_monthly_mobility("data/raw/covariates/mobility_lad_weekly.csv")


# ------------------------------------------------------------------------------
# Create distance between LTLAs and if they are neighbours
source("src/R/functions/function_create_distance_between_ltlas_and_neighbours.R")
create_distance_between_ltlas()
