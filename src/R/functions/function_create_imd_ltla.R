create_imd_ltla <- function(filename){
  ltla_imd <- read_excel(filename,
                             sheet = "IMD")
  
  # Use continuous IMD
  ltla_imd <- ltla_imd %>%
    # mutate(IMD_decile = ntile(`IMD - Rank of average score`, 10)) %>%
    mutate(imd = `IMD - Average score`) %>%
    rename(location_fine = `Local Authority District code (2019)`) %>%
    select(location_fine, imd)
  
  saveRDS(ltla_imd, file = "data/processed/covariates/ltla_imd_processed.rds")
  
  return(ltla_imd)
  
}