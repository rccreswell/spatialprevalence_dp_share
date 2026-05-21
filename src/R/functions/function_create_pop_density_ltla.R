create_pop_density_ltla <- function(filename){
  pop_density <- read_excel(filename,
                              sheet = "MYE 5", skip = 7) 
  
  pop_density <- pop_density %>%
    select(Code, `2019 people per sq. km`) %>%
    filter(startsWith(Code, "E0")) %>%
    rename(location_fine = Code,
           pop_density = `2019 people per sq. km`) %>%
    mutate(location_fine = str_trim(location_fine)) # two LTLAs have blanks at the end

  saveRDS(pop_density, file = "data/processed/covariates/ltla_pop_density_processed.rds")
  
  return(pop_density)
  
}
