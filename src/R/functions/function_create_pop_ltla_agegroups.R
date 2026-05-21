create_pop_ltla_agegroups <- function(filename){
  pop_agegroups <- read_excel(filename,
                             sheet = "MYE2 - Persons", skip = 7) 
  
  pop_agegroups <- pop_agegroups %>%
    mutate(`00-04` = rowSums(across(c("0", "1", "2", "3", "4")), na.rm = TRUE),,
           `05-11` = rowSums(across(c("5", "6", "7", "8", "9", "10", "11")), na.rm = TRUE),
           `12-17` = rowSums(across(c("12", "13", "14", "15", "16", "17")), na.rm = TRUE),
           `18-24` = rowSums(across(c("18", "19", "20", "21", "22", "23", "24")), na.rm = TRUE),
           `25-34` = rowSums(across(c("25", "26", "27", "28", "29", "30", "31", "32", "33", "34")), na.rm = TRUE),
           `35-44` = rowSums(across(c("35", "36", "37", "38", "39", "40", "41", "42", "43", "44")), na.rm = TRUE),
           `45-54` = rowSums(across(c("45", "46", "47", "48", "49", "50", "51", "52", "53", "54")), na.rm = TRUE),
           `55-64` = rowSums(across(c("55", "56", "57", "58", "59", "60", "61", "62", "63", "64")), na.rm = TRUE),
           `65-74` = rowSums(across(c("65", "66", "67", "68", "69", "70", "71", "72", "73", "74")), na.rm = TRUE),
           `75+` = rowSums(across(c("75", "76", "77", "78", "79", 
                                    "80", "81", "82", "83", "84", "85", "86", "87", "88", "89",
                                    "90+")), na.rm = TRUE)) %>%
    dplyr::select(Code, `00-04`, `05-11`, `12-17`, `18-24`, `25-34`, `35-44`, `45-54`, `55-64`, `65-74`, `75+`) %>%
    pivot_longer(!Code, names_to = "age_group",
                 values_to = "pop",values_drop_na = TRUE) %>%
    filter(startsWith(Code, "E0")) %>%
    rename(location_fine = Code)
  
  saveRDS(pop_agegroups, file = "data/processed/covariates/ltla_pop_agegroups_processed.rds")
  
  return(pop_agegroups)
  
}
