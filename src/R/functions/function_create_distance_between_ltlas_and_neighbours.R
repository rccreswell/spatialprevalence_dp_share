library(sf) # to read in the shapefile
library(tmap) # to plot the map
library(geosphere) # for distnaces between LTLAs
library(tidyverse)

create_distance_between_ltlas <- function(){
# England map
Eng_map <- st_read(dsn = "data/raw/Local_Authority_Districts_(December_2022)_Boundaries_UK_BFC/",
                   layer = "LAD_DEC_2022_UK_BFC")

Eng_map <- subset(Eng_map, startsWith(LAD22CD, "E"))

# If the CRS is not EPSG:4326, transform it
if (st_crs(Eng_map)$epsg != 4326) {
  Eng_map <- st_transform(Eng_map, crs = 4326)
}

# Ensure all geometries are valid
Eng_map <- st_make_valid(Eng_map)

# DISTANCE
# Calculate centroids of each LTLA
centroids <- st_centroid(Eng_map)

# Extract coordinates of the centroids
coords <- st_coordinates(centroids)

# Calculate pairwise distances
distance_matrix <- distm(coords, fun = distHaversine)

# Optionally, assign row and column names to the distance matrix for clarity
rownames(distance_matrix) <- Eng_map$LAD22CD  # assuming LAD22NM is the name column in your data
colnames(distance_matrix) <- Eng_map$LAD22CD


distance_long <- distance_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "location_fine1") %>%
  pivot_longer(cols = -location_fine1, names_to = "location_fine2", 
               values_to = "distance_m", values_drop_na = TRUE)

saveRDS(distance_long, file = "data/processed/distance_between_ltlas.rds")

# return(distance_long)
# NEIGHBOURS
# Identify neighboring LTLAs using st_touches
neighbors_matrix <- st_touches(Eng_map, sparse = FALSE)

# Optionally, assign row and column names to the distance matrix for clarity
rownames(neighbors_matrix) <- Eng_map$LAD22CD  # assuming LAD22NM is the name column in your data
colnames(neighbors_matrix) <- Eng_map$LAD22CD

neighbours_long <- neighbors_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "location_fine1") %>%
  pivot_longer(cols = -location_fine1, names_to = "location_fine2", 
               values_to = "is_neighbour", values_drop_na = TRUE)

saveRDS(neighbours_long, file = "data/processed/is_neighbours_ltlas.rds")


}
