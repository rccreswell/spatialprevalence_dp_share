# create_prob_clustering <- function(mmyy, maxIters){
#   results <- readRDS(paste0("data/processed/all_clustertrend_assignment_", mmyy, ".rds"))
#   # Cluster summary
#   # Keep only columns of cluster assignments
#   assignments <- results %>% select(-1)
# 
#   # Create an empty matrix
#   matrix_size <- nrow(assignments)
#   matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
#   
#   # 
#   # Iterate over assignments in the second half
#   for (col_i in ((maxIters %/% 2) + 1):maxIters) {
#     assignment <- assignments[, col_i]
#     
#     # Iterate over characters and indices in the assignment
#     for (ii in seq_along(assignment)) {
#       for (jj in seq_along(assignment)) {
#         # Check if indices are equal
#         if (ii == jj) {
#           matrix[ii, jj] <- matrix[ii, jj] + 1
#         } else if (ii != jj && assignment[ii] == assignment[jj]) {
#           matrix[ii, jj] <- matrix[ii, jj] + 1
#         }
#       }
#     }
#   }
#   
#   # Similarity matrix between 0 and 1
#   matrix <- matrix/(maxIters-(maxIters %/% 2))
#   
#   # Keep only the upper triangle (excluding the diagonal)
#   triangle_matrix <- matrix
#   triangle_matrix[lower.tri(triangle_matrix, diag = TRUE)] <- NA
#   
#   colnames(triangle_matrix) <- results[ ,1]
#   rownames(triangle_matrix) <- results[ ,1]
#   
#   df_long <- triangle_matrix %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "location_fine1") %>%
#     pivot_longer(cols = -location_fine1, names_to = "location_fine2", 
#                  values_to = "prob", values_drop_na = TRUE)
#   
#   saveRDS(df_long, file = paste0("data/processed/prob_clustertrend_assignment_",
#                                  mmyy, ".rds"))
#   
# }

create_prob_clustering <- function(mmyy, maxIters, data_location_fine) {
  # Read in all chains (results is now a list of chains)
  results_list <- readRDS(paste0("data/processed/all_chains_all_clustertrend_assignment_", mmyy, ".rds"))
  
  # Number of data points from first chain
  n <- nrow(results_list[[1]])
  
  # Initialize global co-clustering matrix
  co_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Loop over chains
  for (chain in results_list) {
    # Drop first column if it contains IDs
    assignments <- chain
    
    # Loop over second half of iterations (burn-in removed)
    for (col_i in ((maxIters %/% 2) + 1):maxIters) {
      assignment <- assignments[, col_i]
      
      # Update co-clustering counts
      for (ii in seq_along(assignment)) {
        for (jj in seq_along(assignment)) {
          if (ii == jj || assignment[ii] == assignment[jj]) {
            co_matrix[ii, jj] <- co_matrix[ii, jj] + 1
          }
        }
      }
    }
  }
  
  # Normalize to get similarity matrix (values between 0 and 1)
  total_iterations <- length(results_list) * (maxIters - (maxIters %/% 2))
  similarity_matrix <- co_matrix / total_iterations
  
  # Keep only upper triangle (excluding diagonal)
  triangle_matrix <- similarity_matrix
  triangle_matrix[lower.tri(triangle_matrix, diag = TRUE)] <- NA
  
  # Set row/column names from first chain (assuming first column was IDs)
  colnames(triangle_matrix) <- data_location_fine$location_fine
  rownames(triangle_matrix) <- data_location_fine$location_fine
  
  # Convert to long format for saving
  df_long <- triangle_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "location_fine1") %>%
    pivot_longer(
      cols = -location_fine1,
      names_to = "location_fine2",
      values_to = "prob",
      values_drop_na = TRUE
    )
  
  # Save as RDS
  saveRDS(df_long, file = paste0("data/processed/prob_clustertrend_assignment_", mmyy, ".rds"))
}

