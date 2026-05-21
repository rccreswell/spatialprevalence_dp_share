# Cluster summary
summary_cluster_assign <- function(results, maxIters, data){
  
  results <- readRDS("data/processed/all_chains_all_clustertrend_assignment_10-2021.rds")
  
  maxIters = 2000
  
  dat_og <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")
  weeks_monthyear <- tibble(weeks = unique(dat_og$week_date) ) %>%
    mutate(monthyear = format(as.Date(weeks), "%Y-%m"))
  temp_monthyear <- "2021-10"
  temp_weeks <- weeks_monthyear %>% 
    filter(monthyear == temp_monthyear) %>%
    select(weeks) %>%
    pull() %>%
    as.Date()
  mmyy <- paste0(month(temp_weeks[1]), "-", year(temp_weeks[1]))
  data <- dat_og %>% 
    filter(week_date %in% temp_weeks) %>% 
    dplyr::select(location_fine, week_date, mean_prev)
  
  # Keep only columns of cluster assignments
  assignments <- results[[1]]
  
  # # Fake assignment matrix to see if code is doing the right thing
  # assignments <- matrix(c(10, 10, 21,
  #                         23, 5, 23), nrow = 3, ncol = 2)
  # maxIters <- 2
  
  # Create an empty matrix
  matrix_size <- nrow(assignments)
  matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  
  # 
  # Iterate over assignments in the second half
  for (col_i in ((maxIters %/% 2) + 1):maxIters) {
  # for (col_i in (maxIters-9):maxIters) {
    assignment <- assignments[, col_i]
    
    # Iterate over characters and indices in the assignment
    for (ii in seq_along(assignment)) {
      for (jj in seq_along(assignment)) {
        # Check if indices are equal
        if (ii == jj) {
          matrix[ii, jj] <- matrix[ii, jj] + 1
        } else if (ii != jj && assignment[ii] == assignment[jj]) {
          matrix[ii, jj] <- matrix[ii, jj] + 1
        }
      }
    }
  }
  
  # Similarity matrix between 0 and 1
  matrix <- matrix/(maxIters-(maxIters %/% 2))
  # matrix <- matrix/10
  
  # number of clusters at each iteration
  len_unique <- function(x) length(unique(x))
  num_clusters <- apply(results[[1]], MARGIN = 2, FUN = len_unique)
  mode_num_clusters <- names(sort(-table(num_clusters)))[1] %>% as.numeric() # MODE!
  mode_num_clusters
  
  matrix_kernel <- as.kernelMatrix(matrix)
  # Perform spectral clustering on this similarity or affinity matrix
  final_clusters <- specc(matrix_kernel, centers = mode_num_clusters, iterations = 10000,
                          nystrom.red = FALSE)
  final_clusters <- as.numeric(final_clusters)
  
  length(unique(final_clusters)) == mode_num_clusters # check final no. of clusters is as expected

  # assign final cluster assignments to prevalence data
  data_cluster <- cbind(data, final_clusters)
  
  return(data_cluster)
}


consensus_clustering <- function(mmyy, maxIters, results_list) {
  n <- nrow(results_list[[1]])
  total_iters <- (maxIters/2)*length(results_list)
  
  # First pass: build probability connectivity matrix
  connectivity <- matrix(0, nrow = n, ncol = n)
  
  for(chain in results_list) {
    assignments <- chain
    for(col_i in ((maxIters %/% 2) + 1):maxIters) {
      assignment <- assignments[, col_i]
      
      for(ii in seq_along(assignment)) {
        for(jj in seq_along(assignment)) {
          if(assignment[ii] == assignment[jj]) {
            connectivity[ii, jj] <- connectivity[ii, jj] + 1
          }
        }
      }
    }
  }
  
  # Probability of co-clustering
  prob_matrix <- connectivity/total_iters
  
  # Second pass: compute SSD for each iteration
  ssd_vec <- c()
  clustering_list <- list()
  
  for (chain in results_list) {
    assignments <- chain
    for (col_i in ((maxIters %/% 2) + 1):maxIters) {
      assignment <- assignments[, col_i]
      
      # Binary connectivity for this iteration
      binary_mat <- outer(assignment, assignment, FUN = function(a, b) as.integer(a == b))
      
      # SSD against probability matrix
      ssd <- sum((binary_mat - prob_matrix)^2)
      
      ssd_vec <- c(ssd_vec, ssd)
      clustering_list[[length(ssd_vec)]] <- assignment
    }
  }
  
  # Pick clustering with minimum SSD
  best_idx <- which.min(ssd_vec)
  best_clustering <- clustering_list[[best_idx]]
  
  return(cluster = best_clustering)
}



