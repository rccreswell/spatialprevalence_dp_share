# function to standardize a vector
standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# function to calculate number of unique values in a vector
length_unique <- function(x){
  return(length(unique(x)))
}

# Function to calculate the number of unique values in the last column of each data frame
length_unique_last_column <- function(df) {
  unique_values <- unique(df[, ncol(df)])
  return(length(unique_values))
}

# Function to calculate mode of the number of unique values per column
mode_unique_per_column <- function(df) {
  num_unique <- apply(df, 2, function(col) length(unique(col)))
  mode_val <- as.integer(names(sort(table(num_unique), decreasing = TRUE)[1]))
  return(mode_val)
}

# function to create initial cluster assignmnets for LTLAs based on user 
# specified number of clusters
create_c_init <- function(cl, data){
  cl_vec <- seq(1:cl) 
  c_init <- rep(NA, nrow(data))
  rep_len <- nrow(data)%/%cl
  extra <- nrow(data) - (rep_len*cl)
  c_init_r <- c(rep(cl_vec, rep_len), sample(cl_vec, extra)) 
}

run_multiple_chains <- function(n_more_chains, 
                                data = data, alpha = alpha, mu0 = mu0,
                                sigma0 = sigma0, sigma = sigma, c_init = c_init_r,
                                maxIters = maxIters,
                                results_list){
  # vector of number of clusters
  cl_temp_vec <- sample(1:nrow(data), n_more_chains)
  
  # loop through the different n_more_chains initialisations and save MCMC samples 
  # c_init values in the dp_gibbs function changes for each new initialisation
  for(i in 1:length(cl_temp_vec)){
    cl_temp <- cl_temp_vec[i]
    print(paste0("Initial number of clusters = ", cl_temp))
    
    results_temp <- dp_gibbs(data = data, alpha = alpha, mu0 = mu0,
                             sigma0 = sigma0, sigma = sigma, c_init = create_c_init(cl_temp, data),
                             maxIters = maxIters)
    
    results_list[[i]] <- results_temp
  }
  return(results_list)
}
