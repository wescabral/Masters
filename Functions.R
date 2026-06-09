library(Rcpp)

generate_image <- function(grid_size,
                           n_colors,
                           params,
                           centroids,
                           n_refreshes) {
  ### Parameters
  
  n_centroids <- nrow(centroids)
  max_dist <- sqrt(2 * (grid_size / 2)^2)          # Maximum distance in the image
  alpha <- params[1]                               # Uniformity bias
  betas <- c(0, params[2:n_colors])                # Color global penalization
  gammas <- matrix(
    # Color local penalization
    data = c(rep(0, n_centroids), params[(n_colors + 1):length(params)]),
    nrow = n_centroids,
    ncol = n_colors,
    byrow = FALSE
  )
  
  
  ### Functions
  
  # Counts uniform neighbors for a pixel
  count_matches_by_color <- function(grid, i, j, color) {
    count <- 0
    if (i > 1)
      count <- count + (grid[i - 1, j] == color)
    if (i < grid_size)
      count <- count + (grid[i + 1, j] == color)
    if (j > 1)
      count <- count + (grid[i, j - 1] == color)
    if (j < grid_size)
      count <- count + (grid[i, j + 1] == color)
    
    return(count)
  }
  
  # Generates a matrix of uniform neighbors count
  compute_all_match_counts <- function(grid) {
    match_counts <- matrix(0, grid_size, grid_size)
    for (i in 1:grid_size) {
      for (j in 1:grid_size) {
        match_counts[i, j] <- count_matches_by_color(grid, i, j, grid[i, j])
      }
    }
    return(match_counts)
  }
  
  # Updates the uniform neighbors count matrix if needed
  update_match_counts <- function(match_counts,
                                  grid,
                                  i,
                                  j,
                                  old_color,
                                  new_color) {
    neighbors <- list(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j +
                                                                          1))
    for (neighbor in neighbors) {
      ni <- neighbor[1]
      nj <- neighbor[2]
      if (ni >= 1 &&
          ni <= grid_size && nj >= 1 && nj <= grid_size) {
        match_counts[ni, nj] <- count_matches_by_color(grid, ni, nj, grid[ni, nj])
      }
    }
    return(match_counts)
  }
  
  
  ### Image generation
  
  # Random grid
  grid <- matrix(sample(0:(n_colors - 1), grid_size^2, replace = TRUE), nrow = grid_size)
  match_counts <- compute_all_match_counts(grid)
  
  for (refresh in 1:n_refreshes) {
    # Pixel order
    pixel_order <- sample(1:(grid_size^2))
    
    for (step in 1:(grid_size^2)) {
      # Refreshing pixel
      pixel_idx <- pixel_order[step]
      i <- ((pixel_idx - 1) %% grid_size) + 1
      j <- ((pixel_idx - 1) %/% grid_size) + 1
      
      old_color <- grid[i, j]
      old_matches <- match_counts[i, j]
      
      # Nearest centroid and the distance to it
      dists_norm <- (sqrt((i - centroids[, 1])^2 + (j - centroids[, 2])^2) + 1) / max_dist
      nrst_centroid <- which.min(dists_norm)
      nrst_dist <- dists_norm[nrst_centroid]
      
      # Energies calculation
      energies <- rep(0, n_colors)
      for (color in 0:(n_colors - 1)) {
        if (color == old_color) {
          matches <- old_matches
        } else {
          matches <- count_matches_by_color(grid, i, j, color)
        }
        energies[color + 1] <- -alpha * matches - betas[color + 1] - gammas[nrst_centroid, color + 1] * nrst_dist
      }
      
      # Gibbs sampler
      prob_colors <- exp(energies) / sum(exp(energies))
      new_color <- sample(0:(n_colors - 1), 1, prob = prob_colors)
      grid[i, j] <- new_color
      
      # Criteria for updating the uniform neighbors count matrix
      if (new_color != old_color) {
        match_counts <- update_match_counts(match_counts, grid, i, j, old_color, new_color)
      }
    }
  }
  
  return(list(
    image = image(grid, col = topo.colors(n_colors), main = ""),
    grid = grid
  ))
}



infer_parameters <- function(grid, n_centroids, n_steps) {
  ### Parameters
  grid_size <- nrow(grid)
  max_dist <- sqrt(2 * (grid_size / 2)^2)
  n_colors <- length(unique(as.vector(grid)))
  
  
  ### Functions
  
  # Calculate distances
  calc_dists_matrix <- function(centroids) {
    i_seq <- 1:grid_size
    j_seq <- 1:grid_size
    distances_list <- lapply(1:nrow(centroids), function(k) {
      dist_matrix <- matrix(NA, grid_size, grid_size)
      for (i in i_seq) {
        for (j in j_seq) {
          dist_matrix[i, j] <- sqrt((i - centroids[k, 1])^2 + (j - centroids[k, 2])^2) + 1
        }
      }
      dist_matrix
    })
    all_dists <- simplify2array(distances_list) / max_dist
    dist_min_matrix <- apply(all_dists, c(1, 2), min)
    centroid_idx_matrix <- apply(all_dists, c(1, 2), which.min)
    return(list(dists = dist_min_matrix, indices = centroid_idx_matrix))
  }
  
  # Calculates log-prior
  log_prior <- function(params) {
    return(sum(dnorm(
      params,
      mean = 0,
      sd = 10,
      log = TRUE
    )))
  }
  
  # C++ function for log pseudo-likelihood
  cppFunction(
    '
    double log_pl(NumericVector params,
                  NumericMatrix dists_matrix,
                  IntegerMatrix centroid_idx_matrix,
                  IntegerMatrix grid_cpp,
                  int grid_size,
                  int n_colors) {

      //Parâmetros
      double alpha_est = params[0];
      NumericVector betas_est = {0.0, params[1]};
      NumericMatrix gammas_est(2, 2);
      gammas_est(0, 0) = 0.0;
      gammas_est(0, 1) = params[2];
      gammas_est(1, 0) = 0.0;
      gammas_est(1, 1) = params[3];

      double log_pl_total = 0.0;


      // Contando matches pra não chamar função global
      for(int i = 0; i < grid_size; i++) {
        for(int j = 0; j < grid_size; j++) {
          int color_obs = grid_cpp(i, j);
          int n_viz_obs = 0;
          if(i > 0 && grid_cpp(i-1, j) == color_obs) n_viz_obs++;
          if(i < grid_size-1 && grid_cpp(i+1, j) == color_obs) n_viz_obs++;
          if(j > 0 && grid_cpp(i, j-1) == color_obs) n_viz_obs++;
          if(j < grid_size-1 && grid_cpp(i, j+1) == color_obs) n_viz_obs++;

          double dist_norm = dists_matrix(i, j);
          int centroid_idx = centroid_idx_matrix(i, j) - 1;
          double num = -alpha_est * n_viz_obs - betas_est[color_obs] -
                       (gammas_est(centroid_idx, color_obs) * dist_norm);

          double max_den = -1e10;
          NumericVector den(n_colors);
          for(int m = 0; m < n_colors; m++) {
            int n_viz_m = 0;
            if(i > 0 && grid_cpp(i-1, j) == m) n_viz_m++;
            if(i < grid_size-1 && grid_cpp(i+1, j) == m) n_viz_m++;
            if(j > 0 && grid_cpp(i, j-1) == m) n_viz_m++;
            if(j < grid_size-1 && grid_cpp(i, j+1) == m) n_viz_m++;

            den[m] = -alpha_est * n_viz_m - betas_est[m] -
                     (gammas_est(centroid_idx, m) * dist_norm);
            if(den[m] > max_den) max_den = den[m];
          }

          double sum_exp = 0.0;
          for(int m = 0; m < n_colors; m++) {
            sum_exp += std::exp(den[m] - max_den);
          }
          double log_pl_local = max_den + std::log(sum_exp);
          log_pl_total += num - log_pl_local;
        }
      }
      return -log_pl_total;
    }'
  )
  
  
  ### First guessing
  
  # Centroid
  kmeans_result <- kmeans(which(grid == 1, arr.ind = TRUE), centers = n_centroids)
  centroids_0 <- kmeans_result$centers
  centroid_order <- order(centroids_0[, 1], centroids_0[, 2])
  centroids_0 <- centroids_0[centroid_order, ]
  
  # Parameters
  dists_info_0 <- calc_dists_matrix(centroids_0)
  fit_pl <- optim(
    par = c(0, 0, 0, 0),
    fn = function(params)
    {
      log_pl(params,
             dists_info_0[[1]],
             dists_info_0[[2]],
             grid,
             grid_size,
             n_colors)
    }
  )
  
  
  ### Metropolis-Hastings
  
  # Current parameters
  param_samples <- matrix(0, nrow = n_steps, ncol = 4)
  centroid_samples <- matrix(0, nrow = n_steps, 4)
  current_params <- fit_pl$par
  current_centroid <- centroids_0
  
  # Current distances
  current_dists_info <- calc_dists_matrix(current_centroid)
  current_dists <- current_dists_info[[1]]
  current_idx <- current_dists_info[[2]]
  
  # Current log-likelihood
  current_log_lik <- -log_pl(current_params,
                             current_dists,
                             current_idx,
                             grid,
                             grid_size,
                             n_colors) + log_prior(current_params)
  
  # Random walk
  for (s in 1:n_steps)
  {
    proposed_params <- current_params + rnorm(4, 0, 0.01)
    proposed_log_lik <- -log_pl(proposed_params,
                                current_dists,
                                current_idx,
                                grid,
                                grid_size,
                                n_colors) + log_prior(proposed_params)
    
    if (log(runif(1)) < (proposed_log_lik - current_log_lik))
    {
      current_params <- proposed_params
      current_log_lik <- proposed_log_lik
    }
    
    proposed_centroid <- current_centroid + round(rnorm(4, 0, 1))
    proposed_centroid <- matrix(pmax(1, pmin(grid_size, proposed_centroid)), nrow = 2)
    proposed_dists_info <- calc_dists_matrix(proposed_centroid)
    proposed_dists <- proposed_dists_info[[1]]
    proposed_idx <- proposed_dists_info[[2]]
    proposed_log_lik <- -log_pl(current_params,
                                proposed_dists,
                                proposed_idx,
                                grid,
                                grid_size,
                                n_colors) + log_prior(current_params)
    
    if (log(runif(1)) < (proposed_log_lik - current_log_lik))
    {
      current_centroid <- proposed_centroid
      current_dists <- proposed_dists
      current_idx <- proposed_idx
      current_log_lik <- proposed_log_lik
    }
    
    param_samples[s, ] <- current_params
    centroid_samples[s, ] <- current_centroid
    
  }
  
  return(
    list(
      parameters = param_samples,
      centroids = centroid_samples,
      parameters_0 = fit_pl$par,
      centroids_0 = centroids_0
    )
  )
}



generate_results <- function(par_names, param_samples, burnin) {
  par_values <- c(
    mean(param_samples[burnin:nrow(param_samples), 1]),
    mean(param_samples[burnin:nrow(param_samples), 2]),
    mean(param_samples[burnin:nrow(param_samples), 3]),
    mean(param_samples[burnin:nrow(param_samples), 4])
  )
  
  plot_obj <- recordPlot()
  
  par(
    mfrow = c(4, 2),
    mar = c(3, 4, 2, 1),
    mgp = c(2, 0.7, 0)
  )
  
  for (i in 1:ncol(param_samples)) {
    # Amostra do parâmetro i
    sample_current <- param_samples[, i]
    
    # Densidade da distribuição a posteriori
    plot(
      density(sample_current),
      breaks = 30,
      prob = TRUE,
      main = paste("Densidade:", par_names[i], "=", par_values[i]),
      ,
      ylab = "Densidade",
      col = "black"
    )
    
    # Valor real
    abline(v = par_values[i],
           col = "red",
           lwd = 2)
    
    # IC90
    ic_90 <- quantile(sample_current, probs = c(0.05, 0.95))
    abline(
      v = ic_90[2],
      col = "blue",
      lwd = 2,
      lty = 2
    )
    abline(
      v = ic_90[1],
      col = "blue",
      lwd = 2,
      lty = 2
    )
    
    # "Série temporal"
    ts.plot(
      sample_current,
      type = "l",
      col = "black",
      main = paste("Passeio aleatório:", par_names[i], "=", par_values[i]),
      xlab = "Iteração",
      ylab = "Valor"
    )
    
    # Valor real
    abline(h = par_values[i],
           col = "red",
           lwd = 2)
    
  }
  
  plot_obj <- recordPlot()
  
  par(mfrow = c(1, 1))
  
  return(list(plot = plot_obj, plot_values = par_values))
  
}

#image = generate_image(50, 2, c(-1, -1, 5, 3), matrix(data = c(10, 10, 35, 35), 2, byrow = TRUE), 200)
#params_infered = infer_parameters(image$grid, n_steps = 100000, n_centroids = 2)
#results = generate_results(c("Alpha", "Beta 1", "Gamma1 1", "Gamma2 1"),
#                           params_infered$parameters,
#                           30000)
