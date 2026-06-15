library(Rcpp)
library(reshape2)
library(ggplot2)
library(ggforce)

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
  centroids_0 <- round(kmeans_result$centers, 0)

  # Sort by distance from origin
  distances_from_origin <- rowSums(centroids_0^2)
  centroid_order <- order(distances_from_origin)
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

  # Labeling parameters
  param_names <- c("alpha")
  
  for (i in 1:(n_colors - 1)) {
    param_names <- c(param_names, paste0("beta_", i))
  }
  
  for (j in 1:n_colors) {
    for (i in 1:n_centroids) {
      param_names <- c(param_names, paste0("gamma_", i, "_", j))
    }
  }

  param_samples_df <- as.data.frame(param_samples)
  colnames(param_samples_df) <- param_names[1:ncol(param_samples_df)]


  ### Centroids heatmap

  coords_x <- centroid_samples[, 1:n_centroids]
  coords_y <- centroid_samples[, (n_centroids + 1):ncol(centroid_samples)]

  heatmap_matrix <- matrix(0, nrow = grid_size, ncol = grid_size)

  for (i in 1:nrow(centroid_samples)) {
  for (j in 1:n_centroids) {
  x_idx <- coords_x[i, j] 
  y_idx <- coords_y[i, j] 

  # Validating indexes
  if (x_idx > 0 && x_idx <= ncol(heatmap_matrix) &&
  y_idx > 0 && y_idx <= nrow(heatmap_matrix)) {
  heatmap_matrix[y_idx, x_idx] <- heatmap_matrix[y_idx, x_idx] + 1
  }
  }
  }

  
  return(
    list(
      parameters = param_samples_df,
      centroids = centroid_samples,
      parameters_0 = fit_pl$par,
      centroids_0 = centroids_0,
      centroids_heatmap = heatmap_matrix
    )
  )
}



generate_results <- function(par_names, param_samples, centroid_samples, heatmap_matrix, burnin) {
  ### Parameters plot

  par_values <- c(
    mean(param_samples[burnin:nrow(param_samples), 1]),
    mean(param_samples[burnin:nrow(param_samples), 2]),
    mean(param_samples[burnin:nrow(param_samples), 3]),
    mean(param_samples[burnin:nrow(param_samples), 4])
  )
  
  plot_params <- recordPlot()
  
  par(
    mfrow = c(4, 2),
    mar = c(3, 4, 2, 1),
    mgp = c(2, 0.7, 0)
  )
  
  for (i in 1:ncol(param_samples)) {
    sample_current <- param_samples[, i]
    
    # Posterior density
    plot(
      density(sample_current),
      breaks = 30,
      prob = TRUE,
      main = paste("Densidade:", par_names[i], "=", par_values[i]),
      ,
      ylab = "Densidade",
      col = "black"
    )
    
    # Posterior average
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
    
    # "Time series"
    ts.plot(
      sample_current,
      type = "l",
      col = "black",
      main = paste("Passeio aleatório:", par_names[i], "=", par_values[i]),
      xlab = "Iteração",
      ylab = "Valor"
    )
    
    # Posterior average
    abline(h = par_values[i],
           col = "red",
           lwd = 2)
    
  }
  
  plot_params <- recordPlot()
    
  par(mfrow = c(1, 1))

  
  ### Centroid influence area plot

  calculate_beta_gamma_ratios <- function(param_samples, burnin) {
  
    # Extract means after burnin
    param_means <- colMeans(param_samples[burnin:nrow(param_samples), ])
    
    # Extract beta and gamma columns using names
    col_names <- names(param_means)
    
    # Identify beta and gamma columns
    beta_cols <- col_names[grepl("^beta_", col_names)]
    gamma_cols <- col_names[grepl("^gamma_", col_names)]
    
    # Create results list
    ratios_list <- list()
    
    # For each gamma column, find its corresponding beta
    for (gamma_col in gamma_cols) {
      # Parse gamma_i_j to get color j
      gamma_parts <- strsplit(gamma_col, "_")[[1]]
      centroid_idx <- as.numeric(gamma_parts[2])
      color_idx <- as.numeric(gamma_parts[3])
      
      # Find corresponding beta (beta_j where j is the color)
      # Note: beta_1 corresponds to color 1, beta_2 to color 2, etc.
      beta_col <- paste0("beta_", color_idx)
      
      if (beta_col %in% names(param_means)) {
        beta_val <- param_means[beta_col]
        gamma_val <- param_means[gamma_col]
        
        # Calculate ratio
        if (gamma_val != 0) {
          ratio <- beta_val / gamma_val
        } else {
          ratio <- NA
        }
        
        ratios_list[[gamma_col]] <- data.frame(
          Column = gamma_col,
          Centroid = centroid_idx,
          Color = color_idx,
          Beta = beta_val,
          Gamma = gamma_val,
          Ratio_Beta_Gamma = ratio,
          row.names = NULL
        )
      }
    }
    
    # Combine all into one data frame
    ratios_df <- do.call(rbind, ratios_list)    
    return(ratios_df)
  }

  # Circle creating function
  create_circle <- function(idx) {
    center_x <- infarea_df$x[idx]
    center_y <- infarea_df$y[idx]
    radius <- infarea_df$radius[idx]
    
    angles <- seq(0, 2*pi, length.out = 100)
    data.frame(
      x = center_x + radius * cos(angles),
      y = center_y + radius * sin(angles),
      centroid_id = as.character(idx),
      stringsAsFactors = FALSE
    )
  }

  infarea_df <- data.frame(
    x = round(colMeans(centroid_samples[burnin:nrow(centroid_samples), 1:(ncol(centroid_samples)/2)])),
    y = round(colMeans(centroid_samples[burnin:nrow(centroid_samples), (1 + ncol(centroid_samples)/2):ncol(centroid_samples)])),
    radius = calculate_beta_gamma_ratios(param_samples, burnin)
  )

  circles_data_list <- lapply(1:nrow(infarea_df), create_circle)
  circles_data <- do.call(rbind, circles_data_list)
  rownames(circles_data) <- NULL

  plot_infarea <- ggplot() +
    geom_path(data = circles_data, aes(x = x, y = y, group = centroid_id),
              color = "blue", linewidth = 1) +
    geom_point(data = infarea_df, aes(x = x, y = y),
               color = "red", size = 3) +
    xlim(0, 50) +
    ylim(0, 50) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Grid 50x50 com Círculos de Influência")
  
  

  ### Centroid heatmap plot

  heatmap_df <- melt(heatmap_matrix)
  colnames(heatmap_df) <- c("Y", "X", "Frequency")
      
  heatmap <- ggplot(heatmap_df, aes(x = X, y = Y, fill = Frequency)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "darkred") +
        labs(title = paste("Centroids heatmap"),
             x = "X",
             y = "Y",
             fill = "Frequency") +
        theme_minimal() +
        theme(panel.grid = element_blank())

  
  return(list(plot_params = plot_params, 
    plot_values = par_values,
    plot_infarea = plot_infarea, 
    plot_heatmap = heatmap))
}

image = generate_image(50, 2, c(-1, -1, 5, 3), matrix(data = c(10, 10, 35, 35), 2, byrow = TRUE), 200)
params_infered = infer_parameters(image$grid, n_steps = 100000, n_centroids = 2)
results = generate_results(c("Alpha", "Beta 1", "Gamma1 1", "Gamma2 1"),
                           params_infered$parameters,
                           params_infered$centroids,
                           params_infered$centroids_heatmap,
                           30000)