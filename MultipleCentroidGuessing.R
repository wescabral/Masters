library(tictoc)
tic()

# Parâmetros do Modelo
grid_size <- 50
n_colors <- 2
centroides <- matrix(c(10, 15, 35, 25), ncol = 2, byrow = TRUE)
n_centroides <- nrow(centroides) 
alpha <- -1
betas <- c(0, -1.2)
gammas <- matrix(c(0, 5.50, 0, 3.50), nrow = n_centroides, ncol = n_colors, byrow = TRUE)
max_dist <- sqrt(2 * (grid_size/2)^2)
n_iterations <- 500000

### GERAÇÃO DA IMAGEM
grid <- matrix(sample(0:(n_colors-1), grid_size^2, replace = TRUE), nrow = grid_size)

count_matches_by_color <- function(g, i, j, color) {
  res <- 0
  if(i > 1) res <- res + (g[i-1, j] == color)
  if(i < grid_size) res <- res + (g[i+1, j] == color)
  if(j > 1) res <- res + (g[i, j-1] == color)
  if(j < grid_size) res <- res + (g[i, j+1] == color)
  return(res)
}

compute_all_match_counts <- function(g) {
  match_counts <- matrix(0, grid_size, grid_size)
  for(i in 1:grid_size) {
    for(j in 1:grid_size) {
      match_counts[i, j] <- count_matches_by_color(g, i, j, g[i, j])
    }
  }
  return(match_counts)
}

match_counts <- compute_all_match_counts(grid)

update_match_counts <- function(match_counts, grid, i, j, old_color, new_color) {
  neighbors <- list(c(i, j), c(i-1, j), c(i+1, j), c(i, j-1), c(i, j+1))
  for(neighbor in neighbors) {
    ni <- neighbor[1]
    nj <- neighbor[2]
    if(ni >= 1 && ni <= grid_size && nj >= 1 && nj <= grid_size) {
      match_counts[ni, nj] <- count_matches_by_color(grid, ni, nj, grid[ni, nj])
    }
  }
  return(match_counts)
}

for(step in 1:n_iterations) {
  i <- sample(1:grid_size, 1)
  j <- sample(1:grid_size, 1)
  old_color <- grid[i, j]
  old_matches <- match_counts[i, j]
  energies <- rep(0, n_colors)
  dists_norm <- (sqrt((i - centroides[,1])^2 + (j - centroides[,2])^2) + 1) / max_dist
  centroide_prx <- which.min(dists_norm)
  dist_prx <- dists_norm[centroide_prx]
  for(color in 0:(n_colors-1)) {
    if(color == old_color) {
      matches <- old_matches
    } else {
      matches <- count_matches_by_color(grid, i, j, color)
    }
    energies[color + 1] <- - alpha * matches - betas[color + 1] - gammas[centroide_prx, color + 1] * dist_prx
  }
  prob_colors <- exp(energies)/sum(exp(energies))
  new_color <- sample(0:(n_colors-1), 1, prob = prob_colors)
  grid[i, j] <- new_color
  if(new_color != old_color) {
    match_counts <- update_match_counts(match_counts, grid, i, j, old_color, new_color)
  }
}

image(grid, col = topo.colors(n_colors), main = "")

### CHUTE DO CENTROIDE
kmeans_result <- kmeans(which(grid == 1, arr.ind = TRUE), centers = n_centroides)
centroides_0 <- kmeans_result$centers

### PSEUDOVEROSSIMILHANÇA
calc_dists_matrix <- function(centroides) {
  i_seq <- 1:grid_size
  j_seq <- 1:grid_size
  distances_list <- lapply(1:nrow(centroides), function(k) {
    dist_matrix <- matrix(NA, grid_size, grid_size)
    for(i in i_seq) {
      for(j in j_seq) {
        dist_matrix[i, j] <- sqrt((i - centroides[k, 1])^2 + (j - centroides[k, 2])^2) + 1
      }
    }
    dist_matrix
  })
  all_dists <- simplify2array(distances_list) / max_dist
  dist_min_matrix <- apply(all_dists, c(1, 2), min)
  centroid_idx_matrix <- apply(all_dists, c(1, 2), which.min)
  return(list(dists = dist_min_matrix, indices = centroid_idx_matrix))
}

# C++ function - grid is now a parameter
cppFunction('
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
      int centroide_idx = centroid_idx_matrix(i, j) - 1;
      double num = -alpha_est * n_viz_obs - betas_est[color_obs] - 
                   (gammas_est(centroide_idx, color_obs) * dist_norm);
      
      double max_den = -1e10;
      NumericVector den(n_colors);
      for(int m = 0; m < n_colors; m++) {
        int n_viz_m = 0;
        if(i > 0 && grid_cpp(i-1, j) == m) n_viz_m++;
        if(i < grid_size-1 && grid_cpp(i+1, j) == m) n_viz_m++;
        if(j > 0 && grid_cpp(i, j-1) == m) n_viz_m++;
        if(j < grid_size-1 && grid_cpp(i, j+1) == m) n_viz_m++;
        
        den[m] = -alpha_est * n_viz_m - betas_est[m] - 
                 (gammas_est(centroide_idx, m) * dist_norm);
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
}
')

# Otimização
dists_info_0 <- calc_dists_matrix(centroides_0)
fit_pl <- optim(par = c(0, 0, 0, 0), 
                fn = function(params) {
                  log_pl(params, dists_info_0[[1]], dists_info_0[[2]], grid, grid_size, n_colors)
                })

### Metropolis-Hastings
log_priori <- function(params) {
  return(sum(dnorm(params, mean = 0, sd = 10, log = TRUE)))
}

n_steps <- 1000  # Reduced for profvis
param_samples <- matrix(0, nrow = n_steps, ncol = 4)
centroid_samples <- matrix(0, nrow = n_steps, 4)
current_params <- fit_pl$par
current_centroid <- centroides_0

current_dists_info <- calc_dists_matrix(current_centroid)
current_dists <- current_dists_info[[1]]
current_idx <- current_dists_info[[2]]

current_log_lik <- -log_pl(current_params, current_dists, current_idx, grid, grid_size, n_colors) + log_priori(current_params)

for(s in 1:n_steps) {
  proposed_params <- current_params + rnorm(4, 0, 0.01)
  proposed_log_lik <- -log_pl(proposed_params, current_dists, current_idx, grid, grid_size, n_colors) + log_priori(proposed_params)
  
  if(log(runif(1)) < (proposed_log_lik - current_log_lik)) {
    current_params <- proposed_params
    current_log_lik <- proposed_log_lik
  }
  
  proposed_centroid <- current_centroid + round(rnorm(4, 0, 1))
  proposed_centroid <- matrix(pmax(1, pmin(grid_size, proposed_centroid)), nrow = 2)
  proposed_dists_info <- calc_dists_matrix(proposed_centroid)
  proposed_dists <- proposed_dists_info[[1]]
  proposed_idx <- proposed_dists_info[[2]]
  proposed_log_lik <- -log_pl(current_params, proposed_dists, proposed_idx, grid, grid_size, n_colors) + log_priori(current_params)
  
  if(log(runif(1)) < (proposed_log_lik - current_log_lik)) {
    current_centroid <- proposed_centroid
    current_dists <- proposed_dists
    current_idx <- proposed_idx
    current_log_lik <- proposed_log_lik
  }
  
  param_samples[s, ] <- current_params
  centroid_samples[s, ] <- current_centroid
}

results <- data.frame(
  Parâmetro = c("Alpha", "Beta Ref", "Beta 1", "Gamma1 Ref", "Gamma2 Ref", "Gamma1 1", "Gamma2 1"),
  Real = c(alpha, betas, gammas),
  PL = c(fit_pl$par[1], "-", fit_pl$par[2], "-", "-", fit_pl$par[3], fit_pl$par[4]),
  MH = c(mean(param_samples[,1]), "-", mean(param_samples[,2]), "-", "-", mean(param_samples[,3]), mean(param_samples[,4]))
)
print(results)

toc()
