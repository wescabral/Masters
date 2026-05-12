# Parâmetros do Modelo
grid_size <- 50      # Tamanho da imagem (50x50)
n_colors <- 2         # Cores
centroides <- matrix(c(10, 15,   # Coordenadas dos centroides
                       35, 25), 
                     ncol = 2, byrow = TRUE)
n_centroides <- nrow(centroides) 
alpha <- -1       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, -1.2)   # Penalização para cada cor
gammas <- matrix(c(0, 5.50, # Dummy da direção da penalização por centroide: negativo penaliza longe
                   0, 3.50),
                 nrow = n_centroides, ncol = n_colors, byrow = TRUE)
max_dist <- sqrt(2 * (grid_size/2)^2)
n_iterations <- 500000 # Sugestão do livro


### GERAÇÃO DA IMAGEM

# Aleatorizando o grid
grid <- matrix(sample(0:(n_colors-1), grid_size^2, replace = TRUE), nrow = grid_size)

# Contagem de vizinhos iguais
count_matches_by_color <- function(g, i, j, color) {
  res <- 0
  if(i > 1) res <- res + (g[i-1, j] == color)
  if(i < grid_size) res <- res + (g[i+1, j] == color)
  if(j > 1) res <- res + (g[i, j-1] == color)
  if(j < grid_size) res <- res + (g[i, j+1] == color)
  return(res)
}

# Troca de cores por reversibilidade
for(step in 1:n_iterations) {
  # Seleciona pixel aleatório
  i <- sample(1:grid_size, 1)
  j <- sample(1:grid_size, 1)
  
  energies <- rep(0, n_colors)
  dists_norm <- (sqrt((i - centroides[,1])^2 + (j - centroides[,2])^2) + 1) / max_dist
  centroide_prx <- which.min(dists_norm)
  dist_prx <- dists_norm[centroide_prx]
  
  for(color in 0:(n_colors-1)) {
    # Calcula a energia de cada cor
    energies[color + 1] <- - alpha * count_matches_by_color(grid, i, j, color) - betas[color + 1] - gammas[centroide_prx, color + 1] * dist_prx
  }
  
  # Transforma as energias em probabilidade
  prob_colors <- exp(energies)/sum(exp(energies))
  
  # Critério de troca
  grid[i, j] <- sample(0:(n_colors-1), 1, prob = prob_colors)
  
}

# Gráfico
image(grid, col = topo.colors(n_colors), main = "")
legend("topright", 
       legend = paste("Cor", 0:(n_colors-1)), 
       fill = topo.colors(n_colors),
       cex = 0.6,           
       pt.cex = 0.8, 
       seg.len = 0,         
       x.intersp = 0.2,     
       y.intersp = 0.4,     
       text.width = 0.05,   
       bg = rgb(1,1,1,0.7), 
       box.lwd = 0.6, 
       inset = 0.01)


### CHUTE DO CENTROIDE

kmeans = kmeans(which(grid == 1, arr.ind = T), centers = n_centroides)
centroides_0 = kmeans$centers


### PSEUDOVEROSSIMILHANÇA

# Matriz de distâncias normalizadas
calc_dists_matrix <- function(centroides) {
  dist_min_matrix <- matrix(0, grid_size, grid_size)
  centroid_idx_matrix <- matrix(0, grid_size, grid_size)
  
  for(i in 1:grid_size) {
    for(j in 1:grid_size) {
      dist_matrix <- (sqrt((i - centroides[,1])^2 + (j - centroides[,2])^2) + 1) / max_dist
      
      dist_min_matrix[i,j] <- min(dist_matrix)
      centroid_idx_matrix[i,j] <- which.min(dist_matrix)
    }
  }
  return(list(dists = dist_min_matrix, indices = centroid_idx_matrix))
}



# Função de Log-PL que opera diretamente sobre a matriz 'grid'
log_pl <- function(params, dists_matrix, centroid_idx_matrix) {
  
  # Recuperando os parâmetros do vetor 'par' de tamanho 3
  alpha_est <- params[1] 
  
  # Fixando o b0 em 0 para identificabilidade
  betas_est <- c(0, params[2])
  gammas_est <- matrix(c(0, params[3],
                         0, params[4]), 
                       nrow = n_centroides, ncol = n_colors, byrow = TRUE)
  
  log_pl_total <- 0
  
  # Percorrendo cada pixel da imagem 50x50
  for(i in 1:grid_size) {
    for(j in 1:grid_size) {
      
      # Cor observada no pixel atual
      color_obs <- grid[i, j]
      
      # Número de matches
      n_viz_obs <- count_matches_by_color(grid, i, j, color_obs)
      
      # Distância normalizada
      dist_norm <- dists_matrix[i,j]
      centroide_idx = centroid_idx_matrix[i,j]
      
      # Numerador
      num <- -alpha_est * n_viz_obs - betas_est[color_obs + 1] - (gammas_est[centroide_idx, color_obs + 1] * dist_norm)
      
      # Denominador
      den <- sapply(0:(n_colors-1), function(m) {
        n_viz_m <- count_matches_by_color(grid, i, j, m)
        return(-alpha_est * n_viz_m - betas_est[m + 1] - (gammas_est[centroide_idx, m + 1] * dist_norm))
      })
      
      # Acumulando o logaritmo da probabilidade condicional
      log_pl_local <- max(den) + log(sum(exp(den - max(den))))
      log_pl_total <- log_pl_total + num - log_pl_local
    }
  }
  
  # Retorna negativo para que o 'optim' minimize
  return(-log_pl_total)
}

# Otimização
fit_pl <- optim(par = c(0, 0, 0, 0), 
                fn = log_pl, 
                dists_matrix = calc_dists_matrix(centroides_0)[[1]],
                centroid_idx_matrix = calc_dists_matrix(centroides_0)[[2]])


### Metropolis-Hastings

# Função para calcular as estatísticas suficientes t(omega)
t_omega <- function(grid) {
  
  # Contagem de matches independente da cor
  matches_hor <- sum(grid[, -ncol(grid)] == grid[, -1])
  matches_ver <- sum(grid[-nrow(grid), ] == grid[-1, ])
  
  # Contagem de cada cor
  color_counts <- as.vector(table(factor(grid, levels = 0:(n_colors-1))))
  return(c(matches = matches_hor + matches_ver, counts = color_counts))
}

# Distribuição log priori
log_priori <- function(params) {
  return(sum(dnorm(params, mean = 0, sd = 10, log = TRUE)))
}

# Algoritmo
n_steps <- 10000
param_samples <- matrix(0, nrow = n_steps, ncol = 4)
centroid_samples <- matrix(0, nrow = n_steps, 4)
current_params <- fit_pl$par
current_centroid <- centroides_0
current_dists <- calc_dists_matrix(current_centroid)[[1]]
current_idx <- calc_dists_matrix(current_centroid)[[2]]
current_log_lik <- -log_pl(current_params, current_dists, current_idx) + log_priori(current_params)

for(s in 1:n_steps) {
  # Random Walk nos parâmetros
  proposed_params <- current_params + rnorm(4, 0, 0.01)
  proposed_log_lik <- -log_pl(proposed_params, current_dists, current_idx) + log_priori(proposed_params)
  
  # Razão de aceitação dos parâmetros
  if(log(runif(1)) < (proposed_log_lik - current_log_lik)) {
    current_params <- proposed_params
    current_log_lik <- proposed_log_lik
  }
  
  # Random Walk no centroide
  proposed_centroid <- current_centroid + round(rnorm(4, 0, 1))
  proposed_centroid <- matrix(pmax(1, pmin(grid_size, proposed_centroid)), nrow = 2)
  proposed_dists <- calc_dists_matrix(proposed_centroid)[[1]]
  proposed_idx <- calc_dists_matrix(proposed_centroid)[[2]]
  proposed_log_lik <- -log_pl(current_params, proposed_dists, proposed_idx) + log_priori(current_params)
  
  # Razão de aceitação dos parâmetros
  if(log(runif(1)) < (proposed_log_lik - current_log_lik)) {
    current_centroid <- proposed_centroid
    current_dists <- proposed_dists
    current_log_lik <- proposed_log_lik
  }
  
  param_samples[s, ] <- current_params
  centroid_samples[s, ] <- current_centroid
  
}

### RESULTADOS

# Tabela
results <- data.frame(
  Parâmetro = c("Alpha", "Beta Ref", "Beta 1", 
                "Gamma1 Ref", "Gamma1 1", 
                "Gamma2 Ref", "Gamma2 1",
                "Centroide1 X", "Centroide2 X",
                "Centroide1 Y", "Centroide2 Y"),
  Real = c(alpha, betas, gammas, centroides),
  PL = c(fit_pl$par[1], "-", fit_pl$par[2], "-", fit_pl$par[3], "-", fit_pl$par[4], 
         centroides_0[1], centroides_0[2], centroides_0[3], centroides_0[4]),
  MH = c(mean(param_samples[1001:10000,1]), "-", mean(param_samples[1001:10000,2]), 
         "-", mean(param_samples[1001:10000,3]),
         "-", mean(param_samples[1001:10000,4]),
         mean(centroid_samples[1001:10000,1]), mean(centroid_samples[1001:10000,2]),
         mean(centroid_samples[1001:10000,3]), mean(centroid_samples[1001:10000,4]))
)

print(results)

# Gráficos
plot_list <- list()
par_names <- c("Alpha", "Beta 1", "Gamma1 1", "Gamma2 1")
par_values <- c(mean(param_samples[1001:10000,1]), mean(param_samples[1001:10000,2]), 
                mean(param_samples[1001:10000,3]), mean(param_samples[1001:10000,4]))

par(mfrow = c(4, 2), mar = c(3, 4, 2, 1), mgp = c(2, 0.7, 0))

for(i in 1:ncol(param_samples)) {
  
  # Amostra do parâmetro i
  sample_current <- param_samples[, i]
  
  # Densidade da distribuição a posteriori
  plot(density(sample_current), breaks = 30, prob = TRUE,
       main = paste("Densidade:", par_names[i]),, 
       ylab = "Densidade",
       col = "black")
  
  # Valor real
  abline(v = par_values[i], col = "red", lwd = 2)
  
  # IC90
  ic_90 <- quantile(sample_current, probs = c(0.05, 0.95))
  abline(v = ic_90[2], col = "blue", lwd = 2, lty = 2)
  abline(v = ic_90[1], col = "blue", lwd = 2, lty = 2)
  
  # "Série temporal"
  ts.plot(sample_current, type = "l", col = "black",
          main = paste("Passeio aleatório:", par_names[i]),
          xlab = "Iteração", ylab = "Valor")
  
  # Valor real
  abline(h = par_values[i], col = "red", lwd = 2)
  
}

par(mfrow = c(1, 1))
