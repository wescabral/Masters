tic()

# Parâmetros do Modelo
grid_size <- 50      # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- -1.2       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, -1.2, 0)   # Penalização para cada cor
gammas <- c(0, 3, 0) # Dummy da direção da penalização por centroide: negativo penaliza longe
centroide <- c(25, 25) # Coordenada do centroide
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
  dist_norm <- (sqrt((i - centroide[1])^2 + (j - centroide[2])^2) + 1) / max_dist
  
  for(color in 0:(n_colors-1)) {
    # Calcula a energia de cada cor
    energies[color + 1] <- - alpha * count_matches_by_color(grid, i, j, color) - betas[color + 1] - gammas[color + 1] * dist_norm
  }
  
  # Transforma as energias em probabilidade
  prob_colors <- exp(energies)/sum(exp(energies))
  
  # Critério de troca
  grid[i, j] <- sample(0:(n_colors-1), 1, prob = prob_colors)
  
}

# Gráfico
image(grid, col = topo.colors(n_colors), main = "Simulação MRF Multi-Cores")
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



### PSEUDOVEROSSIMILHANÇA

# Matriz de distâncias normalizadas
dist_matrix <- matrix(0, grid_size, grid_size)
for(i in 1:grid_size) {
  for(j in 1:grid_size) {
    dist_matrix[i,j] <- (sqrt((i - centroide[1])^2 + (j - centroide[2])^2) + 1) / max_dist
  }
}

# Função de Log-PL que opera diretamente sobre a matriz 'grid'
log_pl <- function(params) {
  
  # Recuperando os parâmetros do vetor 'par' de tamanho 3
  alpha_est <- params[1] 
  
  # Fixando o b0 em 0 para identificabilidade
  b0 <- 0
  b1 <- params[2]
  b2 <- params[3]
  g0 <- 0
  g1 <- params[4]
  g2 <- params[5]
  betas_est <- c(b0, b1, b2)
  gammas_est <- c(g0, g1, g2)
  
  log_pl_total <- 0
  
  # Percorrendo cada pixel da imagem 50x50
  for(i in 1:grid_size) {
    for(j in 1:grid_size) {
      
      # Cor observada no pixel atual
      color_obs <- grid[i, j]
      
      # Número de matches
      n_viz_obs <- count_matches_by_color(grid, i, j, color_obs)
      
      # Distância normalizada
      dist_norm <- dist_matrix[i,j]
      
      # Numerador
      num <- -alpha_est * n_viz_obs - betas_est[color_obs + 1] - gammas_est[color_obs + 1] * dist_norm
      
      # Denominador
      den <- 0
      for(m in 0:(n_colors-1)) {
        n_viz_m <- count_matches_by_color(grid, i, j, m)
        den <- den + exp(-alpha_est * n_viz_m - betas_est[m + 1] - gammas_est[m + 1] * dist_norm)
      }
      
      # Acumulando o logaritmo da probabilidade condicional
      log_pl_total <- log_pl_total + num - log(den)
    }
  }
  
  # Retorna negativo para que o 'optim' minimize
  return(-log_pl_total)
}

# Otimização
fit_pl <- optim(par = c(0, 0, 0, 0, 0), fn = log_pl)



### VEROSSIMILHANÇA VIA MONTE CARLO

# Função para calcular as estatísticas suficientes t(omega)
t_omega <- function(grid) {
  
  # Contagem de matches independente da cor
  matches_hor <- sum(grid[, -ncol(grid)] == grid[, -1])
  matches_ver <- sum(grid[-nrow(grid), ] == grid[-1, ])
  
  # Contagem de cada cor
  color_counts <- as.vector(table(factor(grid, levels = 0:(n_colors-1))))
  
  # Contagem de distância por cor
  color_dists <- c(sum(dist_matrix[grid == 0]), sum(dist_matrix[grid == 1]), sum(dist_matrix[grid == 2]))
  
  return(c(matches = matches_hor + matches_ver, counts = color_counts, dists = color_dists))
}

t_obs = t_omega(grid)

# Parâmetros Psi obtidos através da pseudoverossimilhança
psi_alpha <- fit_pl$par[1]
psi_betas <- c(0, fit_pl$par[2], fit_pl$par[3]) # b0 fixo em 0 como referência
psi_gammas <- c(0, fit_pl$par[4], fit_pl$par[5])

n_samples <- 500  # Quantidade de amostras para aproximar a verossimilhança
t_omega_sims <- matrix(0, nrow = n_samples, ncol = 1 + n_colors*2)
temp_grid <- grid

for(s in 1:n_samples) {
  
  # Realiza algumas iterações de Gibbs para cada amostra
  for(rep in 1:50000) {
    
    i <- sample(1:grid_size, 1) 
    j <- sample(1:grid_size, 1)
    
    # Distância normalizada
    dist_norm <- dist_matrix[i,j]
    
    energies_sim <- sapply(0:(n_colors-1), function(c) {
      -psi_alpha * count_matches_by_color(temp_grid, i, j, c) - psi_betas[c + 1] - psi_gammas[c + 1] * dist_norm
    })
    
    probs <- exp(energies_sim) / sum(exp(energies_sim))
    
    temp_grid[i,j] <- sample(0:(n_colors-1), 1, prob = probs)
  }
  
  t_omega_sims[s, ] <- t_omega(temp_grid)
}

# Função de log-verossimilhança via Monte Carlo
log_mcml <- function(params) {
  alpha_est <- params[1]
  betas_est <- c(0, params[2], params[3])
  gammas_est <- c(0, params[4], params[5])
  
  # Energia da imagem observada sob o novo parâmetro theta
  u_theta_obs <- -alpha_est * t_obs[1] - sum(betas_est * t_obs[2:4]) - sum(gammas_est * t_obs[5:7])
  
  # Diferença de energia para as amostras simuladas
  diff_alpha <- alpha_est - psi_alpha
  diff_betas <- betas_est - psi_betas
  diff_gammas <- gammas_est - psi_gammas
  
  # Diferença u_diff para cada amostra simulada
  u_diffs <- -diff_alpha * t_omega_sims[, 1] - (t_omega_sims[, 2:4] %*% diff_betas) - (t_omega_sims[, 5:7] %*% diff_gammas)
  max_diff <- max(u_diffs)
  
  # Log-verossimilhança
  log_ratio_Z <- max_diff + log(mean(exp(u_diffs - max_diff)))
  l_theta <- u_theta_obs - log_ratio_Z
  
  return(-l_theta) # Negativo para minimizar
}

# Otimização final MCML
fit_mcml <- optim(par = fit_pl$par, fn = log_mcml)



### RESULTADOS

results <- data.frame(
  Parâmetro = c("Alpha", "Beta Ref", "Beta 1", "Beta 2", "Gamma Ref", "Gamma 1", "Gamma 2"),
  Real = c(alpha, betas, gammas),
  PL = c(fit_pl$par[1], "-", fit_pl$par[2], fit_pl$par[3], "-", fit_pl$par[4], fit_pl$par[5]),
  MCML = c(fit_mcml$par[1], "-", fit_mcml$par[2], fit_mcml$par[3], "-", fit_mcml$par[4], fit_mcml$par[5])
)
print(results)

toc()