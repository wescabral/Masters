# Parâmetros do Modelo
grid_size <- 50       # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- 0.3        # Negativo para ser valorizar vizinhos iguais
betas <- c(0.5, 0, 0.7)   # Penalização para cada cor
n_iterations <- 500000 # Sugestão do livro

# Aleatorizando o grid
grid <- matrix(sample(0:(n_colors-1), grid_size^2, replace = TRUE), nrow = grid_size)

# Contagem de vizinhos iguais
count_matches <- function(g, i, j, color) {
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
  
  old_color <- grid[i, j]
  new_color <- sample((0:(n_colors-1))[- (old_color + 1)], 1)
  

  # U = -alpha * matches - beta_k
  energy_old <- -alpha * count_matches(grid, i, j, old_color) - betas[old_color + 1]
  energy_new <- -alpha * count_matches(grid, i, j, new_color) - betas[new_color + 1]
  
  delta_u <- energy_new - energy_old
  
  # Critério de troca
  if(delta_u > 0) {
    grid[i, j] <- new_color
  }
}

# Gráfico
image(grid, col = topo.colors(n_colors), main = "Simulação MRF Multi-Cores")


# Pseudoverossimilhança


# Função de Log-PL que opera diretamente sobre a matriz 'grid'
log_pl <- function(params) {
  # Recuperando os parâmetros do vetor 'par' de tamanho 3
  alpha_est <- params[1] 
  
  # Fixando o b0 em 0 para identificabilidade
  b0 <- 0
  b1 <- params[2]
  b2 <- params[3]
  betas_est <- c(b0, b1, b2)
  
  log_pl_total <- 0
  
  # Percorrendo cada pixel da imagem 50x50
  for(i in 1:grid_size) {
    for(j in 1:grid_size) {
      
      # Cor observada no pixel atual
      k_obs <- grid[i, j]
      
      # Número de matches
      n_viz_obs <- count_matches(grid, i, j, k_obs)
      
      # Numerador
      num <- exp(-alpha_est * n_viz_obs - betas_est[k_obs + 1])
      
      # Denominador
      den <- 0
      for(m in 0:(n_colors-1)) {
        n_viz_m <- count_matches(grid, i, j, m)
        den <- den + exp(-alpha_est * n_viz_m - betas_est[m + 1])
      }
      
      # Acumulando o logaritmo da probabilidade condicional
      log_pl_total <- log_pl_total + log(num / den)
    }
  }
  
  # Retornamos o negativo para que o 'optim' minimize
  return(-log_pl_total)
}

# Otimização
fit <- optim(par = c(0.1, 0.1, 0.1), fn = log_pl)

# Exibição dos resultados comparando com os parâmetros reais
cat("--- Parâmetros Reais ---\n")
cat("Alpha:", alpha, "| Betas:", betas, "\n\n")

cat("--- Parâmetros Estimados (Besag via Grid) ---\n")
cat("Alpha Estimado:", fit$par[1], "\n")
cat("Beta_1 Estimado:", fit$par[2], "\n")
cat("Beta_2 Estimado:", fit$par[3], "\n")
