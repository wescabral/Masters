library(tictoc)
tic()

# Parâmetros do Modelo
grid_size <- 50       # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- -0.8       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, 0.01, 0.05)   # Penalização para cada cor
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
  
  energies = rep(0, n_colors)
  
  for(color in 0:(n_colors-1)) {
    # Calcula a energia de cada cor
    energies[color + 1] <- -alpha * count_matches_by_color(grid, i, j, color) - betas[color + 1]
  }
  
  # Transforma as energias em probabilidade
  prob_colors <- exp(energies)/sum(exp(energies))
  
  # Critério de troca
  grid[i, j] <- sample(0:(n_colors-1), 1, prob = prob_colors)
  
}

# Gráfico
image(grid, col = topo.colors(n_colors), main = "Simulação MRF Multi-Cores")


### PSEUDOVEROSSIMILHANÇA

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
      color_obs <- grid[i, j]
      
      # Número de matches
      n_viz_obs <- count_matches_by_color(grid, i, j, color_obs)
      
      # Numerador
      num <- -alpha_est * n_viz_obs - betas_est[color_obs + 1]
      
      # Denominador
      den <- sapply(0:(n_colors-1), function(m) {
        n_viz_m <- count_matches_by_color(grid, i, j, m)
        return(-alpha_est * n_viz_m - betas_est[m + 1])
      })
      
      # Acumulando o logaritmo da probabilidade condicional
      max_den = max(den)
      log_pl_total <- log_pl_total + (num - max(den) - log(sum(exp(den - max(den)))))
    }
  }
  
  # Retorna negativo para que o 'optim' minimize
  return(-log_pl_total)
}

# Otimização
fit_pl <- optim(par = c(0, 0, 0), fn = log_pl)


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
param_samples <- matrix(0, nrow = n_steps, ncol = 3)
current_params <- fit_pl$par
current_log_lik <- -log_pl(current_params) + log_priori(current_params)

for(s in 1:n_steps) {
  # Random Walk
  proposed_params <- current_params + rnorm(3, 0, 1)
  proposed_log_lik <- -log_pl(proposed_params) + log_priori(proposed_params)
  
  # Razão de aceitação
  if(log(runif(1)) < (proposed_log_lik - current_log_lik)) {
    current_params <- proposed_params
    current_log_lik <- proposed_log_lik
  }
  param_samples[s, ] <- current_params
}

### RESULTADOS

results <- data.frame(
  Parâmetro = c("Alpha", "Beta Ref", "Beta 1", "Beta 2"),
  Real = c(alpha, betas),
  PL = c(fit_pl$par[1], "-", fit_pl$par[2], fit_pl$par[3]),
  MH = c(mean(param_samples[,1]), "-", mean(param_samples[,2]), mean(param_samples[,3]))
)

print(results)

toc()