library(tictoc)
tic()

# Parâmetros do Modelo
grid_size <- 50       # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- -1.3       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, 0.7, 0.2)   # Penalização para cada cor
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
      num <- exp(-alpha_est * n_viz_obs - betas_est[color_obs + 1])
      
      # Denominador
      den <- 0
      for(m in 0:(n_colors-1)) {
        n_viz_m <- count_matches_by_color(grid, i, j, m)
        den <- den + exp(-alpha_est * n_viz_m - betas_est[m + 1])
      }
      
      # Acumulando o logaritmo da probabilidade condicional
      log_pl_total <- log_pl_total + log(num / den)
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

t_obs = t_omega(grid)
psi_alpha <- fit_pl$par[1]
psi_betas <- c(0, fit_pl$par[2], fit_pl$par[3])
n_samples <- 500
t_sims <- matrix(0, nrow = n_samples, ncol = length(t_obs))
temp_grid <- grid

# Simulações de Monte Carlo
for(s in 1:n_samples) {
  for(rep in 1:50000) { # Gibbs para convergir sob psi
    i <- sample(1:grid_size, 1); j <- sample(1:grid_size, 1)
    energies_sim <- sapply(0:(n_colors-1), function(c) {
      -psi_alpha * count_matches_by_color(temp_grid, i, j, c) - psi_betas[c + 1]
    })
    temp_grid[i,j] <- sample(0:(n_colors-1), 1, prob = exp(energies_sim)/sum(exp(energies_sim)))
  }
  t_sims[s, ] <- t_omega(temp_grid)
}

# Log-verossimilhança via Monte Carlo
log_lik_mcml <- function(theta) {
  alpha_est <- theta[1]
  betas_est <- c(0, theta[2], theta[3])
  
  # Numerador
  u_theta_obs <- -alpha_est * t_obs[1] - sum(betas_est * t_obs[2:4])
  
  # Denominador
  diff_alpha <- alpha_est - psi_alpha
  diff_betas <- betas_est - psi_betas
  u_diffs <- -diff_alpha * t_sims[, 1] - (t_sims[, 2:4] %*% diff_betas)
  
  # Log da razão
  max_diff <- max(u_diffs)
  log_ratio_Z <- max_diff + log(mean(exp(u_diffs - max_diff)))
  
  return(u_theta_obs - log_ratio_Z) # Log-verossimilhança aproximada
}

### 4. ALGORITMO METROPOLIS-HASTINGS
n_steps <- 10000
param_samples <- matrix(0, nrow = n_steps, ncol = 3)
current_params <- fit_pl$par
current_log_lik <- log_lik_mcml(current_params)

for(s in 1:n_steps) {
  # Random Walk
  proposed_params <- current_params + rnorm(3, 0, 0.05)
  proposed_log_lik <- log_lik_mcml(proposed_params)
  
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