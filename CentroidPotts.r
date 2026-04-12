library(tictoc)
tic()

# Parâmetros do Modelo
grid_size <- 50      # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- -1.2       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, 0, -1.2)   # Penalização para cada cor
gammas <- c(0, 0, 3) # Dummy da direção da penalização por centroide: negativo penaliza longe
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
  betas_est <- c(0, params[2], params[3])
  gammas_est <- c(0, params[4], params[5])
  
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
      den <- sapply(0:(n_colors-1), function(m) {
        n_viz_m <- count_matches_by_color(grid, i, j, m)
        return(-alpha_est * n_viz_m - betas_est[m + 1] - gammas_est[m + 1] * dist_norm)
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
fit_pl <- optim(par = c(0, 0, 0, 0, 0), fn = log_pl)


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
param_samples <- matrix(0, nrow = n_steps, ncol = 5)
current_params <- fit_pl$par
current_log_lik <- -log_pl(current_params) + log_priori(current_params)

for(s in 1:n_steps) {
  # Random Walk
  proposed_params <- current_params + rnorm(5, 0, 1)
  proposed_log_lik <- -log_pl(proposed_params) + log_priori(proposed_params)
  
  # Razão de aceitação
  if(log(runif(1)) < (proposed_log_lik - current_log_lik)) {
    current_params <- proposed_params
    current_log_lik <- proposed_log_lik
  }
  param_samples[s, ] <- current_params
}

### RESULTADOS

# Tabela
results <- data.frame(
Parâmetro = c("Alpha", "Beta Ref", "Beta 1", "Beta 2", 
              "Gamma Ref", "Gamma 1", "Gamma 2"),
Real = c(alpha, betas, gammas),
PL = c(fit_pl$par[1], "-", fit_pl$par[2], fit_pl$par[3],
       "-", fit_pl$par[4], fit_pl$par[5]),
MH = c(mean(param_samples[,1]), "-", mean(param_samples[,2]), mean(param_samples[,3]),
       "-", mean(param_samples[,4]), mean(param_samples[,5]))
)

print(results)

# Gráficos
plot_list <- list()
par_names <- c("Alpha", "Beta 1", "Beta 2", "Gamma 1", "Gamma 2")
par_values <- c(mean(param_samples[,1]), mean(param_samples[,2]), 
                mean(param_samples[,3]), mean(param_samples[,4]), 
                mean(param_samples[,5]))

for(i in 1:ncol(param_samples)) {
  
  pdf(NULL) 
  dev.control("enable")
  par(mfrow = c(1, 2))

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
       main = paste("Série Temporal:", par_names[i]),
       xlab = "Iteração", ylab = "Valor")
  
  # Valor real
  abline(h = par_values[i], col = "red", lwd = 2)
  
  plot_list[[i]] <- recordPlot()
  dev.off() 
}

par(mfrow = c(1, 1))

toc()
