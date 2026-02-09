# Parâmetros do Modelo
grid_size <- 50       # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- 0.7        # Negativo para ser valorizar vizinhos iguais
betas <- c(0, 0, 0.7)   # Penalização para cada cor
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

