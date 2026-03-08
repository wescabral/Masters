# Parâmetros do Modelo
grid_size <- 50      # Tamanho da imagem (50x50)
n_colors <- 3         # Cores
alpha <- -0.3       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, 0.5, 0.2)   # Penalização para cada cor
gammas <- c(0, 5, 0) # Dummy da direção da penalização por centroide: negativo penaliza longe
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
    energies[color + 1] <- - alpha * count_matches_by_color(grid, i, j, color) - betas[color + 1] - gammas[color + 1] * (0.5 - dist_norm)
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