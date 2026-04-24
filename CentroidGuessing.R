
# Parâmetros do Modelo
grid_size <- 50      # Tamanho da imagem (50x50)
n_colors <- 2         # Cores
alpha <- -1.2       # Negativo para ser valorizar vizinhos iguais
betas <- c(0, -1.2)   # Penalização para cada cor
gammas <- c(0, 3) # Dummy da direção da penalização por centroide: negativo penaliza longe
centroide <- c(42, 10) # Coordenada do centroide
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


### Chute do centroide
c0_mean <- matrix(nrow = n_colors, ncol = 2)
c0_sd <- matrix(nrow = n_colors, ncol = 2)

for (c in 0:n_colors-1) {
  hor_matches = which(grid == c) %% 50
  ver_matches = which(t(grid) == c) %% 50
  c0_mean[c+1, 1] = mean(replace(hor_matches, which(hor_matches == c), 50))
  c0_mean[c+1, 2] = mean(replace(ver_matches, which(ver_matches == c), 50))
  c0_sd[c+1, 1] = sd(replace(hor_matches, which(hor_matches == c), 50))
  c0_sd[c+1, 2] = sd(replace(ver_matches, which(ver_matches == c), 50))
}

print(paste0("Centroide relacionado a cor ", which.min(rowMeans(c0_sd))-1, 
             " é (", round(c0_mean[which.min(rowMeans(c0_sd)),1]), 
             ", ", round(c0_mean[which.min(rowMeans(c0_sd)),2]), ")"))
