library(R6)
library(Rcpp)

source("./imageData.R")
sourceCpp("./sampleImage.cpp")

seedModel <- R6Class(
  "seedModel",
  public = list(
    initialize = function(alpha, betaVector, seeds, ncolors){
      private$.alpha <- alpha
      private$.betaVector <- betaVector
      private$.seeds <- seeds
      private$.ncolors <- ncolors
      private$.valueSet <- seq_len(ncolors) - 1
    },
    
    sampleImage = function(dim = c(150, 150), steps = 60){
      z <- sample(private$.valueSet, replace = TRUE, size = prod(dim))
      z <- matrix(z, nrow = dim[1], ncol = dim[2])
      
      # Preparar matriz de seeds (n_seeds x 2)
      seedMatrix <- do.call(rbind, lapply(private$.seeds, function(s) s$position))
      
      # Preparar matriz de deltas (n_seeds x n_colors)
      deltaMatrix <- do.call(rbind, lapply(private$.seeds, function(s) s$delta))
      
      z <- gibbsSampler(z, 
                        alpha = private$.alpha,
                        betaVector = private$.betaVector,
                        valueSet = private$.valueSet,
                        seeds = seedMatrix,
                        deltaMatrix = deltaMatrix,
                        steps = steps)
      
      imageData$new(z)
    }
    
  ),
  active = list(
    alpha = function(){
      private$.alpha
    },
    seeds = function(){
      private$.seeds
    }
  ),
  private = list(
    .alpha = NULL,
    .betaVector = NULL,
    .seeds = NULL,
    .ncolors = NULL,
    .valueSet = NULL
  )
)

# Modelo 1

seeds <- list(
  list(position = c(50, 50), delta = c(10, 0)),
  list(position = c(100, 100), delta = c(20, 0))
)

mod <- seedModel$new(
  alpha = 1, 
  betaVector = c(0, 0.9), 
  seeds = seeds,
  ncolors = 2
)

img <- mod$sampleImage(c(150,150), steps = 60)
img$plot()

# Modelo 2 (impossível na parametrização "original")

seeds2 <- list(
  list(position = c(50, 50), delta = c(0, 30, 0)),
  list(position = c(100, 100), delta = c(0, 0, 20))
)

mod2 <- seedModel$new(
  alpha = 1.4, 
  betaVector = c(1, 0, 0), 
  seeds = seeds2,
  ncolors = 3
)

img2 <- mod2$sampleImage(c(150,150), steps = 60)
img2$plot()
