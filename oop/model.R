library(R6)
library(Rcpp)

source("./imageData.R")
sourceCpp("./sampleImage.cpp")
source("./inference.R")

seedModel <- R6Class(
  "seedModel",
  public = list(
    initialize = function(alpha, beta, seeds, ncolors){
      private$.alpha <- alpha
      private$.beta <- beta
      private$.seeds <- seeds
      private$.ncolors <- ncolors
    },
    
    sampleImage = function(dim = c(150, 150), steps = 60){
      z <- matrix(as.integer(sample(0L:(private$.ncolors - 1L), replace = TRUE, size = prod(dim))),
                  nrow = dim[1], ncol = dim[2])
      
      # Preparar matriz de seeds (n_seeds x 2)
      seedMatrix <- do.call(rbind, lapply(private$.seeds, function(s) s$position))
      
      # Preparar vetores de cores e deltas das seeds
      seedColors <- as.integer(sapply(private$.seeds, function(s) s$color))
      seedDeltas <- as.numeric(sapply(private$.seeds, function(s) s$delta))
      
      z <- gibbsSampler(z, 
                        alpha = private$.alpha,
                        beta = private$.beta,
                        ncolors = private$.ncolors,
                        seeds = seedMatrix,
                        seedColors = seedColors,
                        seedDeltas = seedDeltas,
                        steps = steps)
      
      imageData$new(z)
    }
    
  ),
  active = list(
    alpha = function(){
      private$.alpha
    },
    beta = function(){
      private$.beta
    },
    seeds = function(){
      private$.seeds
    }
  ),
  private = list(
    .alpha = NULL,
    .beta = NULL,
    .seeds = NULL,
    .ncolors = NULL
  )
)