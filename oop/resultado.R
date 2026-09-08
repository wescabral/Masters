#rtools_path <- "C:/rtools45/usr/bin"
#Sys.setenv(PATH = paste(rtools_path, Sys.getenv("PATH"), sep = ";"))
#Sys.which("make")

#rtools_bin <- "C:/rtools45/x86_64-w64-mingw32.static.posix/bin"
#Sys.setenv(PATH = paste(rtools_bin, Sys.getenv("PATH"), sep = ";"))
#Sys.which("g++")

library(R6)
library(Rcpp)
library(ggplot2)
library(mrf2d)

source("./imageData.R")
source("./model.R")
sourceCpp("./sampleImage.cpp")
sourceCpp("./pseudoLikelihood.cpp")
source("./inference.R")
source("./bayesianInference.R")

y <- readRDS("aplicacao.rds")
img <- imageData$new(y, continuous = TRUE)
img$plot()

seeds_info <- list(
  list(color = 1, n_seeds = 3)
)

priors <- list(
  shape_alpha = 2, rate_alpha = 1,
  shape_beta  = 2, rate_beta  = 1,
  shape_delta = 10, rate_delta = 0.5,
  m = 0, tau = 0.001,
  a = 10, b = 10
)

sampler <- seedBayesianCont$new(image = img, seeds_info = seeds_info, ncolors = 2, priors = priors)
sampler$plotZ()   # Z0

set.seed(1)
# init = NULL: inicializa automaticamente no MLE
res <- sampler$run(n_iter = 10000, step_size = 0.01, init = list(alpha = 2, beta = 2, deltas = c(30, 30, 30)))
res$plot_logpl_trace() + 
  labs(y = "Log PL", x = "Iteration", title = NULL)         # diagnóstico de convergência
res$plot_image(n_samples = 300, burn_in = 3000) # imagem + 300 círculos (centros variáveis)

# Burn-in manual
res$plot_traces(burn_in = 3000) + labs(title = NULL)
res$plot_densities(burn_in = 3000) + labs(title = NULL)
print(res, burn_in = 3000)
res$posterior_mean(burn_in = 3000)
res$credible_interval(level = 0.95, burn_in = 3000)
res$plot_positions() + labs(title = NULL)
