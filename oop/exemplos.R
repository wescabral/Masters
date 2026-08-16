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

# 1. Modelo verdadeiro e imagem simulada

seeds_true <- list(
  list(position = c(50, 50),   color = 1, delta = 20),
  list(position = c(100, 100), color = 1, delta = 20)
)

mod <- seedModel$new(alpha = 1.5, beta = 1, seeds = seeds_true, ncolors = 2)
img <- mod$sampleImage(dim = c(150, 150), steps = 80)
img$plot()
Y <- img$matrix + rnorm(prod(dim(img$matrix)), sd = 0.5)
img_cont <- img$setMatrix(Y, continuous = TRUE)
img_cont$plot()

fit <- fit_ghm(Y, mrfi(1), theta = expand_array(0, "onepar", mrfi(1), 1))
img <- imageData$new(fit$Z_pred)

# 2. Máxima pseudo-likelihood

seeds_info <- list(
  list(color = 1, n_seeds = 2)
)

fitter <- seedModelFitter$new(image = img, seeds_info = seeds_info, ncolors = 2)
mle    <- fitter$fit(alpha_init = 1, beta_init = 1)

mle              # estimativas + convergência
mle$estimates    # vetor nomeado: alpha, beta, delta1, delta2
mle$seeds        # seeds kmeans
mle$plot()       # imagem + círculos de influência nos valores MLE

# Simular nova imagem com os parâmetros ajustados
mle$model$sampleImage(c(150, 150), steps = 80)$plot()

# 3. Inferência Bayesiana via MCMC

prioris <- list(
  shape_alpha = 2, rate_alpha = 1, 
  shape_beta  = 5, rate_beta  = 0.5,
  shape_delta = 5, rate_delta = 0.5
)

sampler <- seedBayesian$new(image = img, seeds_info = seeds_info, ncolors = 2, priors = prioris)
sampler   # estrutura e hiperparâmetros

# init = NULL: inicializa automaticamente no MLE
res <- sampler$run(n_iter = 5000, step_size = 0.1)

res$plot_logpl_trace()          # diagnóstico de convergência
res$plot_image(n_samples = 300, burn_in = 3000) # imagem + 300 círculos (centros variáveis)

# Burn-in manual
res$plot_traces(burn_in = 000)
res$plot_densities(burn_in = 3000)
print(res, burn_in = 3000)
res$posterior_mean(burn_in = 3000)
res$credible_interval(level = 0.95, burn_in = 3000)
res$plot_positions()
