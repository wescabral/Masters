library(R6)
library(Rcpp)
library(ggplot2)

source("./imageData.R")
source("./model.R")
sourceCpp("./sampleImage.cpp")
sourceCpp("./pseudoLikelihood.cpp")
source("./inference.R")
source("./bayesianInference.R")

# 1. Modelo verdadeiro e imagem simulada

seeds_true <- list(
  list(position = c(50, 50),   color = 1, delta = 10),
  list(position = c(100, 100), color = 2, delta = 20)
)

mod <- seedModel$new(alpha = 1.2, beta = 1, seeds = seeds_true, ncolors = 3)
img <- mod$sampleImage(dim = c(150, 150), steps = 80)

# 2. Máxima pseudo-likelihood

seeds_info <- list(
  list(color = 1, n_seeds = 1),
  list(color = 2, n_seeds = 1)
)

fitter <- seedModelFitter$new(image = img, seeds_info = seeds_info, ncolors = 3)
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
  shape_beta  = 2, rate_beta  = 1,
  shape_delta = 10, rate_delta = 1
)

sampler <- seedBayesian$new(image = img, seeds_info = seeds_info, ncolors = 3, priors = prioris)
sampler   # estrutura e hiperparâmetros

# init = NULL: inicializa automaticamente no MLE
res <- sampler$run(n_iter = 5000, step_size = 0.05)

res$plot_logpl_trace()          # diagnóstico de convergência
res$plot_image(n_samples = 300, burn_in = 3000) # imagem + 300 círculos (centros variáveis)

# Burn-in manual
res$plot_traces(burn_in = 3000)
res$plot_densities(burn_in = 3000)
print(res, burn_in = 3000)
res$posterior_mean(burn_in = 3000)
res$credible_interval(level = 0.95, burn_in = 3000)