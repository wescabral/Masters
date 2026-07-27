library(R6)
sourceCpp("./pseudoLikelihood.cpp")

# classe para o resultado do ajuste
seedFitResult <- R6Class(
  "seedFitResult",
  public = list(
    initialize = function(model, logPL, convergence, estimates, image) {
      private$.model       <- model
      private$.logPL       <- logPL
      private$.convergence <- convergence
      private$.estimates   <- estimates
      private$.image       <- image
    },

    print = function(...) {
      status <- if (private$.convergence == 0) "converged" else paste0("failed (code ", private$.convergence, ")")
      cat(sprintf("<seedFitResult> [%s]\n", status))
      cat(sprintf("  log-PL    : %.4f\n", private$.logPL))
      cat("  estimates :\n")
      for (nm in names(private$.estimates)) {
        cat(sprintf("    %-10s  %.4f\n", nm, private$.estimates[[nm]]))
      }
      invisible(self)
    },

    plot = function(...) {
      beta  <- private$.model$beta
      seeds <- private$.model$seeds
      theta <- seq(0, 2 * pi, length.out = 300)

      p <- private$.image$plot() + coord_fixed()

      for (i in seq_along(seeds)) {
        seed  <- seeds[[i]]
        px    <- seed$position[1]
        py    <- seed$position[2]
        delta <- seed$delta

        p <- p + annotate("point", x = px, y = py,
                          shape = 3, size = 3, color = "white", stroke = 1)

        r <- if (beta > 0) delta / beta - 1 else NA_real_
        if (!is.na(r) && r > 0) {
          p <- p + annotate("path",
                            x = px + r * cos(theta),
                            y = py + r * sin(theta),
                            color = "white", linewidth = 0.7)
        }
      }

      p
    }
  ),
  active = list(
    model       = function() private$.model,
    logPL       = function() private$.logPL,
    convergence = function() private$.convergence,
    estimates   = function() private$.estimates,
    image       = function() private$.image
  ),
  private = list(
    .model       = NULL,
    .logPL       = NULL,
    .convergence = NULL,
    .estimates   = NULL,
    .image       = NULL
  )
)

# Ajuste por maximização da pseudo-likelihood com seeds fixas (posição e cor)
seedModelFitter <- R6Class(
  "seedModelFitter",
  public = list(
    initialize = function(image, seeds, ncolors) {
      private$.image    <- image
      private$.seeds    <- seeds
      private$.ncolors  <- ncolors
      nseeds     <- length(seeds)
      seedMatrix <- do.call(rbind, lapply(seeds, function(s) s$position))
      private$.seedMatrix <- seedMatrix
      private$.seedColors <- as.integer(sapply(seeds, function(s) s$color))

      # Pré-calcular matriz de distâncias: n_pixels x n_seeds
      nrows <- image$dim[1]
      ncols <- image$dim[2]
      n     <- nrows * ncols

      pos_x <- (seq_len(n) - 1L) %% nrows
      pos_y <- (seq_len(n) - 1L) %/% nrows

      distMatrix <- matrix(0.0, nrow = n, ncol = nseeds)
      for (s in seq_len(nseeds)) {
        seed_x <- seedMatrix[s, 1] - 1
        seed_y <- seedMatrix[s, 2] - 1
        distMatrix[, s] <- sqrt((pos_x - seed_x)^2 + (pos_y - seed_y)^2)
      }
      private$.zMatrix    <- matrix(as.integer(image$matrix), nrow = nrows, ncol = ncols)
      private$.distMatrix <- distMatrix
    },

    logPL = function(alpha, beta, deltas) {
      computeLogPL(
        z          = private$.zMatrix,
        alpha      = alpha,
        beta       = beta,
        ncolors    = private$.ncolors,
        distMatrix = private$.distMatrix,
        seedColors = private$.seedColors,
        seedDeltas = as.numeric(deltas)
      )
    },

    # Maximiza a log-PL sobre (alpha > 0, beta, deltas > 0)
    # alpha e deltas são log-reparametrizados internamente para otimização irrestrita
    fit = function(alpha_init = 1.0, beta_init = 0.0, deltas_init = NULL) {
      nseeds <- length(private$.seeds)
      if (is.null(deltas_init)) deltas_init <- rep(1.0, nseeds)

      par_init <- c(log(alpha_init), beta_init, log(deltas_init))

      neg_logPL <- function(par) {
        alpha  <- exp(par[1])
        beta   <- par[2]
        deltas <- exp(par[seq(3, 2 + nseeds)])
        -self$logPL(alpha, beta, deltas)
      }

      result <- optim(par_init, neg_logPL, method = "BFGS")

      alpha_hat  <- exp(result$par[1])
      beta_hat   <- result$par[2]
      deltas_hat <- exp(result$par[seq(3, 2 + nseeds)])

      fitted_seeds <- mapply(
        function(s, d) list(position = s$position, color = s$color, delta = d),
        private$.seeds, deltas_hat,
        SIMPLIFY = FALSE
      )

      fitted_model <- seedModel$new(
        alpha   = alpha_hat,
        beta    = beta_hat,
        seeds   = fitted_seeds,
        ncolors = private$.ncolors
      )

      estimates <- c(
        alpha = alpha_hat,
        beta  = beta_hat,
        setNames(deltas_hat, paste0("delta", seq_along(deltas_hat)))
      )

      seedFitResult$new(
        model       = fitted_model,
        logPL       = -result$value,
        convergence = result$convergence,
        estimates   = estimates,
        image       = private$.image
      )
    },

    print = function(...) {
      nseeds <- length(private$.seeds)
      dim    <- private$.image$dim
      cat("<seedModelFitter>\n")
      cat(sprintf("  image   : %d x %d  (%d colors)\n", dim[1], dim[2], private$.ncolors))
      cat(sprintf("  seeds   : %d\n", nseeds))
      for (i in seq_len(nseeds)) {
        s <- private$.seeds[[i]]
        cat(sprintf("    [%d] position = (%d, %d),  color = %d\n",
                    i, s$position[1], s$position[2], s$color))
      }
      invisible(self)
    }
  ),
  private = list(
    .image      = NULL,
    .seeds      = NULL,
    .ncolors    = NULL,
    .zMatrix    = NULL,
    .seedMatrix = NULL,
    .seedColors = NULL,
    .distMatrix = NULL
  )
)