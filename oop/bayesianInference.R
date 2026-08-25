sourceCpp("./pseudoLikelihood.cpp")

# Resultado da inferência Bayesiana via MH
seedBayesianResult <- R6Class(
  "seedBayesianResult",
  public = list(
    initialize = function(chain, seeds, image) {
      private$.chain  <- chain
      private$.seeds  <- seeds
      private$.image  <- image
      nseeds <- length(seeds)
      private$.param_names <- c("alpha", "beta", paste0("delta", seq_len(nseeds)))
      private$.pos_names   <- c(rbind(paste0("x_", seq_len(nseeds)),
                                      paste0("y_", seq_len(nseeds))))
      if (nrow(chain) > 1L) {
        par_diff <- diff(chain[, private$.param_names, drop = FALSE])
        pos_diff <- diff(chain[, private$.pos_names,   drop = FALSE])
        private$.accept_rate_param <- mean(rowSums(par_diff != 0) > 0)
        private$.accept_rate_pos   <- mean(rowSums(pos_diff != 0) > 0)
      } else {
        private$.accept_rate_param <- NA_real_
        private$.accept_rate_pos   <- NA_real_
      }
    },

    print = function(burn_in = NULL, ...) {
      n <- nrow(private$.chain)
      if (is.null(burn_in)) burn_in <- floor(0.2 * n)
      cat("<seedBayesianResult>\n")
      cat(sprintf("  chain               : %d iterations  (summarising from %d)\n", n, burn_in + 1L))
      cat(sprintf("  acceptance (params) : %.1f%%\n", 100 * private$.accept_rate_param))
      cat(sprintf("  acceptance (pos)    : %.1f%%\n", 100 * private$.accept_rate_pos))
      cat("  posterior (mean  [95% CI]):\n")
      ci    <- self$credible_interval(level = 0.95, burn_in = burn_in)
      means <- self$posterior_mean(burn_in = burn_in)
      for (nm in private$.param_names) {
        row <- ci[ci$param == nm, ]
        cat(sprintf("    %-10s  %.4f  [%.4f, %.4f]\n",
                    nm, means[[nm]], row$lower, row$upper))
      }
      invisible(self)
    },

    posterior_mean = function(burn_in = NULL) {
      colMeans(private$get_param_chain(burn_in))
    },

    credible_interval = function(level = 0.95, burn_in = NULL) {
      a      <- (1 - level) / 2
      bounds <- apply(private$get_param_chain(burn_in), 2, quantile, probs = c(a, 1 - a))
      data.frame(
        param = private$.param_names,
        lower = bounds[1, ],
        upper = bounds[2, ],
        row.names = NULL
      )
    },

    plot_traces = function(burn_in = NULL) {
      n <- nrow(private$.chain)
      if (is.null(burn_in)) burn_in <- floor(0.2 * n)
      long_df <- private$chain_long(burn_in)
      ggplot(long_df, aes(x = iteration, y = value)) +
        geom_line(linewidth = 0.3, alpha = 0.8) +
        facet_wrap(~param, scales = "free_y", ncol = 1) +
        labs(x = "Iteration", y = NULL, title = "MCMC Traces") +
        theme_bw(base_size = 11) +
        theme(strip.text = element_text(face = "bold"))
    },

    plot_densities = function(burn_in = NULL) {
      n <- nrow(private$.chain)
      if (is.null(burn_in)) burn_in <- floor(0.2 * n)
      long_df <- private$chain_long(burn_in)
      mean_df <- data.frame(
        param    = private$.param_names,
        mean_val = as.numeric(self$posterior_mean(burn_in))
      )
      mean_df$param <- factor(mean_df$param, levels = private$.param_names)
      ggplot(long_df, aes(x = value)) +
        geom_density(fill = "steelblue", alpha = 0.4) +
        geom_vline(data = mean_df, aes(xintercept = mean_val),
                   linetype = "dashed", color = "firebrick") +
        facet_wrap(~param, scales = "free", ncol = 2) +
        labs(x = NULL, y = "Density", title = "Posterior Distributions") +
        theme_bw(base_size = 11) +
        theme(strip.text = element_text(face = "bold"))
    },

    plot_radii = function(burn_in = NULL) {
      chain  <- private$get_param_chain(burn_in)
      nseeds <- length(private$.seeds)
      radii_list <- lapply(seq_len(nseeds), function(s) {
        r <- chain[, paste0("delta", s)] / chain[, "beta"] - 1
        data.frame(seed = paste0("seed ", s), radius = r)
      })
      radii_df       <- do.call(rbind, radii_list)
      radii_df$seed  <- factor(radii_df$seed)
      valid_df       <- radii_df[radii_df$radius > 0, ]

      if (nrow(valid_df) == 0) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No positive radii in posterior samples") +
          theme_void())
      }

      frac_valid <- tapply(radii_df$radius > 0, radii_df$seed, mean)
      caption <- paste(
        sapply(names(frac_valid), function(s)
          sprintf("%s: %.0f%% positive", s, 100 * frac_valid[[s]])),
        collapse = "  |  "
      )

      ggplot(valid_df, aes(x = radius)) +
        geom_density(fill = "darkorange", alpha = 0.4) +
        geom_vline(xintercept = 0, linetype = "dotted") +
        facet_wrap(~seed, scales = "free") +
        labs(x = "\u03b4/\u03b2 - 1", y = "Density",
             title = "Posterior: influence radius per seed",
             caption = caption) +
        theme_bw(base_size = 11) +
        theme(strip.text = element_text(face = "bold"))
    },

    plot_image = function(n_samples = 200, burn_in = NULL) {
      chain  <- private$get_chain(burn_in)
      n_kept <- nrow(chain)
      idx    <- sample(seq_len(n_kept), min(n_samples, n_kept))
      sub    <- chain[idx, , drop = FALSE]
      nseeds <- length(private$.seeds)
      theta  <- seq(0, 2 * pi, length.out = 100)

      circle_rows <- list()
      for (i in seq_along(idx)) {
        beta_i <- sub[i, "beta"]
        for (s in seq_len(nseeds)) {
          r <- sub[i, paste0("delta", s)] / beta_i - 1
          if (!is.na(r) && r > 0) {
            px <- sub[i, paste0("x_", s)]
            py <- sub[i, paste0("y_", s)]
            circle_rows <- c(circle_rows, list(data.frame(
              x     = px + r * cos(theta),
              y     = py + r * sin(theta),
              group = paste(i, s, sep = "_")
            )))
          }
        }
      }

      p <- private$.image$plot() + coord_fixed()

      if (length(circle_rows) > 0) {
        circle_df <- do.call(rbind, circle_rows)
        p <- p + geom_path(
          data        = circle_df,
          aes(x = x, y = y, group = group),
          color       = "white",
          alpha       = 0.12,
          linewidth   = 0.3,
          inherit.aes = FALSE
        )
      }

      center_list <- lapply(seq_len(nseeds), function(s) {
        data.frame(
          x = mean(sub[, paste0("x_", s)]),
          y = mean(sub[, paste0("y_", s)])
        )
      })
      center_df <- do.call(rbind, center_list)
      p + geom_point(
        data        = center_df,
        aes(x = x, y = y),
        shape = 3, size = 3, color = "white", stroke = 1,
        inherit.aes = FALSE
      )
    },

    plot_positions = function(burn_in = NULL) {
      chain  <- private$get_chain(burn_in)
      nseeds <- length(private$.seeds)
      nrows  <- private$.image$dim[1]
      ncols  <- private$.image$dim[2]

      pos_list <- lapply(seq_len(nseeds), function(s) {
        data.frame(
          x    = chain[, paste0("x_", s)],
          y    = chain[, paste0("y_", s)],
          seed = paste0("seed ", s)
        )
      })
      pos_df      <- do.call(rbind, pos_list)
      pos_df$seed <- factor(pos_df$seed)

      init_df <- data.frame(
        x    = sapply(private$.seeds, function(s) s$position[1]),
        y    = sapply(private$.seeds, function(s) s$position[2]),
        seed = factor(paste0("seed ", seq_len(nseeds)))
      )

      ggplot(pos_df, aes(x = x, y = y)) +
        geom_bin2d(binwidth = 1) +
        geom_point(data = init_df, aes(x = x, y = y),
                   shape = 3, size = 3, color = "firebrick", stroke = 1,
                   inherit.aes = FALSE) +
        annotate("rect",
                 xmin = 0.5, xmax = nrows + 0.5,
                 ymin = 0.5, ymax = ncols + 0.5,
                 fill = NA, color = "grey50", linetype = "dashed") +
        scale_fill_viridis_c(option = "magma") +
        coord_fixed() +
        facet_wrap(~seed) +
        labs(x = "x", y = "y", fill = "Count",
             title = "Posterior: seed positions") +
        theme_bw(base_size = 11) +
        theme(strip.text = element_text(face = "bold"))
    },

    plot_logpl_trace = function(burn_in = NULL) {
      n <- nrow(private$.chain)
      if (is.null(burn_in)) burn_in <- floor(0.2 * n)
      ch <- private$get_chain(burn_in)
      df <- data.frame(
        iteration = seq(burn_in + 1L, burn_in + nrow(ch)),
        log_pl    = ch[, "log_pl"]
      )
      ggplot(df, aes(x = iteration, y = log_pl)) +
        geom_line(linewidth = 0.3, alpha = 0.8, color = "steelblue") +
        labs(x = "Iteração", y = "log pseudo-verossimilhança",
             title = "Trace: log pseudo-verossimilhança") +
        theme_bw(base_size = 11)
    }
  ),
  active = list(
    chain             = function() private$.chain,
    image             = function() private$.image,
    seeds             = function() private$.seeds,
    accept_rate_param = function() private$.accept_rate_param,
    accept_rate_pos   = function() private$.accept_rate_pos
  ),
  private = list(
    get_chain = function(burn_in) {
      n <- nrow(private$.chain)
      if (is.null(burn_in)) burn_in <- floor(0.2 * n)
      if (burn_in >= n) stop("burn_in must be less than the chain length")
      private$.chain[seq(burn_in + 1L, n), , drop = FALSE]
    },
    get_param_chain = function(burn_in) {
      private$get_chain(burn_in)[, private$.param_names, drop = FALSE]
    },
    chain_long = function(burn_in) {
      nm     <- private$.param_names
      ch     <- private$get_param_chain(burn_in)
      n_kept <- nrow(ch)
      df <- do.call(rbind, lapply(nm, function(p) {
        data.frame(
          iteration = seq(burn_in + 1L, burn_in + n_kept),
          param     = p,
          value     = ch[, p]
        )
      }))
      df$param <- factor(df$param, levels = nm)
      df
    },
    .chain             = NULL,
    .seeds             = NULL,
    .image             = NULL,
    .accept_rate_param = NULL,
    .accept_rate_pos   = NULL,
    .param_names       = NULL,
    .pos_names         = NULL
  )
)

# Sampler Metropolis-Hastings com posições das seeds como variáveis latentes
seedBayesian <- R6Class(
  "seedBayesian",
  public = list(
    initialize = function(image, seeds_info, ncolors, priors) {
      private$.image      <- image
      private$.seeds      <- self$seeds_kmeans(seeds_info)
      private$.ncolors    <- ncolors
      private$.priors     <- priors
      private$.seedColors <- as.integer(sapply(private$.seeds, function(s) s$color))

      nrows <- image$dim[1]
      ncols <- image$dim[2]
      n     <- nrows * ncols

      private$.pos_x   <- (seq_len(n) - 1L) %% nrows
      private$.pos_y   <- (seq_len(n) - 1L) %/% nrows
      private$.zMatrix <- matrix(as.integer(image$matrix), nrow = nrows, ncol = ncols)
    },

    seeds_kmeans = function(seeds_info) {
      all_seeds <- list()
      
      for (config in seeds_info) {
        color <- config$color
        n_seeds <- config$n_seeds
        
        color_coords <- which(private$.image$matrix == color, arr.ind = TRUE)
        
        if (nrow(color_coords) < n_seeds) {
          stop(sprintf("Not enough color %d pixels (%d) for %d seeds. Need at least %d pixels.",
                       color, nrow(color_coords), n_seeds, n_seeds))
        }
        
        # Média pro caso de uma seed, kmeans pra mais de uma
        if (n_seeds == 1) {
          centroid <- colMeans(color_coords)
          centroids <- matrix(round(centroid, 0), nrow = 1)
        } else {
          kmeans_result <- kmeans(color_coords, centers = n_seeds, iter.max = 100)
          centroids <- round(kmeans_result$centers, 0)
          
          # Ordena as seeds
          distances_from_origin <- rowSums(centroids^2)
          centroid_order <- order(distances_from_origin)
          centroids <- centroids[centroid_order, ]
        }
        
        for (i in seq_len(n_seeds)) {
          all_seeds <- c(all_seeds, list(list(
            position = as.integer(centroids[i, ]),
            color = as.integer(color)
          )))
        }
      }
      
      return(all_seeds)
    },

    # Roda n_iter iterações MH alternando passo de parâmetros e passo de posições
    run = function(n_iter, step_size = 0.1, init = NULL) {
      nseeds     <- length(private$.seeds)
      n_par_cont <- 2L + nseeds
      nrows      <- private$.image$dim[1]
      ncols      <- private$.image$dim[2]
      n_pos_end  <- n_par_cont + 2L * nseeds
      n_total    <- n_pos_end + 1L  # +1 para a coluna log_pl
      p          <- private$.priors

      # inicializar parâmetros contínuos
      if (is.null(init)) {
        message("Initializing from maximum pseudo-likelihood estimate...")
        fitter <- seedModelFitter$new(
          image   = private$.image,
          seeds_info   = seeds_info,
          ncolors = private$.ncolors
        )
        est <- fitter$fit()$estimates
        est <- est[names(est) != "seeds"]
        par <- c(est[["alpha"]],
                 est[["beta"]],
                 est[paste0("delta", seq_len(nseeds))])
      } else {
        par <- c(init$alpha, init$beta, init$deltas)
      }

      # inicializar posições das seeds
      seed_pos <- matrix(
        as.integer(do.call(rbind, lapply(private$.seeds, function(s) s$position))),
        nrow = nseeds, ncol = 2L
      )

      # construir distMatrix inicial
      distMatrix <- matrix(0.0, nrow = length(private$.pos_x), ncol = nseeds)
      for (s in seq_len(nseeds)) {
        distMatrix[, s] <- private$compute_dist_col(seed_pos[s, 1L], seed_pos[s, 2L])
      }

      log_post        <- private$log_posterior(par, distMatrix)
      kept            <- matrix(NA_real_, nrow = n_iter, ncol = n_total)
      update_interval <- max(1L, n_iter %/% 200L)
      pb              <- txtProgressBar(min = 0, max = n_iter, style = 3)

      for (i in seq_len(n_iter)) {

        # Passo A: parâmetros contínuos
        par_prop <- par + rnorm(n_par_cont, 0, step_size)
        lp_prop  <- private$log_posterior(par_prop, distMatrix)

        if (log(runif(1L)) < lp_prop - log_post) {
          par      <- par_prop
          log_post <- lp_prop
        }

        # Passo B: posições das seeds
        for (s in seq_len(nseeds)) {
          z        <- rnorm(2L)
          step     <- as.integer(sign(z) * ceiling(abs(z)))
          pos_prop <- seed_pos[s, ] + step

          # rejeitar propostas fora dos limites (prior uniforme sobre pixels)
          if (pos_prop[1L] < 1L || pos_prop[1L] > nrows ||
              pos_prop[2L] < 1L || pos_prop[2L] > ncols) next

          dist_col_prop        <- private$compute_dist_col(pos_prop[1L], pos_prop[2L])
          distMatrix_prop      <- distMatrix
          distMatrix_prop[, s] <- dist_col_prop
          lp_prop              <- private$log_posterior(par, distMatrix_prop)

          if (log(runif(1L)) < lp_prop - log_post) {
            seed_pos[s, ]   <- pos_prop
            distMatrix[, s] <- dist_col_prop
            log_post        <- lp_prop
          }
        }

        # registrar estado
        kept[i, seq_len(n_par_cont)]             <- c(par[1L],
                                                      par[2L],
                                                      par[seq(3L, 2L + nseeds)])
        kept[i, seq(n_par_cont + 1L, n_pos_end)] <- as.numeric(t(seed_pos))
        kept[i, n_total] <- log_post - (
          dgamma(par[1L], p$shape_alpha, p$rate_alpha, log = TRUE) +
          dgamma(par[2L], p$shape_beta,  p$rate_beta,  log = TRUE) +
          sum(dgamma(par[seq(3L, 2L + nseeds)], p$shape_delta, p$rate_delta, log = TRUE))
        )

        if (i %% update_interval == 0L) setTxtProgressBar(pb, i)
      }
      close(pb)

      colnames(kept) <- c(
        "alpha", "beta", paste0("delta", seq_len(nseeds)),
        as.vector(rbind(paste0("x_", seq_len(nseeds)),
                        paste0("y_", seq_len(nseeds)))),
        "log_pl"
      )

      seedBayesianResult$new(
        chain = kept,
        seeds = private$.seeds,
        image = private$.image
      )
    },

    print = function(...) {
      nseeds <- length(private$.seeds)
      dim    <- private$.image$dim
      p      <- private$.priors
      cat("<seedBayesian>\n")
      cat(sprintf("  image   : %d x %d  (%d colors)\n", dim[1], dim[2], private$.ncolors))
      cat(sprintf("  seeds   : %d  (posições são amostradas)\n", nseeds))
      for (i in seq_len(nseeds)) {
        s <- private$.seeds[[i]]
        cat(sprintf("    [%d] posição inicial = (%d, %d),  color = %d\n",
                    i, s$position[1], s$position[2], s$color))
      }
      cat("  priors  :\n")
      cat(sprintf("    alpha ~ Gamma(%.2f, %.2f)\n", p$shape_alpha, p$rate_alpha))
      cat(sprintf("    beta  ~ Gamma(%.2f, %.2f)\n", p$shape_beta,  p$rate_beta))
      cat(sprintf("    delta ~ Gamma(%.2f, %.2f)\n", p$shape_delta, p$rate_delta))
      cat("    positions  ~ Uniform\n")
      invisible(self)
    }
  ),
  private = list(
    compute_dist_col = function(sx, sy) {
      sqrt((private$.pos_x - (sx - 1L))^2 + (private$.pos_y - (sy - 1L))^2)
    },
    log_posterior = function(par, distMatrix) {
      nseeds <- length(private$.seeds)
      alpha  <- par[1]
      beta   <- par[2]
      deltas <- par[seq(3L, 2L + nseeds)]
      p      <- private$.priors

      lpl <- computeLogPL(
        z          = private$.zMatrix,
        alpha      = alpha,
        beta       = beta,
        ncolors    = private$.ncolors,
        distMatrix = distMatrix,
        seedColors = private$.seedColors,
        seedDeltas = deltas
      )

      lpl +
        dgamma(par[1], p$shape_alpha, p$rate_alpha, log = TRUE) +
        dgamma(par[2], p$shape_beta,  p$rate_beta,  log = TRUE) +
        sum(dgamma(par[seq(3L, 2L + nseeds)], p$shape_delta, p$rate_delta, log = TRUE))
    },
    .image      = NULL,
    .seeds      = NULL,
    .ncolors    = NULL,
    .priors     = NULL,
    .zMatrix    = NULL,
    .seedColors = NULL,
    .pos_x      = NULL,
    .pos_y      = NULL
  )
)

# Sampler Metropolis-Hastings para imagens continuas
seedBayesianCont <- R6Class(
  "seedBayesian",
  public = list(
    initialize = function(image, seeds_info, ncolors, priors) {
      private$.image      <- image
      private$.seeds_info  <- seeds_info
      private$.seeds       <- NULL
      private$.ncolors    <- ncolors
      private$.priors     <- priors
      private$.seedColors  <- NULL
      private$.zMatrix     <- NULL
      private$.mus         <- NULL
      private$.sigmas      <- NULL

      nrows <- image$dim[1]
      ncols <- image$dim[2]
      n     <- nrows * ncols
      private$.pos_x   <- (seq_len(n) - 1L) %% nrows
      private$.pos_y   <- (seq_len(n) - 1L) %/% nrows 
      p <- private$.priors
      alpha_init <- p$shape_alpha / p$rate_alpha
      beta_init  <- p$shape_beta  / p$rate_beta
      
      img_values <- as.numeric(private$.image$matrix)
      km1d <- kmeans(img_values, centers = ncolors, iter.max = 100)
      centers_ord <- order(km1d$centers)
      mus_0 <- as.numeric(km1d$centers[centers_ord])
      
      cl <- km1d$cluster
      sigmas_0 <- numeric(ncolors)
      for (k in seq_len(ncolors)) {
        members <- img_values[cl == which(centers_ord == k)]
        if (length(members) <= 1) {
          sigmas_0[k] <- var(img_values, na.rm = TRUE)
        } else {
          sigmas_0[k] <- var(members, na.rm = TRUE)
        }
        if (is.na(sigmas_0[k]) || sigmas_0[k] <= 0) sigmas_0[k] <- 1e-6
      }

      private$.zMatrix <- private$update_Z(
        Y = private$.image$matrix,
        Z = NULL,
        alpha = alpha_init,
        beta  = beta_init,
        deltas = NULL,
        seeds = NULL,
        mus = mus_0,
        sigmas = sigmas_0
      )
      
      
      private$.mus <- mus_0
      private$.sigmas <- sigmas_0
      private$ensure_seeds_computed()
    },

    seeds_kmeans = function(seeds_info) {
      all_seeds <- list()
      mrf2d::dplot(private$.zMatrix)
      
      for (config in seeds_info) {
        color <- config$color
        n_seeds <- config$n_seeds
        
        z_source <- if (!is.null(private$.zMatrix)) private$.zMatrix else private$.image$matrix
        color_coords <- which(z_source == color, arr.ind = TRUE)

        if (nrow(color_coords) < n_seeds) {
          stop(sprintf("Not enough color %d pixels (%d) for %d seeds. Need at least %d pixels.",
                       color, nrow(color_coords), n_seeds, n_seeds))
        }
        
        # Média pro caso de uma seed, kmeans pra mais de uma
        if (n_seeds == 1) {
          centroid <- colMeans(color_coords)
          centroids <- matrix(round(centroid, 0), nrow = 1)
        } else {
          kmeans_result <- kmeans(color_coords, centers = n_seeds, iter.max = 100)
          centroids <- round(kmeans_result$centers, 0)
          
          # Ordena as seeds
          distances_from_origin <- rowSums(centroids^2)
          centroid_order <- order(distances_from_origin)
          centroids <- centroids[centroid_order, ]
        }
        
        for (i in seq_len(n_seeds)) {
          all_seeds <- c(all_seeds, list(list(
            position = as.integer(centroids[i, ]),
            color = as.integer(color)
          )))
        }
      }
      
      return(all_seeds)
    },
    
    plotZ = function(){
      mrf2d::dplot(private$.zMatrix)
    },

    # Roda n_iter iterações MH alternando passo de parâmetros e passo de posições
    run = function(n_iter, step_size = 0.1, init = NULL) {
      nseeds     <- length(private$.seeds)
      n_par_cont <- 2L + nseeds
      nrows      <- private$.image$dim[1]
      ncols      <- private$.image$dim[2]
      n_pos_end  <- n_par_cont + 2L * nseeds
      n_total    <- n_pos_end + 1L  # +1 para a coluna log_pl
      p          <- private$.priors

      # inicializar parâmetros contínuos
      if (is.null(init)) {
        message("Initializing from maximum pseudo-likelihood estimate...")
        fitter <- seedModelFitter$new(
          image   = private$.zMatrix,
          seeds_info   = private$.seeds_info,
          ncolors = private$.ncolors
        )
        est <- fitter$fit()$estimates
        est <- est[names(est) != "seeds"]
        par <- c(est[["alpha"]],
                 est[["beta"]],
                 est[paste0("delta", seq_len(nseeds))])
      } else {
        par <- c(init$alpha, init$beta, init$deltas)
      }

      # inicializar posições das seeds
      seed_pos <- matrix(
        as.integer(do.call(rbind, lapply(private$.seeds, function(s) s$position))),
        nrow = nseeds, ncol = 2L
      )

      # construir distMatrix inicial
      distMatrix <- matrix(0.0, nrow = length(private$.pos_x), ncol = nseeds)
      for (s in seq_len(nseeds)) {
        distMatrix[, s] <- private$compute_dist_col(seed_pos[s, 1L], seed_pos[s, 2L])
      }

      log_post        <- private$log_posterior(par, distMatrix)
      kept            <- matrix(NA_real_, nrow = n_iter, ncol = n_total)
      update_interval <- max(1L, n_iter %/% 200L)
      pb              <- txtProgressBar(min = 0, max = n_iter, style = 3)
      emission <- list()

      for (i in seq_len(n_iter)) {

        # Passo A: parâmetros contínuos
        par_prop <- par + rnorm(n_par_cont, 0, step_size)
        lp_prop  <- private$log_posterior(par_prop, distMatrix)

        if (log(runif(1L)) < lp_prop - log_post) {
          par      <- par_prop
          log_post <- lp_prop
        }

        # Passo B: posições das seeds
        for (s in seq_len(nseeds)) {
          z        <- rnorm(2L)
          step     <- as.integer(sign(z) * ceiling(abs(z)))
          pos_prop <- seed_pos[s, ] + step

          # rejeitar propostas fora dos limites (prior uniforme sobre pixels)
          if (pos_prop[1L] < 1L || pos_prop[1L] > nrows ||
              pos_prop[2L] < 1L || pos_prop[2L] > ncols) next

          dist_col_prop        <- private$compute_dist_col(pos_prop[1L], pos_prop[2L])
          distMatrix_prop      <- distMatrix
          distMatrix_prop[, s] <- dist_col_prop
          lp_prop              <- private$log_posterior(par, distMatrix_prop)

          if (log(runif(1L)) < lp_prop - log_post) {
            seed_pos[s, ]   <- pos_prop
            distMatrix[, s] <- dist_col_prop
            log_post        <- lp_prop
          }
        }

        # registrar estado
        kept[i, seq_len(n_par_cont)]             <- c(par[1L],
                                                      par[2L],
                                                      par[seq(3L, 2L + nseeds)])
        kept[i, seq(n_par_cont + 1L, n_pos_end)] <- as.numeric(t(seed_pos))
        kept[i, n_total] <- log_post - (
          dgamma(par[1L], p$shape_alpha, p$rate_alpha, log = TRUE) +
          dgamma(par[2L], p$shape_beta,  p$rate_beta,  log = TRUE) +
          sum(dgamma(par[seq(3L, 2L + nseeds)], p$shape_delta, p$rate_delta, log = TRUE))
        )

        # Atualizar Z, mus e sigmas ao final de cada passo, MENOS na última iteração
        if (i < n_iter) {
          seeds_updated <- data.frame(
            x = seed_pos[, 1L],
            y = seed_pos[, 2L],
            cor = sapply(private$.seeds, function(s) s$color)
          )

          # Atualizar Z
          private$.zMatrix <- private$update_Z(
            Y = private$.image$matrix,
            Z = private$.zMatrix,
            alpha = par[1L],
            beta = par[2L],
            deltas = par[seq(3L, 2L + nseeds)],
            seeds = seeds_updated,
            mus = private$.mus,
            sigmas = private$.sigmas
          )
          
          # Atualizar mus e sigmas
          noise_params <- private$update_noise_params(
            Y = private$.image$matrix,
            Z = private$.zMatrix,
            mus = private$.mus,
            sigmas = private$.sigmas,
            priors = p
          )
          private$.mus <- noise_params$mus
          private$.sigmas <- noise_params$sigmas
          emission <- c(emission, list(mu = private$.mus, sigmas = private$.sigmas))

          log_post <- private$log_posterior(par, distMatrix)
        }

        if (i %% update_interval == 0L) setTxtProgressBar(pb, i)
      }
      close(pb)

      colnames(kept) <- c(
        "alpha", "beta", paste0("delta", seq_len(nseeds)),
        as.vector(rbind(paste0("x_", seq_len(nseeds)),
                        paste0("y_", seq_len(nseeds)))),
        "log_pl"
      )
      
      print(emission)

      seedBayesianResult$new(
        chain = kept,
        seeds = private$.seeds,
        image = imageData$new(private$.zMatrix)
      )
    },

    print = function(...) {
      nseeds <- length(private$.seeds)
      dim    <- private$.image$dim
      p      <- private$.priors
      cat("<seedBayesian>\n")
      cat(sprintf("  image   : %d x %d  (%d colors)\n", dim[1], dim[2], private$.ncolors))
      cat(sprintf("  seeds   : %d  (posições são amostradas)\n", nseeds))
      for (i in seq_len(nseeds)) {
        s <- private$.seeds[[i]]
        cat(sprintf("    [%d] posição inicial = (%d, %d),  color = %d\n",
                    i, s$position[1], s$position[2], s$color))
      }
      cat("  priors  :\n")
      cat(sprintf("    alpha ~ Gamma(%.2f, %.2f)\n", p$shape_alpha, p$rate_alpha))
      cat(sprintf("    beta  ~ Gamma(%.2f, %.2f)\n", p$shape_beta,  p$rate_beta))
      cat(sprintf("    delta ~ Gamma(%.2f, %.2f)\n", p$shape_delta, p$rate_delta))
      cat("    positions  ~ Uniform\n")
      invisible(self)
    }
  ),
  private = list(
    compute_dist_col = function(sx, sy) {
      sqrt((private$.pos_x - (sx - 1L))^2 + (private$.pos_y - (sy - 1L))^2)
    },
    update_Z = function(Y, Z = NULL, alpha, beta, deltas = NULL, seeds = NULL,
                             mus, sigmas) {
      n1 <- nrow(Y)
      n2 <- ncol(Y)
      K <- length(mus)
      if(is.null(Z)) Z <- matrix(sample(0:(K-1), n1*n2, replace = TRUE), nrow = n1, ncol = n2)
      if(is.null(deltas)) deltas <- 0.0
      if(is.null(seeds)) seeds <- list(list(position = c(1,1), color = 0))
      
      if(is.data.frame(seeds)){
        seeds_pos <- as.matrix(seeds[,1:2], drop = FALSE)
        seed_colors <- as.vector(seeds[,3])
      } else {
        seeds_pos <- lapply(seeds, \(x) x[["position"]]) |> simplify2array() |> t()
        seed_colors <- sapply(seeds, \(x) x[["color"]])
      }
      
      conditionalGibbsSampler(Y, Z, alpha, beta, K, seeds = seeds_pos, seedColors = seed_colors, seedDeltas = deltas, mus = mus, sigmas = sigmas, steps = 1)
    },
    update_noise_params = function(Y, Z, mus, sigmas, priors) {
      n_colors <- length(mus)
      
      new_mus <- mus
      new_sigmas <- sigmas
      
      for(k in 1:n_colors) {
        y_k <- Y[Z == (k-1)]
        n_k <- length(y_k)
        
        if(n_k == 0) next 
        
        # Atualização do mu
        num_m <- (priors$m[k] / priors$tau[k]) + (sum(y_k) / new_sigmas[k])
        den_m <- (1 / priors$tau[k]) + (n_k / new_sigmas[k])
        m <- num_m / den_m
        tau <- 1 / den_m
                
        new_mus[k] <- rnorm(1, mean = m, sd = sqrt(tau))
        
        # Atualização do sigma
        a <- priors$a[k] + (n_k / 2)
        b <- priors$b[k] + (sum((y_k - new_mus[k])^2) / 2)
        
        new_sigmas[k] <- 1 / rgamma(1, shape = a, rate = b)
      }
      
      return(list(mus = new_mus, sigmas = new_sigmas))
    },
    ensure_seeds_computed = function() {
      if (is.null(private$.seeds) && !is.null(private$.seeds_info)) {
        private$.seeds <- self$seeds_kmeans(private$.seeds_info)
        private$.seedColors <- as.integer(sapply(private$.seeds, function(s) s$color))
      }
    },
    log_posterior = function(par, distMatrix) {
      nseeds <- length(private$.seeds)
      alpha  <- par[1]
      beta   <- par[2]
      deltas <- par[seq(3L, 2L + nseeds)]
      p      <- private$.priors

      lpl <- computeLogPL(
        z          = private$.zMatrix,
        alpha      = alpha,
        beta       = beta,
        ncolors    = private$.ncolors,
        distMatrix = distMatrix,
        seedColors = private$.seedColors,
        seedDeltas = deltas
      )

      lpl +
        dgamma(par[1], p$shape_alpha, p$rate_alpha, log = TRUE) +
        dgamma(par[2], p$shape_beta,  p$rate_beta,  log = TRUE) +
        sum(dgamma(par[seq(3L, 2L + nseeds)], p$shape_delta, p$rate_delta, log = TRUE))
    },
    .image      = NULL,
    .seeds_info = NULL,
    .seeds      = NULL,
    .ncolors    = NULL,
    .priors     = NULL,
    .seedColors = NULL,
    .zMatrix    = NULL,
    .mus        = NULL,
    .sigmas     = NULL,
    .pos_x      = NULL,
    .pos_y      = NULL
  )
)
