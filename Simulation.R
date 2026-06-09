#setwd("C:/Users/Wesley/Documents/Unicamp/Mestrado/Masters")
source("Functions.R")

# Simulation parameters
gamma1_seq <- seq(3, 6, by = 0.5)
gamma2_seq <- seq(6, 3, by = -0.5)

# Storage for results
simulation_results <- list()
all_grids <- list()
all_inferences <- list()

# Fixed parameters
alpha <- -1
beta <- -1
grid_size <- 50
n_colors <- 2
n_refreshes <- 200
centroids <- matrix(data = c(10, 10, 35, 35), 2, byrow = TRUE)
n_centroids <- 2
n_inference_steps <- 100000
burnin <- 30000

# Run simulation for each gamma1 and gamma2 pair
sim_idx <- 1
for (idx in seq_along(gamma1_seq)) {
  gamma1 <- gamma1_seq[idx]
  gamma2 <- gamma2_seq[idx]
  
  cat("Running simulation", idx, "of", length(gamma1_seq), 
      "- Gamma1 =", gamma1, ", Gamma2 =", gamma2, "\n")
  
  # Generate image with current gamma1 and gamma2
  params <- c(alpha, beta, gamma1, gamma2)
  image_result <- generate_image(grid_size, n_colors, params, centroids, n_refreshes)
  
  # Infer parameters
  inferred_params <- infer_parameters(image_result$grid, n_centroids = n_centroids, n_steps = n_inference_steps)
  
  # Store results
  simulation_results[[sim_idx]] <- list(
    gamma1_true = gamma1,
    gamma2_true = gamma2,
    param_means = c(
      mean(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 1]),
      mean(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 2]),
      mean(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 3]),
      mean(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 4])
    ),
    param_sd = c(
      sd(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 1]),
      sd(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 2]),
      sd(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 3]),
      sd(inferred_params$parameters[burnin:nrow(inferred_params$parameters), 4])
    ),
    inferred_samples = inferred_params
  )
  
  all_grids[[sim_idx]] <- image_result$grid
  all_inferences[[sim_idx]] <- inferred_params
  sim_idx <- sim_idx + 1
}

# Create summary table
summary_table <- data.frame(
  gamma1_true = gamma1_seq,
  gamma2_true = gamma2_seq,
  alpha_mean = sapply(simulation_results, function(x) x$param_means[1]),
  alpha_sd = sapply(simulation_results, function(x) x$param_sd[1]),
  beta_mean = sapply(simulation_results, function(x) x$param_means[2]),
  beta_sd = sapply(simulation_results, function(x) x$param_sd[2]),
  gamma1_mean = sapply(simulation_results, function(x) x$param_means[3]),
  gamma1_sd = sapply(simulation_results, function(x) x$param_sd[3]),
  gamma2_mean = sapply(simulation_results, function(x) x$param_means[4]),
  gamma2_sd = sapply(simulation_results, function(x) x$param_sd[4])
)

print(summary_table)

# Save results
saveRDS(simulation_results, "simulation_results_gammas.rds")
saveRDS(summary_table, "simulation_summary_table.rds")


# Plot 1: Gamma1 estimation
plot(summary_table$gamma1_true, summary_table$gamma1_mean,
     main = "Gamma1 Estimation Across Simulation Study",
     xlab = "True Gamma1",
     ylab = "Estimated Gamma1 (mean)",
     pch = 19,
     col = "blue",
     ylim = c(min(summary_table$gamma1_mean - 1.96*summary_table$gamma1_sd),
              max(summary_table$gamma1_mean + 1.96*summary_table$gamma1_sd)))
arrows(summary_table$gamma1_true, 
       summary_table$gamma1_mean - 1.96*summary_table$gamma1_sd,
       summary_table$gamma1_true,
       summary_table$gamma1_mean + 1.96*summary_table$gamma1_sd,
       angle = 90, code = 3, length = 0.05, col = "blue")
abline(0, 1, col = "red", lty = 2, lwd = 2)
legend("topleft", c("Estimated", "True"), col = c("blue", "red"), lty = c(0, 2), pch = c(19, NA))

# Plot 2: Gamma2 estimation
plot(summary_table$gamma2_true, summary_table$gamma2_mean,
     main = "Gamma2 Estimation Across Simulation Study",
     xlab = "True Gamma2",
     ylab = "Estimated Gamma2 (mean)",
     pch = 19,
     col = "green",
     ylim = c(min(summary_table$gamma2_mean - 1.96*summary_table$gamma2_sd),
              max(summary_table$gamma2_mean + 1.96*summary_table$gamma2_sd)))
arrows(summary_table$gamma2_true, 
       summary_table$gamma2_mean - 1.96*summary_table$gamma2_sd,
       summary_table$gamma2_true,
       summary_table$gamma2_mean + 1.96*summary_table$gamma2_sd,
       angle = 90, code = 3, length = 0.05, col = "green")
abline(0, 1, col = "red", lty = 2, lwd = 2)
legend("topleft", c("Estimated", "True"), col = c("green", "red"), lty = c(0, 2), pch = c(19, NA))

# Plot 3: All parameters across simulation index
par(mfrow = c(2, 2))
par_names <- c("Alpha", "Beta", "Gamma1", "Gamma2")
par_true <- c(alpha, beta, NA, NA)

for (i in 1:4) {
  plot(1:nrow(summary_table), summary_table[, i*2],
       main = paste(par_names[i], "Estimation"),
       xlab = "Simulation Index",
       ylab = paste("Estimated", par_names[i]),
       pch = 19,
       col = "darkblue",
       type = "b")
  arrows(1:nrow(summary_table),
         summary_table[, i*2] - 1.96*summary_table[, i*2+1],
         1:nrow(summary_table),
         summary_table[, i*2] + 1.96*summary_table[, i*2+1],
         angle = 90, code = 3, length = 0.05, col = "darkblue")
  if (!is.na(par_true[i])) {
    abline(h = par_true[i], col = "red", lty = 2, lwd = 2)
  }
}

dev.off()

cat("\n=== Simulation Study Complete ===\n")
cat("Summary table saved to: simulation_summary_table.csv\n")
cat("Full results saved to: simulation_results_gammas.rds\n")
