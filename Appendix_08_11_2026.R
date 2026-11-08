# Appendix S1. 

local({
  # 1. Design and output directories
  
  seed <- 12345L
  nrep <- 400L
  train_fraction <- 0.70
  n_train_rep <- floor(nrep * train_fraction)
  B_boot <- 1000L
  mu0 <- 0
  mu1 <- 1.5
  sigma <- 1
  mu_tau <- 0.75
  sd_tau_grid <- seq(0.05, 0.90, length.out = 12)
  N_grid <- c(200L, 800L)
  n_raters_grid <- c(4L, 8L, 16L)
  prevalence_grid <- c(0.10, 0.50, 0.90)
  lower <- -4
  upper <- 5
  V_ref <- (upper - lower)^2 / 12
  s_grid <- seq(0, 1.50, length.out = 151L)
  quadrature_nodes <- 30L
  quadrature_check_nodes <- c(60L, 100L)
  
  # SCINDEX_OUTPUT_DIR selects the output directory; a compatible simulation is reused.
  root <- Sys.getenv("SCINDEX_OUTPUT_DIR", unset = "EPM_26_1108_results")
  reuse_simulation <- TRUE
  paths <- as.list(file.path(root, c(
    "01_design", "02_simulation", "03_analytic_recovery",
    "04_prevalence", "05_within_regime_prediction", "06_prevalence_transport",
    "07_SDT_calibration", "08_predictive_sensitivity",
    "09_calibration_diagnostics", "10_latent_target_benchmarks",
    "11_figures", "12_session"
  )))
  names(paths) <- c("design", "simulation", "analytic_recovery", "prevalence",
                    "within_regime", "prevalence_transport", "calibration", "sensitivity",
                    "diagnostics", "latent", "figures", "session")
  for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(x, directory, filename) {
    write.csv(x, file.path(directory, filename), row.names = FALSE, na = "")
  }
  compact_model <- function(fit) {
    attr(fit$terms, ".Environment") <- baseenv()
    if (!is.null(fit$model)) attr(fit$model, "terms") <- fit$terms
    if (inherits(fit$call$formula, "formula")) environment(fit$call$formula) <- baseenv()
    fit
  }
  
  simulation_config <- list(seed = seed, nrep = nrep, mu0 = mu0, mu1 = mu1,
                            sigma = sigma, mu_tau = mu_tau, sd_tau_grid = sd_tau_grid, N_grid = N_grid,
                            n_raters_grid = n_raters_grid, prevalence_grid = prevalence_grid,
                            SCI_reference_interval = c(lower, upper), V_ref = V_ref,
                            RNGkind = RNGkind(), R_version = R.version.string, schema_version = 2L)
  analysis_config <- list(train_fraction = train_fraction, B_boot = B_boot,
                          calibration_sd_grid = s_grid, quadrature_nodes = quadrature_nodes,
                          quadrature_check_nodes = quadrature_check_nodes,
                          calibration_prevalence = c("nominal", "observed_current"))
  
  # 2. Simulation
  
  cohen_kappa <- function(a, b) {
    po <- mean(a == b)
    pa <- mean(a)
    pb <- mean(b)
    pe <- pa * pb + (1 - pa) * (1 - pb)
    if (abs(1 - pe) < 1e-10) return(NA_real_)
    (po - pe) / (1 - pe)
  }
  pairwise_kappa <- function(ratings) {
    idx <- combn(ncol(ratings), 2)
    vals <- apply(idx, 2, function(j) {
      cohen_kappa(ratings[, j[1]], ratings[, j[2]])
    })
    mean(vals, na.rm = TRUE)
  }
  new_rater_kappa <- function(new_rating, panel_ratings) {
    vals <- apply(panel_ratings, 2, function(x) {
      cohen_kappa(new_rating, x)
    })
    mean(vals, na.rm = TRUE)
  }
  estimate_threshold <- function(rating, truth) {
    pos <- truth == 1
    neg <- truth == 0
    H <- (sum(rating[pos] == 1) + 0.5) / (sum(pos) + 1)
    F <- (sum(rating[neg] == 1) + 0.5) / (sum(neg) + 1)
    (mu0 + mu1) / 2 - sigma * 0.5 * (qnorm(H) + qnorm(F))
  }
  estimate_thresholds <- function(ratings, truth) {
    apply(ratings, 2, estimate_threshold, truth = truth)
  }
  compute_sci <- function(thresholds) {
    sci <- 1 - var(thresholds) / V_ref
    pmax(0, pmin(1, sci))
  }
  generate_ratings <- function(truth, thresholds) {
    N <- length(truth)
    n_raters <- length(thresholds)
    m <- ifelse(truth == 1, mu1, mu0)
    ratings <- matrix(0L, nrow = N, ncol = n_raters)
    for (i in seq_len(n_raters)) {
      x <- rnorm(N, mean = m, sd = sigma)
      ratings[, i] <- as.integer(x > thresholds[i])
    }
    ratings
  }
  rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
  mcse <- function(x) sd(x[is.finite(x)]) / sqrt(sum(is.finite(x)))
  # One Monte Carlo replication
  run_one <- function(N, n_raters, prevalence, sd_tau) {
    truth_current <- rbinom(N, 1, prevalence)
    truth_future <- rbinom(N, 1, prevalence)
    if (length(unique(truth_current)) < 2 ||
        length(unique(truth_future)) < 2) return(NULL)
    tau <- rnorm(n_raters, mean = mu_tau, sd = sd_tau)
    current_ratings <- generate_ratings(truth_current, tau)
    tau_hat <- estimate_thresholds(current_ratings, truth_current)
    current_kappa <- pairwise_kappa(current_ratings)
    current_sci <- compute_sci(tau_hat)
    latent_current_sci <- compute_sci(tau)
    future_panel <- generate_ratings(truth_future, tau)
    tau_new <- rnorm(1, mean = mu_tau, sd = sd_tau)
    new_rating <- generate_ratings(truth_future, tau_new)[, 1]
    tau_new_hat <- estimate_threshold(new_rating, truth_future)
    data.frame(
      N = N,
      n_raters = n_raters,
      prevalence = prevalence,
      observed_prevalence_current = mean(truth_current),
      observed_prevalence_future = mean(truth_future),
      positive_cases_current = sum(truth_current),
      positive_cases_future = sum(truth_future),
      sd_tau = sd_tau,
      SCI_G = 1 - sd_tau^2 / V_ref,
      current_kappa = current_kappa,
      current_sci = current_sci,
      latent_current_sci = latent_current_sci,
      future_kappa = new_rater_kappa(new_rating, future_panel),
      future_displacement = abs(tau_new_hat - mean(tau_hat)),
      true_displacement = abs(tau_new - mean(tau)),
      tau_new = tau_new,
      tau_bar = mean(tau),
      tau_bar_hat = mean(tau_hat),
      threshold_sd = sd(tau),
      estimated_threshold_sd = sd(tau_hat)
    )
  }
  
  grid <- expand.grid(N = N_grid, n_raters = n_raters_grid,
                      prevalence = prevalence_grid, sd_tau = sd_tau_grid, rep = seq_len(nrep),
                      KEEP.OUT.ATTRS = FALSE)
  grid$replication_id <- seq_len(nrow(grid))
  cache_file <- file.path(paths$simulation, "simulation_results.rds")
  
  if (reuse_simulation && file.exists(cache_file)) {
    simulation <- readRDS(cache_file)
    if (!identical(simulation$config, simulation_config)) {
      stop("Simulation settings differ from the saved data. Select a new SCINDEX_OUTPUT_DIR.")
    }
    results_all <- simulation$replications
    stopifnot(nrow(results_all) == nrow(grid),
              identical(results_all$replication_id, grid$replication_id))
    assign(".Random.seed", simulation$rng_state, envir = .GlobalEnv)
  } else {
    set.seed(seed)
    res <- vector("list", nrow(grid))
    status <- rep("missing_class", nrow(grid))
    pb <- txtProgressBar(min = 0, max = nrow(grid), style = 3)
    for (i in seq_len(nrow(grid))) {
      value <- run_one(grid$N[i], grid$n_raters[i], grid$prevalence[i], grid$sd_tau[i])
      if (!is.null(value)) {
        value$replication_id <- grid$replication_id[i]
        res[[i]] <- value
        status[i] <- "retained"
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    generated <- do.call(rbind, res)
    if (is.null(generated)) stop("No replication contained both case classes.")
    design_fields <- c("N", "n_raters", "prevalence", "sd_tau")
    generated <- generated[, setdiff(names(generated), design_fields), drop = FALSE]
    results_all <- merge(grid, generated, by = "replication_id", all.x = TRUE, sort = FALSE)
    results_all <- results_all[match(grid$replication_id, results_all$replication_id), ]
    finite_fields <- c("current_kappa", "current_sci", "future_kappa",
                       "future_displacement", "true_displacement")
    finite <- rowSums(!is.finite(as.matrix(results_all[, finite_fields]))) == 0L
    status[status == "retained" & !finite] <- "nonfinite"
    results_all$status <- status
    row.names(results_all) <- NULL
    simulation <- list(config = simulation_config, replications = results_all,
                       rng_state = get(".Random.seed", envir = .GlobalEnv))
    temporary_file <- paste0(cache_file, ".tmp")
    saveRDS(simulation, temporary_file)
    if (!file.rename(temporary_file, cache_file)) stop("Cannot save simulation checkpoint.")
    res <- NULL
    generated <- NULL
  }
  dput(simulation_config, file = file.path(paths$design, "simulation_config.R"))
  dput(analysis_config, file = file.path(paths$design, "analysis_config.R"))
  write_table(grid, paths$design, "simulation_design.csv")
  
  results <- results_all[results_all$status == "retained", , drop = FALSE]
  if (!nrow(results)) stop("No finite replications remain.")
  row.names(results) <- NULL
  results$sample <- ifelse(results$rep <= n_train_rep, "train", "validation")
  results$SCI_transform <- sqrt(pmax(0, 1 - results$current_sci))
  results$prospective_scale <- sqrt(V_ref * pmax(0, 1 - results$current_sci) *
                                      (1 + 1 / results$n_raters))
  audit <- data.frame(status = c("design_total", "missing_class", "nonfinite", "retained"),
                      n = c(nrow(grid), sum(results_all$status == "missing_class"),
                            sum(results_all$status == "nonfinite"), nrow(results)))
  stopifnot(sum(audit$n[2:4]) == nrow(grid), !anyDuplicated(results$replication_id))
  write_table(audit, paths$simulation, "replication_counts.csv")
  write_table(results_all[results_all$status != "retained", ],
              paths$simulation, "excluded_replications.csv")
  write_table(data.frame(metric = c("SCI_equal_zero", "SCI_equal_one", "SCI_between_zero_one"),
                         n = c(sum(results$current_sci == 0), sum(results$current_sci == 1),
                               sum(results$current_sci > 0 & results$current_sci < 1))),
              paths$simulation, "SCI_boundary_counts.csv")
  saveRDS(results, file.path(paths$simulation, "analysis_data.rds"))
  write_table(results, paths$simulation, "analysis_data.csv")
  
  # 3. Analytic recovery
  out_dir <- paths$analytic_recovery
  analytic_cells <- do.call(rbind, lapply(
    split(results, interaction(results$n_raters, results$sd_tau, drop = TRUE)),
    function(d) {
      data.frame(
        n_raters = d$n_raters[1],
        sd_tau = d$sd_tau[1],
        SCI_G = d$SCI_G[1],
        mc_mean_D = mean(d$true_displacement),
        mc_mean_D2 = mean(d$true_displacement^2),
        mcse_D = mcse(d$true_displacement),
        theory_mean_D = d$sd_tau[1] * sqrt(1 + 1 / d$n_raters[1]) * sqrt(2 / pi),
        theory_mean_D2 = d$sd_tau[1]^2 * (1 + 1 / d$n_raters[1])
      )
    }
  ))
  analytic_cells$error_D <- analytic_cells$mc_mean_D - analytic_cells$theory_mean_D
  analytic_cells$error_D2 <- analytic_cells$mc_mean_D2 - analytic_cells$theory_mean_D2
  analytic_cells <- analytic_cells[order(analytic_cells$n_raters, analytic_cells$sd_tau), ]
  write.csv(analytic_cells,
            file.path(out_dir, "analytic_prospective_check.csv"),
            row.names = FALSE
  )
  # Recovery of prospective criterion displacement
  out_dir <- paths$analytic_recovery
  recovery <- do.call(rbind, lapply(
    split(results, interaction(results$N, results$prevalence, drop = TRUE)),
    function(d) {
      err <- d$future_displacement - d$true_displacement
      data.frame(
        N = d$N[1],
        prevalence = d$prevalence[1],
        n = nrow(d),
        correlation = cor(d$future_displacement, d$true_displacement),
        bias = mean(err),
        RMSE = sqrt(mean(err^2))
      )
    }
  ))
  write.csv(recovery,
            file.path(out_dir, "criterion_displacement_recovery.csv"),
            row.names = FALSE
  )
  # 4. Prevalence dissociation
  out_dir <- paths$prevalence
  prev_summary <- do.call(rbind, lapply(
    split(results, results$prevalence),
    function(d) {
      data.frame(
        prevalence = d$prevalence[1],
        future_kappa = mean(d$future_kappa),
        true_displacement = mean(d$true_displacement),
        estimated_displacement = mean(d$future_displacement),
        kappa_mcse = mcse(d$future_kappa),
        true_D_mcse = mcse(d$true_displacement),
        estimated_D_mcse = mcse(d$future_displacement)
      )
    }
  ))
  write.csv(prev_summary,
            file.path(out_dir, "prevalence_dissociation_summary.csv"),
            row.names = FALSE
  )
  # 5. Within-regime prediction
  out_dir <- paths$within_regime
  train <- results[results$sample == "train", , drop = FALSE]
  validation <- results[results$sample == "validation", , drop = FALSE]
  M_kappa <- lm(
    future_displacement ~ current_kappa + factor(N) + factor(n_raters) + factor(prevalence),
    data = train
  )
  M_kappa_SCI <- lm(
    future_displacement ~ current_kappa + current_sci +
      factor(N) + factor(n_raters) + factor(prevalence),
    data = train
  )
  M_kappa_transform <- lm(
    future_displacement ~ current_kappa + SCI_transform +
      factor(N) + factor(n_raters) + factor(prevalence),
    data = train
  )
  M_kappa_theory <- lm(
    future_displacement ~ current_kappa + prospective_scale +
      factor(N) + factor(prevalence),
    data = train
  )
  pred_models <- list(
    kappa = predict(M_kappa, newdata = validation),
    kappa_plus_linear_SCI = predict(M_kappa_SCI, newdata = validation),
    kappa_plus_sqrt_1_minus_SCI = predict(M_kappa_transform, newdata = validation),
    kappa_plus_theory_scale = predict(M_kappa_theory, newdata = validation)
  )
  heldout <- data.frame(
    model = names(pred_models),
    RMSE = vapply(pred_models, function(p) rmse(validation$future_displacement, p), numeric(1))
  )
  heldout$delta_RMSE_vs_kappa <- heldout$RMSE[1] - heldout$RMSE
  write.csv(heldout,
            file.path(out_dir, "heldout_model_comparison.csv"),
            row.names = FALSE
  )
  # Paired fixed-fit bootstrap
  set.seed(97531)
  y_valid <- validation$future_displacement
  base_pred <- pred_models$kappa
  heldout_boot <- do.call(rbind, lapply(names(pred_models)[-1], function(nm) {
    alt_pred <- pred_models[[nm]]
    delta <- replicate(B_boot, {
      ii <- sample.int(length(y_valid), replace = TRUE)
      rmse(y_valid[ii], base_pred[ii]) - rmse(y_valid[ii], alt_pred[ii])
    })
    data.frame(
      comparison = paste0("kappa_minus_", nm),
      delta_RMSE = rmse(y_valid, base_pred) - rmse(y_valid, alt_pred),
      low = unname(quantile(delta, 0.025)),
      high = unname(quantile(delta, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(heldout_boot,
            file.path(out_dir, "heldout_paired_RMSE_differences.csv"),
            row.names = FALSE
  )
  capture.output(summary(M_kappa),
                 file = file.path(out_dir, "model_kappa.txt"))
  capture.output(summary(M_kappa_SCI),
                 file = file.path(out_dir, "model_kappa_linear_SCI.txt"))
  capture.output(summary(M_kappa_transform),
                 file = file.path(out_dir, "model_kappa_sqrt_SCI.txt"))
  capture.output(summary(M_kappa_theory),
                 file = file.path(out_dir, "model_kappa_theory_scale.txt"))
  saveRDS(list(models = lapply(list(kappa = M_kappa, linear_SCI = M_kappa_SCI,
                                    transformed_SCI = M_kappa_transform, prospective_SCI = M_kappa_theory), compact_model),
               replication_id = validation$replication_id, predictions = pred_models),
          file.path(out_dir, "models_and_predictions.rds"))
  
  # 6. Prevalence transport
  out_dir <- paths$prevalence_transport
  primary_transport <- list()
  train_balanced <- results[results$prevalence == 0.50 & results$rep <= n_train_rep, ]
  test_extreme <- results[results$prevalence %in% c(0.10, 0.90) & results$rep > n_train_rep, ]
  train_extreme <- results[results$prevalence %in% c(0.10, 0.90) & results$rep <= n_train_rep, ]
  test_balanced <- results[results$prevalence == 0.50 & results$rep > n_train_rep, ]
  # Paired fixed-fit bootstrap
  set.seed(13579)
  transport_boot <- list()
  transport_delta_boot <- list()
  idx_out <- 1L
  idx_delta <- 1L
  for (direction in c("Balanced to extreme", "Extreme to balanced")) {
    if (direction == "Balanced to extreme") {
      tr <- train_balanced
      te <- test_extreme
    } else {
      tr <- train_extreme
      te <- test_balanced
    }
    forms <- list(
      kappa = future_displacement ~ current_kappa + factor(N) + factor(n_raters),
      SCI = future_displacement ~ prospective_scale + factor(N),
      kappa_plus_SCI = future_displacement ~ current_kappa + prospective_scale +
        factor(N)
    )
    fits <- lapply(forms, lm, data = tr)
    preds <- lapply(fits, predict, newdata = te)
    primary_transport[[direction]] <- list(models = lapply(fits, compact_model),
                                           replication_id = te$replication_id, predictions = preds)
    y <- te$future_displacement
    boot_idx <- replicate(B_boot, sample.int(length(y), replace = TRUE), simplify = FALSE)
    for (nm in names(preds)) {
      boot_rmse <- vapply(boot_idx, function(ii) rmse(y[ii], preds[[nm]][ii]), numeric(1))
      transport_boot[[idx_out]] <- data.frame(
        direction = direction,
        model = nm,
        RMSE = rmse(y, preds[[nm]]),
        low = unname(quantile(boot_rmse, 0.025)),
        high = unname(quantile(boot_rmse, 0.975)),
        stringsAsFactors = FALSE
      )
      idx_out <- idx_out + 1L
    }
    for (nm in c("SCI", "kappa_plus_SCI")) {
      boot_delta <- vapply(boot_idx, function(ii) {
        rmse(y[ii], preds$kappa[ii]) - rmse(y[ii], preds[[nm]][ii])
      }, numeric(1))
      transport_delta_boot[[idx_delta]] <- data.frame(
        direction = direction,
        comparison = paste0("kappa_minus_", nm),
        delta_RMSE = rmse(y, preds$kappa) - rmse(y, preds[[nm]]),
        low = unname(quantile(boot_delta, 0.025)),
        high = unname(quantile(boot_delta, 0.975)),
        stringsAsFactors = FALSE
      )
      idx_delta <- idx_delta + 1L
    }
  }
  transport_ci <- do.call(rbind, transport_boot)
  transport_delta_ci <- do.call(rbind, transport_delta_boot)
  write.csv(transport_ci,
            file.path(out_dir, "prospective_transport_RMSE.csv"),
            row.names = FALSE
  )
  write.csv(transport_delta_ci,
            file.path(out_dir, "prospective_transport_paired_RMSE_differences.csv"),
            row.names = FALSE
  )
  saveRDS(primary_transport, file.path(out_dir, "models_and_predictions.rds"))
  
  # Primary figures
  out_dir <- paths$figures
  png(file.path(out_dir, "Figure1_analytic_prospective_relation.png"),
      width = 1300, height = 850, res = 150)
  
  par(mar = c(5, 5, 1.5, 1), las = 1)
  
  xlim <- range(analytic_cells$SCI_G)
  ylim <- range(c(analytic_cells$mc_mean_D,
                  analytic_cells$theory_mean_D))
  
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       xlab = expression(SCI[G]),
       ylab = "Prospective criterion displacement",
       xaxt = "n")
  
  axis(1, at = pretty(xlim))
  
  lty_vals <- c(1, 2, 3)
  pch_vals <- c(1, 2, 16)
  
  n_raters_vals <- sort(unique(analytic_cells$n_raters))
  
  for (j in seq_along(n_raters_vals)) {
    
    d <- analytic_cells[
      analytic_cells$n_raters == n_raters_vals[j], ]
    
    d <- d[order(d$SCI_G), ]
    
    lines(d$SCI_G,
          d$theory_mean_D,
          lty = lty_vals[j],
          lwd = 1.6)
    
    points(d$SCI_G,
           d$mc_mean_D,
           pch = pch_vals[j])
  }
  
  legend("topright",
         legend = paste0("J = ", n_raters_vals),
         lty = lty_vals,
         pch = pch_vals,
         bty = "n")
  
  box()
  dev.off()
  
  # Figure 2: prospective transport across prevalence shifts
  
  png(file.path(out_dir, "Figure2_prospective_transport.png"),
      width = 1300, height = 850, res = 150)
  
  par(mar = c(6, 5, 1.5, 1), las = 1)
  
  model_order <- c("kappa", "SCI", "kappa_plus_SCI")
  direction_order <- c("Balanced to extreme",
                       "Extreme to balanced")
  
  x <- seq_along(model_order)
  offset <- c(-0.12, 0.12)
  pch_dir <- c(1, 16)
  
  ylim <- range(c(transport_ci$low,
                  transport_ci$high))
  
  plot(NA,
       xlim = c(0.5, 3.5),
       ylim = ylim,
       xaxt = "n",
       xlab = "",
       ylab = "Validation RMSE")
  
  axis(1,
       at = x,
       labels = c(expression(kappa),
                  "SCI prospective scale",
                  "kappa + SCI prospective scale"))
  
  for (j in seq_along(direction_order)) {
    
    d <- transport_ci[
      transport_ci$direction == direction_order[j], ]
    
    d <- d[match(model_order, d$model), ]
    
    xx <- x + offset[j]
    
    points(xx,
           d$RMSE,
           pch = pch_dir[j])
    
    arrows(xx,
           d$low,
           xx,
           d$high,
           angle = 90,
           code = 3,
           length = 0.04)
  }
  
  legend("topright",
         legend = c("Balanced to extreme",
                    "Extreme to balanced"),
         pch = pch_dir,
         bty = "n")
  
  box()
  dev.off()
  
  # 7. SDT calibration
  
  r_squared <- function(y, prediction) {
    denominator <- mean((y - mean(y))^2)
    if (denominator <= 0) return(NA_real_)
    1 - mean((y - prediction)^2) / denominator
  }
  gauss_normal <- function(n) {
    A <- matrix(0, n, n)
    v <- sqrt(seq_len(n - 1L) / 2)
    for (i in seq_len(n - 1L)) A[i, i + 1L] <- A[i + 1L, i] <- v[i]
    e <- eigen(A, symmetric = TRUE)
    list(z = sqrt(2) * e$values, w = e$vectors[1, ]^2)
  }
  
  # Conditional independence of the two evidence draws is assumed.
  build_calibration <- function(prevalences, dispersion_grid, nodes) {
    gh <- gauss_normal(nodes)
    pairs <- expand.grid(i = seq_along(gh$z), j = seq_along(gh$z))
    weights <- gh$w[pairs$i] * gh$w[pairs$j]
    values <- matrix(NA_real_, length(dispersion_grid), length(prevalences))
    for (r in seq_along(dispersion_grid)) {
      thresholds <- mu_tau + dispersion_grid[r] * gh$z
      t1 <- thresholds[pairs$i]
      t2 <- thresholds[pairs$j]
      a0 <- pnorm((t1 - mu0) / sigma, lower.tail = FALSE)
      a1 <- pnorm((t1 - mu1) / sigma, lower.tail = FALSE)
      b0 <- pnorm((t2 - mu0) / sigma, lower.tail = FALSE)
      b1 <- pnorm((t2 - mu1) / sigma, lower.tail = FALSE)
      na0 <- pnorm((t1 - mu0) / sigma)
      na1 <- pnorm((t1 - mu1) / sigma)
      nb0 <- pnorm((t2 - mu0) / sigma)
      nb1 <- pnorm((t2 - mu1) / sigma)
      midpoint <- (mu0 + mu1) / 2
      delta1 <- ifelse(t1 >= midpoint, a1 - a0, na0 - na1)
      delta2 <- ifelse(t2 >= midpoint, b1 - b0, nb0 - nb1)
      for (j in seq_along(prevalences)) {
        p <- prevalences[j]
        m1 <- (1 - p) * a0 + p * a1
        m2 <- (1 - p) * b0 + p * b1
        nm1 <- (1 - p) * na0 + p * na1
        nm2 <- (1 - p) * nb0 + p * nb1
        denominator <- m1 * nm2 + nm1 * m2
        if (any(!is.finite(denominator)) || any(denominator <= 0)) {
          stop("Non-positive SDT kappa denominator.")
        }
        values[r, j] <- sum(weights * 2 * p * (1 - p) * delta1 * delta2 / denominator)
      }
    }
    data.frame(calibration_prevalence = rep(prevalences, each = length(dispersion_grid)),
               sd_tau = rep(dispersion_grid, times = length(prevalences)),
               expected_kappa = as.vector(values))
  }
  
  # Every realised current prevalence has its own curve; no prevalence interpolation.
  p_values <- sort(unique(c(results$prevalence, results$observed_prevalence_current)))
  stopifnot(all(p_values > 0 & p_values < 1))
  message("Constructing SDT calibration curves for ", length(p_values), " prevalence values.")
  calibration <- build_calibration(p_values, s_grid, quadrature_nodes)
  write_table(calibration, paths$calibration, "population_kappa_curves.csv")
  saveRDS(calibration, file.path(paths$calibration, "population_kappa_curves.rds"))
  monotonicity <- do.call(rbind, lapply(split(calibration,
                                              match(calibration$calibration_prevalence, p_values)), function(curve) {
                                                data.frame(calibration_prevalence = curve$calibration_prevalence[1],
                                                           maximum_increase = max(c(0, diff(curve$expected_kappa))),
                                                           minimum_decrease = min(-diff(curve$expected_kappa)))
                                              }))
  write_table(monotonicity, paths$calibration, "monotonicity.csv")
  if (any(monotonicity$maximum_increase > 1e-8)) {
    stop("A calibration curve is not numerically monotone.")
  }
  
  invert_calibration <- function(kappa, prevalence_values, curves, dispersion_grid) {
    sd_estimate <- rep(NA_real_, length(kappa))
    region <- rep(NA_character_, length(kappa))
    for (p in sort(unique(prevalence_values))) {
      ix <- which(prevalence_values == p)
      curve <- curves[curves$calibration_prevalence == p, ]
      curve <- curve[order(curve$sd_tau), ]
      stopifnot(nrow(curve) == length(dispersion_grid),
                all(diff(curve$expected_kappa) <= 1e-8))
      maximum <- curve$expected_kappa[1]
      minimum <- tail(curve$expected_kappa, 1)
      high <- kappa[ix] > maximum
      low <- kappa[ix] < minimum
      clipped <- pmax(minimum, pmin(maximum, kappa[ix]))
      estimate <- approx(rev(curve$expected_kappa), rev(curve$sd_tau),
                         xout = clipped, ties = mean, rule = 2)$y
      estimate[kappa[ix] >= maximum] <- min(dispersion_grid)
      estimate[kappa[ix] <= minimum] <- max(dispersion_grid)
      sd_estimate[ix] <- estimate
      region[ix] <- ifelse(high, "above_kappa_range",
                           ifelse(low, "below_kappa_range", "interior"))
    }
    stopifnot(all(is.finite(sd_estimate)), !anyNA(region))
    data.frame(sd_estimate = sd_estimate, calibration_region = region)
  }
  nominal <- invert_calibration(results$current_kappa, results$prevalence,
                                calibration, s_grid)
  observed <- invert_calibration(results$current_kappa,
                                 results$observed_prevalence_current, calibration, s_grid)
  predictive_data <- results
  predictive_data$kappa_SDT_nominal_scale <- nominal$sd_estimate *
    sqrt(1 + 1 / results$n_raters)
  predictive_data$kappa_SDT_observed_scale <- observed$sd_estimate *
    sqrt(1 + 1 / results$n_raters)
  predictive_data$nominal_calibration_region <- nominal$calibration_region
  predictive_data$observed_calibration_region <- observed$calibration_region
  saveRDS(predictive_data, file.path(paths$calibration, "calibrated_analysis_data.rds"))
  write_table(data.frame(replication_id = results$replication_id,
                         nominal_prevalence = results$prevalence,
                         observed_prevalence_current = results$observed_prevalence_current,
                         sd_nominal = nominal$sd_estimate, sd_observed = observed$sd_estimate,
                         region_nominal = nominal$calibration_region,
                         region_observed = observed$calibration_region),
              paths$calibration, "replication_calibration.csv")
  
  # Quadrature diagnostics cover nominal prevalences and selected realised values.
  observed_values <- sort(unique(results$observed_prevalence_current))
  check_p <- sort(unique(c(prevalence_grid,
                           observed_values[c(1L, ceiling(length(observed_values) / 2), length(observed_values))])))
  quadrature_curves <- list()
  for (nodes in c(quadrature_nodes, quadrature_check_nodes)) {
    curves <- if (nodes == quadrature_nodes) {
      calibration[calibration$calibration_prevalence %in% check_p, ]
    } else build_calibration(check_p, s_grid, nodes)
    curves$nodes <- nodes
    quadrature_curves[[as.character(nodes)]] <- curves
  }
  quadrature_reference <- quadrature_curves[[as.character(max(quadrature_check_nodes))]]
  quadrature_checks <- do.call(rbind, lapply(quadrature_curves, function(curves) {
    paired <- merge(curves, quadrature_reference,
                    by = c("calibration_prevalence", "sd_tau"), suffixes = c("", "_reference"))
    do.call(rbind, lapply(split(paired, paired$calibration_prevalence), function(x) {
      data.frame(nodes = x$nodes[1], reference_nodes = x$nodes_reference[1],
                 calibration_prevalence = x$calibration_prevalence[1],
                 maximum_absolute_kappa_difference = max(abs(x$expected_kappa - x$expected_kappa_reference)))
    }))
  }))
  write_table(quadrature_checks, paths$calibration, "quadrature_comparison.csv")
  
  # 8. Predictive sensitivity
  
  train <- predictive_data[predictive_data$sample == "train", ]
  validation <- predictive_data[predictive_data$sample == "validation", ]
  sets <- list(
    within = list(train = train, test = validation),
    balanced_to_extreme = list(train = subset(train, prevalence == .50),
                               test = subset(validation, prevalence %in% c(.10, .90))),
    extreme_to_balanced = list(train = subset(train, prevalence %in% c(.10, .90)),
                               test = subset(validation, prevalence == .50))
  )
  forms <- list(
    kappa_original = future_displacement ~ current_kappa + factor(N) + factor(n_raters),
    kappa_flexible = future_displacement ~
      splines::ns(current_kappa, df = 4, Boundary.knots = c(-1, 1)) *
      factor(n_raters) + factor(N),
    kappa_SDT_nominal = future_displacement ~ kappa_SDT_nominal_scale + factor(N),
    kappa_SDT_observed = future_displacement ~ kappa_SDT_observed_scale + factor(N),
    SCI_prospective_scale = future_displacement ~ prospective_scale + factor(N)
  )
  forms_within <- forms
  forms_within$kappa_original <- future_displacement ~ current_kappa + factor(N) +
    factor(n_raters) + factor(prevalence)
  forms_within$kappa_flexible <- future_displacement ~
    splines::ns(current_kappa, df = 4, Boundary.knots = c(-1, 1)) *
    factor(n_raters) + factor(N) + factor(prevalence)
  
  score_models <- function(y, predictions, direction, target) {
    do.call(rbind, lapply(names(predictions), function(model) {
      data.frame(direction = direction, target = target, model = model, n = length(y),
                 RMSE = rmse(y, predictions[[model]]), R2 = r_squared(y, predictions[[model]]))
    }))
  }
  paired_bootstrap <- function(y, predictions, bootstrap_seed,
                               reference = "SCI_prospective_scale") {
    set.seed(bootstrap_seed)
    alternatives <- setdiff(names(predictions), reference)
    estimates <- matrix(NA_real_, B_boot, length(alternatives))
    for (b in seq_len(B_boot)) {
      ix <- sample.int(length(y), replace = TRUE)
      reference_error <- rmse(y[ix], predictions[[reference]][ix])
      estimates[b, ] <- vapply(alternatives, function(model) {
        rmse(y[ix], predictions[[model]][ix]) - reference_error
      }, numeric(1))
    }
    reference_error <- rmse(y, predictions[[reference]])
    do.call(rbind, lapply(seq_along(alternatives), function(j) {
      alternative_error <- rmse(y, predictions[[alternatives[j]]])
      data.frame(model = alternatives[j], reference = reference,
                 delta_RMSE = alternative_error - reference_error,
                 percent_RMSE_reduction = 100 * (alternative_error - reference_error) / alternative_error,
                 CI_lower = unname(quantile(estimates[, j], .025)),
                 CI_upper = unname(quantile(estimates[, j], .975)))
    }))
  }
  
  estimated_scores <- list()
  estimated_contrasts <- list()
  subgroups <- list()
  boundary_scores <- list()
  boundary_contrasts <- list()
  latent_scores <- list()
  latent_contrasts <- list()
  model_outputs <- list()
  for (direction in names(sets)) {
    tr <- sets[[direction]]$train
    te <- sets[[direction]]$test
    stopifnot(!any(tr$replication_id %in% te$replication_id))
    spec <- if (direction == "within") forms_within else forms
    fits <- lapply(spec, lm, data = tr)
    predictions <- lapply(fits, function(fit) as.numeric(predict(fit, newdata = te)))
    stopifnot(all(vapply(predictions, function(x) all(is.finite(x)), logical(1))))
    y <- te$future_displacement
    estimated_scores[[direction]] <- score_models(y, predictions, direction, "estimated_displacement")
    contrasts <- paired_bootstrap(y, predictions,
                                  switch(direction, within = 51011L, balanced_to_extreme = 51012L,
                                         extreme_to_balanced = 51013L))
    contrasts$direction <- direction
    contrasts$target <- "estimated_displacement"
    estimated_contrasts[[direction]] <- contrasts
    for (variable in c("N", "n_raters", "sd_tau", "prevalence")) {
      for (level in sort(unique(te[[variable]]))) {
        ix <- which(te[[variable]] == level)
        part <- score_models(y[ix], lapply(predictions, function(x) x[ix]),
                             direction, "estimated_displacement")
        part$group <- variable
        part$level <- level
        subgroups[[length(subgroups) + 1L]] <- part
      }
    }
    
    # 9. Calibration-region diagnostics
    for (source in c("nominal", "observed")) {
      region <- te[[paste0(source, "_calibration_region")]]
      calibrated_model <- paste0("kappa_SDT_", source)
      for (category in c("all", "interior", "above_kappa_range", "below_kappa_range")) {
        ix <- if (category == "all") seq_len(nrow(te)) else which(region == category)
        if (!length(ix)) next
        part_predictions <- lapply(predictions, function(x) x[ix])
        part <- score_models(y[ix], part_predictions, direction, "estimated_displacement")
        part$calibration_source <- source
        part$calibration_region <- category
        boundary_scores[[length(boundary_scores) + 1L]] <- part
        if (category == "all") {
          contrast <- contrasts[contrasts$model == calibrated_model, ]
          contrast$direction <- NULL
          contrast$target <- NULL
        } else {
          pair <- part_predictions[c(calibrated_model, "SCI_prospective_scale")]
          contrast <- paired_bootstrap(y[ix], pair,
                                       62000L + match(direction, names(sets)) * 100L +
                                         match(source, c("nominal", "observed")) * 10L +
                                         match(category, c("all", "interior", "above_kappa_range", "below_kappa_range")))
        }
        contrast$direction <- direction
        contrast$calibration_source <- source
        contrast$calibration_region <- category
        contrast$n <- length(ix)
        boundary_contrasts[[length(boundary_contrasts) + 1L]] <- contrast
      }
    }
    
    # 10. Latent-target models and oracle benchmarks
    latent_spec <- lapply(spec, function(f) {
      f[[2]] <- as.name("true_displacement")
      f
    })
    latent_fits <- lapply(latent_spec, lm, data = tr)
    latent_predictions <- lapply(latent_fits,
                                 function(fit) as.numeric(predict(fit, newdata = te)))
    a <- abs(te$tau_bar - mu_tau)
    z <- a / te$sd_tau
    oracles <- list(
      marginal_oracle = te$sd_tau * sqrt(1 + 1 / te$n_raters) * sqrt(2 / pi),
      conditional_oracle = te$sd_tau * sqrt(2 / pi) * exp(-z^2 / 2) +
        a * (2 * pnorm(z) - 1))
    stopifnot(all(vapply(c(latent_predictions, oracles),
                         function(x) all(is.finite(x)), logical(1))))
    latent_scores[[direction]] <- score_models(te$true_displacement,
                                               c(latent_predictions, oracles), direction, "true_displacement")
    contrast <- paired_bootstrap(te$true_displacement, latent_predictions,
                                 71000L + match(direction, names(sets)))
    contrast$direction <- direction
    contrast$target <- "true_displacement"
    latent_contrasts[[direction]] <- contrast
    model_outputs[[direction]] <- list(replication_id = te$replication_id,
                                       estimated = list(models = lapply(fits, compact_model), predictions = predictions),
                                       latent = list(models = lapply(latent_fits, compact_model),
                                                     predictions = latent_predictions, oracles = oracles))
  }
  estimated_scores <- do.call(rbind, estimated_scores)
  estimated_contrasts <- do.call(rbind, estimated_contrasts)
  latent_scores <- do.call(rbind, latent_scores)
  write_table(estimated_scores, paths$sensitivity, "validation_scores.csv")
  write_table(estimated_contrasts, paths$sensitivity, "paired_RMSE_differences.csv")
  write_table(do.call(rbind, subgroups), paths$sensitivity, "subgroup_scores.csv")
  saveRDS(model_outputs, file.path(paths$sensitivity, "models_and_predictions.rds"))
  
  boundary_counts <- do.call(rbind, lapply(c("nominal", "observed"), function(source) {
    counts <- aggregate(replication_id ~ N + n_raters + prevalence + region,
                        data = transform(predictive_data,
                                         region = predictive_data[[paste0(source, "_calibration_region")]]), FUN = length)
    names(counts)[names(counts) == "replication_id"] <- "n"
    counts$calibration_source <- source
    counts
  }))
  write_table(boundary_counts, paths$diagnostics, "calibration_region_counts.csv")
  write_table(do.call(rbind, boundary_scores), paths$diagnostics, "calibration_region_scores.csv")
  write_table(do.call(rbind, boundary_contrasts), paths$diagnostics, "calibration_region_RMSE_differences.csv")
  write_table(latent_scores, paths$latent, "validation_scores.csv")
  write_table(do.call(rbind, latent_contrasts), paths$latent, "paired_RMSE_differences.csv")
  
  # Robustness
  
  extended_grid <- seq(0, 2.50, by = 0.01)
  stopifnot(isTRUE(all.equal(s_grid, seq(0, 1.50, by = 0.01))))
  extended_curves <- build_calibration(p_values, extended_grid, quadrature_nodes)
  maximum_increase <- vapply(split(extended_curves, extended_curves$calibration_prevalence),
                             function(x) max(c(0, diff(x$expected_kappa))), numeric(1))
  stopifnot(all(maximum_increase <= 1e-8))
  
  nominal_extended <- invert_calibration(results$current_kappa, results$prevalence,
                                         extended_curves, extended_grid)
  observed_extended <- invert_calibration(results$current_kappa,
                                          results$observed_prevalence_current,
                                          extended_curves, extended_grid)
  
  d_robust <- predictive_data
  d_robust$kappa_SDT_nominal_extended_scale <- nominal_extended$sd_estimate *
    sqrt(1 + 1 / results$n_raters)
  d_robust$kappa_SDT_observed_extended_scale <- observed_extended$sd_estimate *
    sqrt(1 + 1 / results$n_raters)
  
  robustness_forms <- list(
    kappa_SDT_nominal = future_displacement ~ kappa_SDT_nominal_scale + factor(N),
    kappa_SDT_observed = future_displacement ~ kappa_SDT_observed_scale + factor(N),
    kappa_SDT_nominal_extended = future_displacement ~ kappa_SDT_nominal_extended_scale + factor(N),
    kappa_SDT_observed_extended = future_displacement ~ kappa_SDT_observed_extended_scale + factor(N),
    SCI_prospective_scale = future_displacement ~ prospective_scale + factor(N)
  )
  
  robustness_sets <- list(
    within = list(train = subset(d_robust, sample == "train"),
                  test = subset(d_robust, sample == "validation")),
    balanced_to_extreme = list(
      train = subset(d_robust, sample == "train" & prevalence == .50),
      test = subset(d_robust, sample == "validation" & prevalence %in% c(.10, .90))),
    extreme_to_balanced = list(
      train = subset(d_robust, sample == "train" & prevalence %in% c(.10, .90)),
      test = subset(d_robust, sample == "validation" & prevalence == .50))
  )
  
  robustness_scores <- list()
  robustness_contrasts <- list()
  for (direction in names(robustness_sets)) {
    tr <- robustness_sets[[direction]]$train
    te <- robustness_sets[[direction]]$test
    stopifnot(!any(tr$replication_id %in% te$replication_id))
    fits <- lapply(robustness_forms, lm, data = tr)
    predictions <- lapply(fits, function(fit) as.numeric(predict(fit, newdata = te)))
    stopifnot(all(vapply(predictions, function(x) all(is.finite(x)), logical(1))))
    y <- te$future_displacement
    robustness_scores[[direction]] <- data.frame(
      direction = direction, model = names(predictions), n = length(y),
      RMSE = vapply(predictions, function(x) rmse(y, x), numeric(1)))
    contrast <- paired_bootstrap(y, predictions,
                                 83000L + match(direction, names(robustness_sets)))
    contrast$direction <- direction
    robustness_contrasts[[direction]] <- contrast
  }
  
  robustness_scores <- do.call(rbind, robustness_scores)
  robustness_contrasts <- do.call(rbind, robustness_contrasts)
  rownames(robustness_scores) <- NULL
  rownames(robustness_contrasts) <- NULL
  
  original_scores <- estimated_scores[
    estimated_scores$model %in% c("kappa_SDT_nominal", "kappa_SDT_observed",
                                  "SCI_prospective_scale"),
    c("direction", "model", "RMSE")]
  ix <- match(paste(original_scores$direction, original_scores$model),
              paste(robustness_scores$direction, robustness_scores$model))
  stopifnot(!anyNA(ix),
            max(abs(robustness_scores$RMSE[ix] - original_scores$RMSE)) < 1e-10)
  
  boundary_robustness <- do.call(rbind, lapply(c("nominal", "observed"), function(src) {
    old <- if (src == "nominal") nominal else observed
    new <- if (src == "nominal") nominal_extended else observed_extended
    data.frame(source = src, original_region = old$calibration_region,
               extended_region = new$calibration_region,
               delta_sd = new$sd_estimate - old$sd_estimate)
  }))
  
  boundary_summary <- aggregate(
    cbind(n = rep(1L, nrow(boundary_robustness)),
          changed = as.integer(abs(boundary_robustness$delta_sd) > 1e-10)),
    by = boundary_robustness[c("source", "original_region", "extended_region")],
    FUN = sum)
  
  sd_change <- do.call(rbind, lapply(c("nominal", "observed"), function(src) {
    z <- boundary_robustness[boundary_robustness$source == src, , drop = FALSE]
    data.frame(source = src, n_changed = sum(abs(z$delta_sd) > 1e-10),
               maximum_absolute_sd_change = max(abs(z$delta_sd)),
               mean_absolute_sd_change = mean(abs(z$delta_sd)))
  }))
  
  stopifnot(all(vapply(c("nominal", "observed"), function(src) {
    sum(boundary_summary$n[boundary_summary$source == src]) == nrow(results) &&
      sum(boundary_summary$changed[boundary_summary$source == src]) ==
      sd_change$n_changed[sd_change$source == src]
  }, logical(1))))
  
  write_table(robustness_scores, paths$diagnostics, "robustness_RMSE_by_grid_endpoint.csv")
  write_table(robustness_contrasts, paths$diagnostics, "robustness_paired_RMSE_differences.csv")
  write_table(boundary_summary, paths$diagnostics, "robustness_boundary_transitions.csv")
  write_table(sd_change, paths$diagnostics, "robustness_inversion_change_summary.csv")
  write_table(data.frame(replication_id = results$replication_id,
                         nominal_old_sd = nominal$sd_estimate,
                         nominal_extended_sd = nominal_extended$sd_estimate,
                         observed_old_sd = observed$sd_estimate,
                         observed_extended_sd = observed_extended$sd_estimate),
              paths$diagnostics, "robustness_replication_inversion_comparison.csv")
  
  # 11. Calibrated transport comparison
  
  png(file.path(paths$figures, "Figure3_calibrated_transport.png"),
      width = 1600, height = 950, res = 150)
  par(mar = c(8, 5, 2, 1), las = 1)
  models <- c("kappa_original", "kappa_SDT_nominal", "kappa_SDT_observed", "SCI_prospective_scale")
  directions <- c("balanced_to_extreme", "extreme_to_balanced")
  transport_scores <- estimated_scores[estimated_scores$direction %in% directions &
                                         estimated_scores$model %in% models, ]
  plot(NA, xlim = c(.5, 4.5), ylim = range(transport_scores$RMSE), xaxt = "n",
       xlab = "", ylab = "Validation RMSE")
  axis(1, at = seq_along(models), labels = c("Kappa", "Kappa-SDT\nnominal prevalence",
                                             "Kappa-SDT\nobserved prevalence", "SCI prospective\nscale"), cex.axis = .85)
  for (j in seq_along(directions)) {
    values <- transport_scores[transport_scores$direction == directions[j], ]
    values <- values[match(models, values$model), ]
    points(seq_along(models) + c(-.10, .10)[j], values$RMSE, pch = c(1, 16)[j])
  }
  legend("topright", legend = c("Balanced to extreme", "Extreme to balanced"),
         pch = c(1, 16), bty = "n")
  dev.off()
  
  # 12. Reproducibility
  
  primary_scores <- data.frame(direction = "within", model = "kappa_original",
                               RMSE = heldout$RMSE[heldout$model == "kappa"])
  transport_reference <- transport_ci[transport_ci$model %in% c("kappa", "SCI"),
                                      c("direction", "model", "RMSE")]
  transport_reference$direction <- ifelse(transport_reference$direction == "Balanced to extreme",
                                          "balanced_to_extreme", "extreme_to_balanced")
  transport_reference$model <- ifelse(transport_reference$model == "kappa",
                                      "kappa_original", "SCI_prospective_scale")
  primary_scores <- rbind(primary_scores, transport_reference)
  shared_scores <- merge(primary_scores, estimated_scores,
                         by = c("direction", "model"), suffixes = c("_primary", "_sensitivity"))
  shared_scores$absolute_RMSE_difference <- abs(shared_scores$RMSE_primary -
                                                  shared_scores$RMSE_sensitivity)
  stopifnot(nrow(shared_scores) == nrow(primary_scores),
            max(shared_scores$absolute_RMSE_difference) < 1e-10)
  write_table(shared_scores, paths$session, "shared_model_checks.csv")
  stopifnot(all(results$observed_prevalence_current == results$positive_cases_current / results$N),
            all(results$observed_prevalence_future == results$positive_cases_future / results$N))
  interior_sci <- results$current_sci > 0 & results$current_sci < 1
  scale_error <- max(abs(results$prospective_scale[interior_sci] -
                           results$estimated_threshold_sd[interior_sci] *
                           sqrt(1 + 1 / results$n_raters[interior_sci])))
  stopifnot(scale_error < 1e-10)
  write_table(data.frame(check = c("analytic_mean_maximum_absolute_error",
                                   "analytic_squared_mean_maximum_absolute_error", "prospective_scale_identity_error"),
                         value = c(max(abs(analytic_cells$error_D)), max(abs(analytic_cells$error_D2)), scale_error)),
              paths$session, "numerical_checks.csv")
  capture.output(sessionInfo(), file = file.path(paths$session, "sessionInfo.txt"))
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[basename(files) != "output_files.csv"]
  write_table(data.frame(file = substring(files, nchar(root) + 2L),
                         bytes = file.info(files)$size), paths$session, "output_files.csv")
  message("Analysis complete: ", normalizePath(root))
})