# Step 4: species-specific background-climate response surfaces --------------
#
# Response: Delta AUC = warmed AUC minus origin AUC, maximum-duration horizon.
# This script fits one isolated model. It does not add traits, experiment or
# gradient effects, duration, colonizers, or a species-specific warming slope.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "Matrix", "metafor", "readr", "tibble",
  "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

horizon <- Sys.getenv("TRANSPLANT_STEP_AUC_HORIZON", "max")
output_dir <- file.path("output", "step_by_step_modeling", "species_context")
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = horizon)
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z)

auc_by_plot <- readRDS(file.path("output", "intermediate", "auc_by_plot.rds"))
sampling <- reconstruct_auc_covariance(auc_data, auc_by_plot)
if (!isTRUE(sampling$audit$full_covariance_valid[[1]])) {
  stop("The full known AUC-change sampling covariance did not validate.",
       call. = FALSE)
}
V_auc <- sampling$V

# Record the literal requested syntax. metafor rejects ':' notation inside the
# random argument and explicitly requests an interaction variable instead.
literal_syntax_attempt <- auc_capture(
  metafor::rma.mv(
    yi = auc_effect,
    V = V_auc,
    mods = ~ T_z * P_z * W_z,
    random = list(
      ~ 1 | species_id,
      ~ 0 + T_z + P_z + T_z:P_z | species_id
    ),
    struct = "GEN",
    data = auc_data,
    method = "REML",
    sparse = TRUE
  )
)

# Algebraically identical syntax accepted by metafor. TP_z is exactly T_z*P_z;
# no term has been dropped or simplified.
fit_attempt <- auc_capture(
  metafor::rma.mv(
    yi = auc_effect,
    V = V_auc,
    mods = ~ T_z * P_z * W_z,
    random = list(
      ~ 1 | species_id,
      ~ 0 + T_z + P_z + TP_z | species_id
    ),
    struct = "GEN",
    data = auc_data,
    method = "REML",
    sparse = TRUE
  )
)
if (is.null(fit_attempt$fit)) {
  stop("Step-4 species-context model failed: ", fit_attempt$error, call. = FALSE)
}
model_climate_species_context <- fit_attempt$fit

optimizer_code <- auc_value_or(model_climate_species_context$opt.res$convergence,
                               NA_integer_)
coefficient_covariance <- stats::vcov(model_climate_species_context)
converged <- !nzchar(fit_attempt$error) && !length(fit_attempt$warnings) &&
  isTRUE(optimizer_code == 0L) &&
  all(is.finite(stats::coef(model_climate_species_context))) &&
  all(is.finite(coefficient_covariance))

# Verify the actual random-effect blocks created by the fitted object.
G_slopes <- as.matrix(model_climate_species_context$G)
slope_names <- colnames(G_slopes)
if (!identical(slope_names, c("T_z", "P_z", "TP_z"))) {
  stop("Unexpected species-context random-slope columns: ",
       paste(slope_names, collapse = ", "), call. = FALSE)
}
G_full <- matrix(
  0, nrow = 4L, ncol = 4L,
  dimnames = list(c("intercept", slope_names), c("intercept", slope_names))
)
G_full["intercept", "intercept"] <-
  as.numeric(model_climate_species_context$sigma2[[1]])
G_full[slope_names, slope_names] <- G_slopes
R_full <- stats::cov2cor(G_full)

covariance_table <- as.data.frame(as.table(G_full), responseName = "covariance") |>
  tibble::as_tibble() |>
  dplyr::rename(row_coefficient = "Var1",
                column_coefficient = "Var2") |>
  dplyr::mutate(
    status = dplyr::if_else(
      xor(.data$row_coefficient == "intercept",
          .data$column_coefficient == "intercept"),
      "fixed to zero by separate random blocks", "estimated"
    )
  )
correlation_table <- as.data.frame(as.table(R_full), responseName = "correlation") |>
  tibble::as_tibble() |>
  dplyr::rename(row_coefficient = "Var1",
                column_coefficient = "Var2") |>
  dplyr::mutate(
    status = dplyr::if_else(
      xor(.data$row_coefficient == "intercept",
          .data$column_coefficient == "intercept"),
      "fixed to zero by separate random blocks", "estimated"
    )
  )

# Profile-likelihood intervals are expensive because each bound reoptimizes all
# other variances and correlations. A signature-checked cache avoids repeating
# unchanged profiles on every rerun.
variance_specifications <- tibble::tribble(
  ~component, ~argument, ~index,
  "species_intercept", "sigma2", 1L,
  "species_T", "tau2", 1L,
  "species_P", "tau2", 2L,
  "species_TP", "tau2", 3L
)
variance_estimates <- c(
  as.numeric(model_climate_species_context$sigma2[[1]]),
  as.numeric(model_climate_species_context$tau2)
)
profile_signature <- paste(
  format(c(variance_estimates, model_climate_species_context$rho),
         digits = 16, scientific = TRUE),
  collapse = "|"
)
profile_cache_file <- file.path(output_dir, "variance_profile_intervals.rds")
profile_cache <- if (file.exists(profile_cache_file)) readRDS(profile_cache_file) else NULL
profile_cache_valid <- !is.null(profile_cache) &&
  identical(profile_cache$signature, profile_signature) &&
  nrow(profile_cache$table) == nrow(variance_specifications)

if (profile_cache_valid) {
  variance_profiles <- profile_cache$table
} else {
  profile_one <- function(i) {
    specification <- variance_specifications[i, ]
    arguments <- list(object = model_climate_species_context)
    arguments[[specification$argument[[1]]]] <- specification$index[[1]]
    captured <- auc_capture(do.call(stats::confint, arguments))
    if (is.null(captured$fit)) {
      return(tibble::tibble(
        component = specification$component[[1]], estimate = variance_estimates[[i]],
        conf_low = NA_real_, conf_high = NA_real_, warnings = "",
        error = captured$error
      ))
    }
    row_name <- if (specification$argument[[1]] == "sigma2") "sigma^2" else "tau^2"
    interval <- captured$fit$random[row_name, ]
    tibble::tibble(
      component = specification$component[[1]],
      estimate = as.numeric(interval[["estimate"]]),
      conf_low = as.numeric(interval[["ci.lb"]]),
      conf_high = as.numeric(interval[["ci.ub"]]),
      warnings = paste(captured$warnings, collapse = " | "), error = ""
    )
  }
  profile_cores <- suppressWarnings(as.integer(
    Sys.getenv("TRANSPLANT_STEP4_PROFILE_CORES", "4")
  ))
  profile_cores <- max(1L, min(profile_cores, parallel::detectCores()))
  profile_results <- if (.Platform$OS.type != "windows" && profile_cores > 1L) {
    parallel::mclapply(seq_len(nrow(variance_specifications)), profile_one,
                       mc.cores = profile_cores)
  } else {
    lapply(seq_len(nrow(variance_specifications)), profile_one)
  }
  variance_profiles <- dplyr::bind_rows(profile_results)
  invisible(write_rds_atomic(
    list(signature = profile_signature, table = variance_profiles),
    profile_cache_file
  ))
}

boundary_tolerance <- 1e-8 * max(1, stats::var(auc_data$auc_effect))
variance_profiles <- variance_profiles |>
  dplyr::mutate(
    boundary_tolerance = boundary_tolerance,
    finite_profile_interval = is.finite(.data$conf_low) &
      is.finite(.data$conf_high),
    effectively_on_boundary = .data$estimate <= .data$boundary_tolerance |
      !.data$finite_profile_interval | .data$conf_low <= .data$boundary_tolerance
  )

G_eigenvalues <- eigen(G_full, symmetric = TRUE, only.values = TRUE)$values
correlations_estimated <- model_climate_species_context$rho
covariance_positive_definite <- all(is.finite(G_eigenvalues)) &&
  min(G_eigenvalues) > boundary_tolerance
covariance_identifiable <- covariance_positive_definite &&
  all(variance_profiles$finite_profile_interval) &&
  !any(variance_profiles$effectively_on_boundary) &&
  all(abs(correlations_estimated) < 1 - 1e-6)

# Species replication and within-species rank across the requested four-column
# random-coefficient design [1, T, P, T*P].
species_groups <- dplyr::group_split(auc_data, .data$species_id, .keep = TRUE)
species_replication <- dplyr::bind_rows(lapply(species_groups, function(species) {
  context <- dplyr::distinct(species, .data$T_z, .data$P_z)
  design <- cbind(1, context$T_z, context$P_z, context$T_z * context$P_z)
  tibble::tibble(
    species_id = as.character(species$species_id[[1]]),
    n_effects = nrow(species),
    n_experiments = dplyr::n_distinct(species$experiment_id),
    n_distinct_TP_contexts = nrow(context),
    context_design_rank = qr(design)$rank,
    full_rank_four_coefficients = qr(design)$rank == 4L,
    T_span_C = diff(range(species$T)),
    P_span = diff(range(species$P_raw))
  )
})) |>
  dplyr::arrange(dplyr::desc(.data$n_experiments), .data$species_id)

replication_summary <- tibble::tibble(
  n_species = nrow(species_replication),
  n_one_experiment = sum(species_replication$n_experiments == 1L),
  n_two_experiments = sum(species_replication$n_experiments == 2L),
  n_three_experiments = sum(species_replication$n_experiments == 3L),
  n_at_least_four_experiments = sum(species_replication$n_experiments >= 4L),
  median_experiments = stats::median(species_replication$n_experiments),
  max_experiments = max(species_replication$n_experiments),
  n_context_design_rank_1 = sum(species_replication$context_design_rank == 1L),
  n_context_design_rank_2 = sum(species_replication$context_design_rank == 2L),
  n_context_design_rank_3 = sum(species_replication$context_design_rank == 3L),
  n_context_design_rank_4 = sum(species_replication$context_design_rank == 4L)
)

# Conditional BLUPs. ranef.rma.mv does not currently handle this multi-column
# numeric inner formula, so compute the standard D Z' M^-1 residual expression
# from the fitted metafor matrices without refitting or altering the model.
residual <- as.numeric(
  model_climate_species_context$yi -
    model_climate_species_context$X %*% model_climate_species_context$beta
)
weighted_residual <- Matrix::solve(model_climate_species_context$M, residual)
species_levels <- levels(model_climate_species_context$mf.g$outer)
n_species <- length(species_levels)

random_intercepts <- as.numeric(
  model_climate_species_context$sigma2[[1]] *
    Matrix::crossprod(model_climate_species_context$Z.S[[1]], weighted_residual)
)
Z_context <- do.call(cbind, lapply(seq_len(ncol(G_slopes)), function(j) {
  Matrix::Diagonal(
    x = as.numeric(model_climate_species_context$Z.G1[, j])
  ) %*% model_climate_species_context$Z.G2
}))
D_context <- Matrix::kronecker(G_slopes, Matrix::Diagonal(n_species))
random_context_vector <- as.numeric(
  D_context %*% Matrix::crossprod(Z_context, weighted_residual)
)
random_context_matrix <- matrix(
  random_context_vector, nrow = n_species, ncol = length(slope_names),
  dimnames = list(species_levels, slope_names)
)
species_random_coefficients <- tibble::tibble(
  species_id = species_levels,
  random_intercept = random_intercepts,
  random_T = random_context_matrix[, "T_z"],
  random_P = random_context_matrix[, "P_z"],
  random_TP = random_context_matrix[, "TP_z"]
) |>
  dplyr::left_join(species_replication, by = "species_id")

# Plot the six most replicated species that have full-rank T-P support. Each
# surface is masked to that species' own observed two-dimensional convex hull.
display_species <- species_replication |>
  dplyr::filter(.data$full_rank_four_coefficients) |>
  dplyr::slice_max(.data$n_experiments, n = 6L, with_ties = FALSE) |>
  dplyr::arrange(dplyr::desc(.data$n_experiments), .data$species_id)

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
W_z_display <- stats::median(experiment_warming$W_z)
W_display <- stats::median(experiment_warming$W)
beta <- stats::coef(model_climate_species_context)

species_surface_one <- function(species_name, points = 65L) {
  observed <- auc_data |>
    dplyr::filter(as.character(.data$species_id) == species_name) |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z,
                    .data$T, .data$P_raw)
  tile_width <- diff(range(observed$T_z)) / (points - 1L) * 1.02
  tile_height <- diff(range(observed$P_z)) / (points - 1L) * 1.02
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  )
  hull <- geometry::convhulln(as.matrix(observed[, c("T_z", "P_z")]),
                              options = "Qt")
  candidate <- candidate |>
    dplyr::mutate(
      tile_width = tile_width,
      tile_height = tile_height,
      supported = geometry::inhulln(
        hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
      )
    ) |>
    dplyr::filter(.data$supported) |>
    dplyr::mutate(W_z = W_z_display, TP_z = .data$T_z * .data$P_z)
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, candidate, names(beta)
  )
  random <- species_random_coefficients |>
    dplyr::filter(.data$species_id == species_name)
  n_experiments <- random$n_experiments[[1]]
  candidate |>
    dplyr::mutate(
      estimate = as.numeric(fixed_design %*% beta) +
        random$random_intercept[[1]] + random$random_T[[1]] * .data$T_z +
        random$random_P[[1]] * .data$P_z + random$random_TP[[1]] * .data$TP_z,
      species_id = species_name,
      n_experiments = n_experiments,
      panel = paste0(species_name, " (", n_experiments, " experiments)"),
      T_raw = .data$T_z *
        scaled_auc$scaling$scale[scaled_auc$scaling$scaled_field == "T_z"] +
        scaled_auc$scaling$center[scaled_auc$scaling$scaled_field == "T_z"],
      P_raw = exp(
        .data$P_z *
          scaled_auc$scaling$scale[scaled_auc$scaling$scaled_field == "P_z"] +
          scaled_auc$scaling$center[scaled_auc$scaling$scaled_field == "P_z"]
      )
    )
}

species_surfaces <- dplyr::bind_rows(lapply(
  display_species$species_id, species_surface_one
)) |>
  dplyr::mutate(
    panel = factor(
      .data$panel,
      levels = paste0(display_species$species_id, " (",
                      display_species$n_experiments, " experiments)")
    )
  )
surface_observations <- auc_data |>
  dplyr::filter(as.character(.data$species_id) %in% display_species$species_id) |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    panel = factor(
      paste0(.data$species_id, " (",
             display_species$n_experiments[
               match(.data$species_id, display_species$species_id)
             ], " experiments)"),
      levels = levels(species_surfaces$panel)
    )
  )

surface_limit <- max(abs(species_surfaces$estimate), na.rm = TRUE)
T_center <- scaled_auc$scaling$center[
  scaled_auc$scaling$scaled_field == "T_z"
]
T_scale <- scaled_auc$scaling$scale[
  scaled_auc$scaling$scaled_field == "T_z"
]
P_center <- scaled_auc$scaling$center[
  scaled_auc$scaling$scaled_field == "P_z"
]
P_scale <- scaled_auc$scaling$scale[
  scaled_auc$scaling$scaled_field == "P_z"
]
surface_plot <- ggplot2::ggplot(
  species_surfaces,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_tile(
    ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
  ) +
  ggplot2::geom_point(
    data = surface_observations,
    ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
    shape = 21, fill = "white", colour = "grey20", size = 1.5, stroke = 0.35
  ) +
  ggplot2::facet_wrap(~ panel, ncol = 3, scales = "free") +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
    limits = c(-surface_limit, surface_limit)
  ) +
  ggplot2::scale_x_continuous(
    labels = function(x) sprintf("%.1f", x * T_scale + T_center)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::labs(
    title = "Species-specific background-climate response surfaces",
    subtitle = paste0(
      "Predicted Delta AUC at median warming (Delta T = ",
      sprintf("%.2f", W_display), " °C); each panel uses only that species' ",
      "observed T-P convex hull."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Predicted Delta AUC"
  )

figure_stem <- file.path(figure_dir, "step_04_species_tp_surfaces")
ggplot2::ggsave(
  paste0(figure_stem, ".pdf"), surface_plot,
  width = 11, height = 7.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(figure_stem, ".png"), surface_plot,
  width = 11, height = 7.5, dpi = 180, bg = "white"
)

model_summary <- tibble::tibble(
  horizon = horizon,
  converged_without_warning = converged,
  optimizer_code = optimizer_code,
  optimizer_message = auc_value_or(
    model_climate_species_context$opt.res$message, ""
  ),
  warnings = paste(fit_attempt$warnings, collapse = " | "),
  error = fit_attempt$error,
  literal_random_formula_error = literal_syntax_attempt$error,
  actual_random_structure = paste0(
    "block diagonal: independent species intercept plus GEN(T_z, P_z, TP_z)"
  ),
  fully_general_four_by_four = FALSE,
  minimum_covariance_eigenvalue = min(G_eigenvalues),
  maximum_covariance_eigenvalue = max(G_eigenvalues),
  covariance_condition_number = kappa(G_full),
  covariance_positive_definite = covariance_positive_definite,
  covariance_identifiable = covariance_identifiable,
  any_component_on_boundary = any(variance_profiles$effectively_on_boundary),
  figure_pdf = paste0(figure_stem, ".pdf")
)

invisible(write_csv_atomic(
  model_summary, file.path(output_dir, "model_summary.csv")
))
invisible(write_csv_atomic(
  variance_profiles, file.path(output_dir, "variance_profile_intervals.csv")
))
invisible(write_csv_atomic(
  covariance_table, file.path(output_dir, "random_covariance_matrix.csv")
))
invisible(write_csv_atomic(
  correlation_table, file.path(output_dir, "random_correlation_matrix.csv")
))
invisible(write_csv_atomic(
  species_replication, file.path(output_dir, "species_tp_replication.csv")
))
invisible(write_csv_atomic(
  replication_summary, file.path(output_dir, "species_tp_replication_summary.csv")
))
invisible(write_csv_atomic(
  species_random_coefficients,
  file.path(output_dir, "species_random_context_coefficients.csv")
))
invisible(write_csv_atomic(
  species_surfaces,
  file.path(output_dir, "high_replication_species_tp_surfaces.csv")
))
invisible(write_rds_atomic(
  model_climate_species_context,
  file.path(output_dir, "auc_change_species_context_model.rds")
))

cat("Step 4: species-specific T-P background surfaces\n")
cat("Converged without warnings: ", converged, "\n", sep = "")
cat("Actual covariance: independent intercept block + 3x3 GEN slope block\n")
cat("Positive definite: ", covariance_positive_definite,
    "; identifiable: ", covariance_identifiable,
    "; boundary component: ", any(variance_profiles$effectively_on_boundary),
    "\n", sep = "")
cat("Figure: ", normalizePath(paste0(figure_stem, ".pdf")), "\n", sep = "")
