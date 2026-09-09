# Step-by-step resident AUC-change modeling -----------------------------------
#
# Response in both steps: Delta AUC = warmed AUC minus origin AUC.
# Step 1: Delta AUC ~ 1 with known sampling covariance and a species intercept.
# Step 2: Delta AUC ~ T_z * P_z * W_z with the same covariance and hierarchy.
#
# This script is intentionally independent of the large visual-report workflow.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "Matrix", "metafor", "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

horizon <- Sys.getenv("TRANSPLANT_STEP_AUC_HORIZON", "max")
output_dir <- file.path("output", "step_by_step_modeling", "auc_intercept_species")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare resident raw warmed-minus-origin AUC changes for one horizon. Climate
# scaling uses the distinct experiments, exactly as in the main analysis.
scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = horizon)
)
auc_data <- scaled_auc$data

# Reconstruct the known sampling covariance from the original per-experiment
# AUC contrasts. Its diagonal is the reported auc_se^2.
auc_by_plot <- readRDS(file.path("output", "intermediate", "auc_by_plot.rds"))
sampling <- reconstruct_auc_covariance(auc_data, auc_by_plot)
if (!isTRUE(sampling$audit$full_covariance_valid[[1]])) {
  stop("The full known AUC sampling covariance did not validate.", call. = FALSE)
}
V_auc <- sampling$V

# Exact model requested for this step.
fit_attempt <- auc_capture(
  metafor::rma.mv(
    yi = auc_effect,
    V = V_auc,
    mods = ~ 1,
    random = ~ 1 | species_id,
    data = auc_data,
    method = "REML",
    sparse = TRUE
  )
)
if (is.null(fit_attempt$fit)) {
  stop("Intercept-only species model failed: ", fit_attempt$error, call. = FALSE)
}
model <- fit_attempt$fit

coefficient_covariance <- stats::vcov(model)
converged <- !nzchar(fit_attempt$error) && !length(fit_attempt$warnings) &&
  all(is.finite(stats::coef(model))) && all(is.finite(coefficient_covariance))

# Profile-likelihood interval for the species variance.
variance_interval_attempt <- auc_capture(
  stats::confint(model, sigma2 = 1, progbar = FALSE)
)
variance_interval <- if (is.null(variance_interval_attempt$fit)) {
  matrix(NA_real_, nrow = 2, ncol = 3,
         dimnames = list(c("sigma^2", "sigma"), c("estimate", "ci.lb", "ci.ub")))
} else {
  variance_interval_attempt$fit$random
}
species_variance <- as.numeric(model$sigma2[[1]])
species_variance_estimable <- is.finite(species_variance) && species_variance > 0 &&
  is.finite(variance_interval["sigma^2", "ci.lb"]) &&
  is.finite(variance_interval["sigma^2", "ci.ub"]) &&
  variance_interval["sigma^2", "ci.lb"] > 0

# Influence diagnostics. Exact case-deletion Cook distances would refit 1,625
# models, so this step reports the conventional leverage/studentized-residual
# approximation and labels it explicitly.
standardized_residual <- as.numeric(metafor::rstandard.rma.mv(model)$z)
leverage <- as.numeric(stats::hatvalues(model))
cooks_distance_approx <- standardized_residual^2 * leverage /
  pmax(1 - leverage, 1e-8)

# For an intercept-only rma.mv model, row-sum weights are the actual percentage
# contributions to the fitted overall mean under the estimated marginal model.
mean_weight_percent <- as.numeric(stats::weights(model, type = "rowsum"))
diagnostics <- auc_data |>
  dplyr::transmute(
    .data$effect_id, species_id = as.character(.data$species_id),
    experiment_id = as.character(.data$experiment_id),
    .data$auc_effect, .data$auc_se,
    mean_weight_percent = .env$mean_weight_percent,
    standardized_residual = .env$standardized_residual,
    leverage = .env$leverage,
    cooks_distance_approx = .env$cooks_distance_approx,
    extreme_abs_z_above_3 = abs(.env$standardized_residual) > 3,
    influential_cook_above_4_over_n = .env$cooks_distance_approx > 4 / nrow(auc_data)
  )

se_10 <- as.numeric(stats::quantile(auc_data$auc_se, 0.10))
se_90 <- as.numeric(stats::quantile(auc_data$auc_se, 0.90))
low_se_median_weight <- stats::median(
  mean_weight_percent[auc_data$auc_se <= se_10]
)
high_se_median_weight <- stats::median(
  mean_weight_percent[auc_data$auc_se >= se_90]
)

summary_table <- tibble::tibble(
  horizon = horizon,
  n_effects = nrow(auc_data),
  n_species = dplyr::n_distinct(auc_data$species_id),
  n_experiments = dplyr::n_distinct(auc_data$experiment_id),
  converged_without_warning = converged,
  overall_mean_auc = as.numeric(stats::coef(model)[[1]]),
  overall_mean_model_se = sqrt(as.numeric(coefficient_covariance[[1, 1]])),
  overall_mean_ci_low = as.numeric(model$ci.lb[[1]]),
  overall_mean_ci_high = as.numeric(model$ci.ub[[1]]),
  species_variance = species_variance,
  species_sd = sqrt(species_variance),
  species_variance_profile_ci_low = variance_interval["sigma^2", "ci.lb"],
  species_variance_profile_ci_high = variance_interval["sigma^2", "ci.ub"],
  species_variance_estimable = species_variance_estimable,
  maximum_absolute_standardized_residual = max(abs(standardized_residual)),
  n_absolute_standardized_residual_above_3 = sum(abs(standardized_residual) > 3),
  maximum_approximate_cooks_distance = max(cooks_distance_approx),
  n_approximate_cook_above_4_over_n = sum(
    cooks_distance_approx > 4 / nrow(auc_data)
  ),
  spearman_auc_se_vs_mean_weight = stats::cor(
    auc_data$auc_se, mean_weight_percent, method = "spearman"
  ),
  median_weight_lowest_SE_decile_percent = low_se_median_weight,
  median_weight_highest_SE_decile_percent = high_se_median_weight,
  high_to_low_SE_median_weight_ratio = high_se_median_weight / low_se_median_weight,
  total_weight_highest_SE_decile_percent = sum(
    mean_weight_percent[auc_data$auc_se >= se_90]
  ),
  sampling_covariance_source = sampling$selected_v_type
)

invisible(write_csv_atomic(summary_table, file.path(output_dir, "model_summary.csv")))
invisible(write_csv_atomic(
  diagnostics |>
    dplyr::arrange(dplyr::desc(.data$cooks_distance_approx)),
  file.path(output_dir, "observation_diagnostics.csv")
))
invisible(write_csv_atomic(
  diagnostics |>
    dplyr::slice_max(.data$cooks_distance_approx, n = 20, with_ties = FALSE),
  file.path(output_dir, "most_influential_observations.csv")
))
invisible(write_rds_atomic(
  model, file.path(output_dir, "auc_intercept_species_model.rds")
))

cat("\nDelta AUC ~ 1 with known sampling covariance and species random intercept\n")
cat("Horizon: ", horizon, "\n", sep = "")
cat("Converged without warnings: ", converged, "\n", sep = "")
cat(
  "Overall mean AUC change: ", sprintf("%.6f", stats::coef(model)[[1]]),
  " (model-based 95% CI ", sprintf("%.6f", model$ci.lb[[1]]),
  " to ", sprintf("%.6f", model$ci.ub[[1]]), ")\n", sep = ""
)
cat(
  "Species variance: ", sprintf("%.6f", species_variance),
  " (profile 95% CI ", sprintf("%.6f", variance_interval["sigma^2", "ci.lb"]),
  " to ", sprintf("%.6f", variance_interval["sigma^2", "ci.ub"]), ")\n",
  sep = ""
)
cat("Species variance estimable: ", species_variance_estimable, "\n", sep = "")
cat(
  "Extreme observations: ", sum(abs(standardized_residual) > 3),
  " with |standardized residual| > 3; maximum approximate Cook distance = ",
  sprintf("%.4f", max(cooks_distance_approx)), "\n", sep = ""
)
cat(
  "SE-weight relationship: Spearman rho = ",
  sprintf("%.3f", stats::cor(auc_data$auc_se, mean_weight_percent, method = "spearman")),
  "; highest-SE-decile median weight is ",
  sprintf("%.1f%%", 100 * high_se_median_weight / low_se_median_weight),
  " of the lowest-SE-decile median weight.\n", sep = ""
)
cat("Outputs: ", normalizePath(output_dir), "\n", sep = "")

# ---------------------------------------------------------------------------
# Step 2: full continuous climate surface, retaining only a species intercept

climate_formula <- ~ T_z * P_z * W_z
climate_attempt <- auc_capture(
  metafor::rma.mv(
    yi = auc_effect,
    V = V_auc,
    mods = ~ T_z * P_z * W_z,
    random = ~ 1 | species_id,
    data = auc_data,
    method = "REML",
    sparse = TRUE
  )
)
if (is.null(climate_attempt$fit)) {
  stop("Step-2 climate model failed: ", climate_attempt$error, call. = FALSE)
}
model_climate <- climate_attempt$fit
climate_vcov <- stats::vcov(model_climate)
climate_converged <- !nzchar(climate_attempt$error) &&
  !length(climate_attempt$warnings) &&
  all(is.finite(stats::coef(model_climate))) && all(is.finite(climate_vcov))

climate_coefficients <- tibble::tibble(
  term = names(stats::coef(model_climate)),
  estimate = as.numeric(stats::coef(model_climate)),
  standard_error = as.numeric(model_climate$se),
  conf_low = as.numeric(model_climate$ci.lb),
  conf_high = as.numeric(model_climate$ci.ub)
)

climate_species_variance <- as.numeric(model_climate$sigma2[[1]])
species_variance_reduction_percent <- 100 *
  (species_variance - climate_species_variance) / species_variance

climate_standardized_residual <- as.numeric(
  metafor::rstandard.rma.mv(model_climate)$z
)
climate_leverage <- as.numeric(stats::hatvalues(model_climate))
climate_cooks_distance_approx <- climate_standardized_residual^2 * climate_leverage /
  pmax(ncol(model_climate$X) * (1 - climate_leverage), 1e-8)

climate_diagnostics <- auc_data |>
  dplyr::transmute(
    .data$effect_id, species_id = as.character(.data$species_id),
    experiment_id = as.character(.data$experiment_id),
    .data$auc_effect, .data$auc_se,
    intercept_standardized_residual = .env$standardized_residual,
    climate_standardized_residual = .env$climate_standardized_residual,
    absolute_z_change = abs(.env$climate_standardized_residual) -
      abs(.env$standardized_residual),
    climate_leverage = .env$climate_leverage,
    intercept_cooks_distance_approx = .env$cooks_distance_approx,
    climate_cooks_distance_approx = .env$climate_cooks_distance_approx,
    climate_extreme_abs_z_above_3 = abs(.env$climate_standardized_residual) > 3,
    climate_influential_cook_above_4_over_n =
      .env$climate_cooks_distance_approx > 4 / nrow(auc_data)
  )

experiment_diagnostics <- climate_diagnostics |>
  dplyr::summarise(
    n_effects = dplyr::n(),
    intercept_max_abs_z = max(abs(.data$intercept_standardized_residual)),
    climate_max_abs_z = max(abs(.data$climate_standardized_residual)),
    intercept_sum_cook = sum(.data$intercept_cooks_distance_approx),
    climate_sum_cook = sum(.data$climate_cooks_distance_approx),
    .by = "experiment_id"
  ) |>
  dplyr::arrange(dplyr::desc(.data$climate_sum_cook))

intercept_top_effects <- diagnostics |>
  dplyr::slice_max(
    .data$cooks_distance_approx, n = 20L,
    with_ties = FALSE
  ) |>
  dplyr::pull(.data$effect_id)
climate_top_effects <- climate_diagnostics |>
  dplyr::slice_max(
    .data$climate_cooks_distance_approx, n = 20L,
    with_ties = FALSE
  ) |>
  dplyr::pull(.data$effect_id)

intercept_top_experiments <- experiment_diagnostics |>
  dplyr::slice_max(.data$intercept_sum_cook, n = 10, with_ties = FALSE) |>
  dplyr::pull(.data$experiment_id)
climate_top_experiments <- experiment_diagnostics |>
  dplyr::slice_max(.data$climate_sum_cook, n = 10, with_ties = FALSE) |>
  dplyr::pull(.data$experiment_id)

# Predict only inside the three-dimensional convex hull of the 42 observed
# experiment climates. Facet warming values are the observed 25th, 50th, and
# 75th percentiles across distinct experiments.
support <- build_auc_climate_grid(auc_data, points = 70L)
prediction_grid <- support$grid |>
  dplyr::filter(.data$supported)
climate_model_result <- list(
  fit = model_climate, formula = climate_formula, data = auc_data, V = V_auc,
  random_kind = "species_intercept"
)
climate_predictions <- predict_auc_fixed(
  climate_model_result, list(model = NULL), prediction_grid
)

scaling_value <- function(field, statistic) {
  scaled_auc$scaling[[statistic]][scaled_auc$scaling$scaled_field == field][[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
W_center <- scaling_value("W_z", "center")
W_scale <- scaling_value("W_z", "scale")

warming_labels <- climate_predictions |>
  dplyr::distinct(.data$warming_level, .data$W_z) |>
  dplyr::arrange(.data$warming_level) |>
  dplyr::mutate(
    W_raw = .data$W_z * .env$W_scale + .env$W_center,
    panel_text = paste0(as.character(.data$warming_level), "\nΔT = ",
                        sprintf("%.2f", .data$W_raw), " °C"),
    panel = factor(.data$panel_text, levels = .data$panel_text)
  )
climate_predictions <- climate_predictions |>
  dplyr::left_join(
    dplyr::select(warming_labels, -"panel_text"),
    by = c("warming_level", "W_z")
  )
observed_experiments <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
surface_limit <- max(abs(climate_predictions$estimate), na.rm = TRUE)
zero_contour_predictions <- climate_predictions |>
  dplyr::summarise(
    crosses_zero = min(.data$estimate) <= 0 && max(.data$estimate) >= 0,
    .by = "panel"
  ) |>
  dplyr::filter(.data$crosses_zero) |>
  dplyr::select("panel") |>
  dplyr::inner_join(climate_predictions, by = "panel")

surface_plot <- ggplot2::ggplot(
  climate_predictions,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_raster() +
  ggplot2::geom_contour(
    data = zero_contour_predictions,
    ggplot2::aes(
      x = .data$T_z, y = .data$P_z, z = .data$estimate
    ),
    inherit.aes = FALSE, breaks = 0,
    colour = "black", linewidth = 0.35
  ) +
  ggplot2::geom_point(
    data = observed_experiments,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey20",
    size = 1.4, stroke = 0.35
  ) +
  ggplot2::facet_wrap(~ panel, nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
    limits = c(-surface_limit, surface_limit),
    guide = ggplot2::guide_colorbar(
      barheight = grid::unit(2.5, "in"),
      barwidth = grid::unit(0.18, "in")
    )
  ) +
  ggplot2::scale_x_continuous(
    labels = function(x) sprintf("%.1f", x * T_scale + T_center)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Step 2: predicted resident AUC change (warmed minus origin)",
    subtitle = paste0(
      "Species-intercept model; masked to observed 3D climate support. ",
      "Points show the 42 experiment climates."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Predicted AUC change"
  )

figure_stem <- file.path(
  "output", "frequentist_hierarchical_auc", "figures",
  "step_02_climate_species_intercept_surface"
)
ggplot2::ggsave(
  paste0(figure_stem, ".pdf"), surface_plot,
  width = 11, height = 4.8, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(figure_stem, ".png"), surface_plot,
  width = 11, height = 4.8, dpi = 180, bg = "white"
)

climate_summary <- tibble::tibble(
  horizon = horizon,
  converged_without_warning = climate_converged,
  warnings = paste(climate_attempt$warnings, collapse = " | "),
  error = climate_attempt$error,
  intercept_only_species_variance = species_variance,
  climate_species_variance = climate_species_variance,
  species_variance_reduction_percent = species_variance_reduction_percent,
  intercept_n_abs_z_above_3 = sum(abs(standardized_residual) > 3),
  climate_n_abs_z_above_3 = sum(abs(climate_standardized_residual) > 3),
  intercept_max_abs_z = max(abs(standardized_residual)),
  climate_max_abs_z = max(abs(climate_standardized_residual)),
  intercept_max_cook_approx = max(cooks_distance_approx),
  climate_max_cook_approx = max(climate_cooks_distance_approx),
  top_20_observation_overlap = length(intersect(
    intercept_top_effects, climate_top_effects
  )),
  top_10_experiment_overlap = length(intersect(
    intercept_top_experiments, climate_top_experiments
  )),
  previous_extremes_still_extreme = sum(
    abs(climate_diagnostics$intercept_standardized_residual) > 3 &
      abs(climate_diagnostics$climate_standardized_residual) > 3
  ),
  new_climate_extremes = sum(
    abs(climate_diagnostics$intercept_standardized_residual) <= 3 &
      abs(climate_diagnostics$climate_standardized_residual) > 3
  ),
  median_absolute_z_change_for_previous_extremes = stats::median(
    climate_diagnostics$absolute_z_change[
      abs(climate_diagnostics$intercept_standardized_residual) > 3
    ]
  ),
  support_method = unique(climate_predictions$support_method)[[1]],
  prediction_rows = nrow(climate_predictions),
  figure_pdf = paste0(figure_stem, ".pdf")
)

invisible(write_csv_atomic(
  climate_summary, file.path(output_dir, "climate_model_summary.csv")
))
invisible(write_csv_atomic(
  climate_coefficients, file.path(output_dir, "climate_fixed_effects.csv")
))
invisible(write_csv_atomic(
  climate_diagnostics |>
    dplyr::arrange(dplyr::desc(.data$climate_cooks_distance_approx)),
  file.path(output_dir, "climate_observation_diagnostics.csv")
))
invisible(write_csv_atomic(
  experiment_diagnostics,
  file.path(output_dir, "climate_experiment_diagnostics.csv")
))
invisible(write_csv_atomic(
  climate_predictions, file.path(output_dir, "climate_prediction_grid.csv")
))
invisible(write_rds_atomic(
  model_climate, file.path(output_dir, "auc_climate_species_intercept_model.rds")
))

cat("\nStep 2: Delta AUC ~ T_z * P_z * W_z with species random intercept\n")
cat("Converged without warnings: ", climate_converged, "\n", sep = "")
cat(
  "Species variance: ", sprintf("%.6f", climate_species_variance),
  "; change from intercept-only = ",
  sprintf("%.2f%%", species_variance_reduction_percent), "\n", sep = ""
)
cat(
  "Extreme |z| > 3: ", sum(abs(standardized_residual) > 3), " before; ",
  sum(abs(climate_standardized_residual) > 3), " after. ",
  "Maximum approximate Cook distance: ", sprintf("%.4f", max(cooks_distance_approx)),
  " before; ", sprintf("%.4f", max(climate_cooks_distance_approx)),
  " after.\n", sep = ""
)
cat(
  "Overlap: ", length(intersect(intercept_top_effects, climate_top_effects)),
  "/20 most influential observations and ",
  length(intersect(intercept_top_experiments, climate_top_experiments)),
  "/10 most influential experiments.\n", sep = ""
)
cat("Prediction figure: ", normalizePath(paste0(figure_stem, ".pdf")), "\n", sep = "")

# ---------------------------------------------------------------------------
# Step 3: continuous climate surface with species intercept and warming slope

step3_attempt <- auc_capture(
  metafor::rma.mv(
    yi = auc_effect,
    V = V_auc,
    mods = ~ T_z * P_z * W_z,
    random = list(~ 1 | species_id, ~ 0 + W_z | species_id),
    struct = "GEN",
    data = auc_data,
    method = "REML",
    sparse = TRUE
  )
)
if (is.null(step3_attempt$fit)) {
  stop("Step-3 species warming-slope model failed: ", step3_attempt$error,
       call. = FALSE)
}
model_climate_species_slope <- step3_attempt$fit
step3_vcov <- stats::vcov(model_climate_species_slope)
step3_converged <- !nzchar(step3_attempt$error) &&
  !length(step3_attempt$warnings) &&
  all(is.finite(stats::coef(model_climate_species_slope))) &&
  all(is.finite(step3_vcov))

step3_coefficients <- tibble::tibble(
  term = names(stats::coef(model_climate_species_slope)),
  estimate = as.numeric(stats::coef(model_climate_species_slope)),
  standard_error = as.numeric(model_climate_species_slope$se),
  conf_low = as.numeric(model_climate_species_slope$ci.lb),
  conf_high = as.numeric(model_climate_species_slope$ci.ub)
)

# confint.rma.mv obtains a profile-likelihood interval for the requested
# variance component. A finite upper bound and a lower bound above the numerical
# boundary are required here before calling the slope variance identifiable.
step3_slope_variance <- as.numeric(model_climate_species_slope$tau2[[1]])
step3_intercept_variance <- as.numeric(model_climate_species_slope$sigma2[[1]])
step3_boundary_tolerance <- 1e-8 * max(1, stats::var(auc_data$auc_effect))
step3_slope_ci_attempt <- auc_capture(
  stats::confint(model_climate_species_slope, tau2 = 1)
)
if (is.null(step3_slope_ci_attempt$fit)) {
  step3_slope_ci <- c(estimate = step3_slope_variance, ci.lb = NA_real_,
                      ci.ub = NA_real_)
} else {
  step3_slope_ci <- step3_slope_ci_attempt$fit$random["tau^2", ]
}
step3_slope_identifiable <- step3_slope_variance > step3_boundary_tolerance &&
  all(is.finite(step3_slope_ci[c("ci.lb", "ci.ub")])) &&
  step3_slope_ci[["ci.lb"]] > step3_boundary_tolerance &&
  step3_slope_ci[["ci.ub"]] > step3_slope_variance

# Compare residual and leverage diagnostics directly with Step 2. Both climate
# models have the same eight fixed-effect columns, so their approximate Cook
# distances use the same denominator and are directly comparable.
step3_standardized_residual <- as.numeric(
  metafor::rstandard.rma.mv(model_climate_species_slope)$z
)
step3_leverage <- as.numeric(stats::hatvalues(model_climate_species_slope))
step3_cooks_distance_approx <- step3_standardized_residual^2 * step3_leverage /
  pmax(ncol(model_climate_species_slope$X) * (1 - step3_leverage), 1e-8)

step3_diagnostics <- climate_diagnostics |>
  dplyr::transmute(
    .data$effect_id, .data$species_id, .data$experiment_id,
    .data$auc_effect, .data$auc_se,
    step2_standardized_residual = .data$climate_standardized_residual,
    step3_standardized_residual = .env$step3_standardized_residual,
    absolute_z_change_from_step2 = abs(.env$step3_standardized_residual) -
      abs(.data$climate_standardized_residual),
    step2_cooks_distance_approx = .data$climate_cooks_distance_approx,
    step3_leverage = .env$step3_leverage,
    step3_cooks_distance_approx = .env$step3_cooks_distance_approx,
    step2_extreme_abs_z_above_3 =
      abs(.data$climate_standardized_residual) > 3,
    step3_extreme_abs_z_above_3 = abs(.env$step3_standardized_residual) > 3
  )

step3_experiment_diagnostics <- step3_diagnostics |>
  dplyr::summarise(
    n_effects = dplyr::n(),
    step2_max_abs_z = max(abs(.data$step2_standardized_residual)),
    step3_max_abs_z = max(abs(.data$step3_standardized_residual)),
    step2_sum_cook = sum(.data$step2_cooks_distance_approx),
    step3_sum_cook = sum(.data$step3_cooks_distance_approx),
    .by = "experiment_id"
  ) |>
  dplyr::arrange(dplyr::desc(.data$step3_sum_cook))

step2_top_effects <- step3_diagnostics |>
  dplyr::slice_max(.data$step2_cooks_distance_approx, n = 20L,
                   with_ties = FALSE) |>
  dplyr::pull(.data$effect_id)
step3_top_effects <- step3_diagnostics |>
  dplyr::slice_max(.data$step3_cooks_distance_approx, n = 20L,
                   with_ties = FALSE) |>
  dplyr::pull(.data$effect_id)
step2_top_experiments <- step3_experiment_diagnostics |>
  dplyr::slice_max(.data$step2_sum_cook, n = 10L, with_ties = FALSE) |>
  dplyr::pull(.data$experiment_id)
step3_top_experiments <- step3_experiment_diagnostics |>
  dplyr::slice_max(.data$step3_sum_cook, n = 10L, with_ties = FALSE) |>
  dplyr::pull(.data$experiment_id)

# Fixed-effect prediction surfaces from Steps 2 and 3, evaluated on exactly the
# same supported grid and with a shared colour scale.
step3_model_result <- list(
  fit = model_climate_species_slope, formula = climate_formula,
  data = auc_data, V = V_auc, random_kind = "species_slope"
)
step3_predictions <- predict_auc_fixed(
  step3_model_result, list(model = NULL), prediction_grid
) |>
  dplyr::left_join(
    dplyr::select(warming_labels, -"panel_text"),
    by = c("warming_level", "W_z")
  )

surface_comparison <- dplyr::bind_rows(
  climate_predictions |>
    dplyr::mutate(model = "Step 2: species intercept"),
  step3_predictions |>
    dplyr::mutate(model = "Step 3: intercept + warming slope")
) |>
  dplyr::mutate(
    model = factor(
      .data$model,
      levels = c("Step 2: species intercept",
                 "Step 3: intercept + warming slope")
    )
  )
surface_comparison_limit <- max(abs(surface_comparison$estimate), na.rm = TRUE)
surface_comparison_contours <- surface_comparison |>
  dplyr::summarise(
    crosses_zero = min(.data$estimate) <= 0 && max(.data$estimate) >= 0,
    .by = c("model", "panel")
  ) |>
  dplyr::filter(.data$crosses_zero) |>
  dplyr::select("model", "panel") |>
  dplyr::inner_join(surface_comparison, by = c("model", "panel"))

step3_surface_plot <- ggplot2::ggplot(
  surface_comparison,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_raster() +
  ggplot2::geom_contour(
    data = surface_comparison_contours,
    ggplot2::aes(x = .data$T_z, y = .data$P_z, z = .data$estimate),
    inherit.aes = FALSE, breaks = 0, colour = "black", linewidth = 0.3
  ) +
  ggplot2::geom_point(
    data = observed_experiments,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey20",
    size = 1.1, stroke = 0.3
  ) +
  ggplot2::facet_grid(model ~ panel) +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
    limits = c(-surface_comparison_limit, surface_comparison_limit),
    guide = ggplot2::guide_colorbar(
      barheight = grid::unit(2.8, "in"),
      barwidth = grid::unit(0.18, "in")
    )
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
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Does a species warming slope change the predicted AUC-change surface?",
    subtitle = paste0(
      "Fixed-effect predictions within observed 3D climate support; ",
      "points are the 42 experiment climates."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Predicted AUC change"
  )

step3_surface_stem <- file.path(
  "output", "frequentist_hierarchical_auc", "figures",
  "step_03_climate_surface_comparison"
)
ggplot2::ggsave(
  paste0(step3_surface_stem, ".pdf"), step3_surface_plot,
  width = 11.5, height = 7.2, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(step3_surface_stem, ".png"), step3_surface_plot,
  width = 11.5, height = 7.2, dpi = 180, bg = "white"
)

step3_surface_change <- step3_predictions |>
  dplyr::select(
    "T_z", "P_z", "W_z", "warming_level", "step3_estimate" = "estimate"
  ) |>
  dplyr::left_join(
    climate_predictions |>
      dplyr::select(
        "T_z", "P_z", "W_z", "warming_level",
        "step2_estimate" = "estimate"
      ),
    by = c("T_z", "P_z", "W_z", "warming_level")
  ) |>
  dplyr::mutate(change = .data$step3_estimate - .data$step2_estimate)

# Species-level random slope deviations and total slopes. Because the fixed
# warming derivative depends on T and P, each total slope is evaluated at that
# species' mean observed climate. Dividing by W_scale expresses slopes per °C.
step3_random_effects <- manual_species_random_effects(step3_model_result)
step3_species_support <- auc_data |>
  dplyr::summarise(
    T_z_mean = mean(.data$T_z),
    P_z_mean = mean(.data$P_z),
    n_effects = dplyr::n(),
    n_experiments = dplyr::n_distinct(.data$experiment_id),
    .by = "species_id"
  )
step3_random_slopes <- step3_random_effects |>
  dplyr::filter(.data$component == "species_W_slope") |>
  dplyr::select(
    "species_id", random_slope_deviation_z = "estimate",
    random_slope_standard_error_z = "standard_error",
    random_slope_conf_low_z = "conf_low",
    random_slope_conf_high_z = "conf_high"
  )
step3_beta <- stats::coef(model_climate_species_slope)
step3_species_slopes <- step3_species_support |>
  dplyr::left_join(step3_random_slopes, by = "species_id") |>
  dplyr::mutate(
    fixed_warming_slope_z =
      .env$step3_beta[["W_z"]] +
      .env$step3_beta[["T_z:W_z"]] * .data$T_z_mean +
      .env$step3_beta[["P_z:W_z"]] * .data$P_z_mean +
      .env$step3_beta[["T_z:P_z:W_z"]] * .data$T_z_mean * .data$P_z_mean,
    total_warming_slope_z = .data$fixed_warming_slope_z +
      .data$random_slope_deviation_z,
    random_slope_deviation_per_C = .data$random_slope_deviation_z /
      .env$W_scale,
    total_warming_slope_per_C = .data$total_warming_slope_z / .env$W_scale
  )

step3_slope_plot_data <- dplyr::bind_rows(
  step3_species_slopes |>
    dplyr::transmute(
      .data$species_id,
      panel = "Random species deviation from the climate-dependent slope",
      slope = .data$random_slope_deviation_per_C
    ),
  step3_species_slopes |>
    dplyr::transmute(
      .data$species_id,
      panel = "Total species slope at its mean observed T and P",
      slope = .data$total_warming_slope_per_C
    )
) |>
  dplyr::mutate(
    panel = factor(
      .data$panel,
      levels = c(
        "Random species deviation from the climate-dependent slope",
        "Total species slope at its mean observed T and P"
      )
    )
  )
step3_slope_distribution_summary <- step3_slope_plot_data |>
  dplyr::summarise(
    n_species = dplyr::n(), min = min(.data$slope),
    median = stats::median(.data$slope), max = max(.data$slope),
    .by = "panel"
  )

step3_slope_plot <- ggplot2::ggplot(
  step3_slope_plot_data,
  ggplot2::aes(x = .data$slope)
) +
  ggplot2::geom_histogram(bins = 45, fill = "#0072B2", colour = "white") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.65, colour = "grey20") +
  ggplot2::geom_vline(
    data = step3_slope_distribution_summary,
    ggplot2::aes(xintercept = .data$median),
    colour = "#D55E00", linewidth = 0.7, linetype = 2
  ) +
  ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::labs(
    title = "Distribution of species-specific warming slopes",
    subtitle = paste0(
      "Conditional BLUPs for all ", nrow(step3_species_slopes),
      " species; orange dashed lines are medians."
    ),
    x = "AUC-change slope per 1 °C warming",
    y = "Number of species"
  )

step3_slope_figure_stem <- file.path(
  "output", "frequentist_hierarchical_auc", "figures",
  "step_03_species_warming_slope_distribution"
)
ggplot2::ggsave(
  paste0(step3_slope_figure_stem, ".pdf"), step3_slope_plot,
  width = 8.5, height = 7, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(step3_slope_figure_stem, ".png"), step3_slope_plot,
  width = 8.5, height = 7, dpi = 180, bg = "white"
)

# Overall observed Delta-AUC distribution: one row is one species-experiment
# effect at the selected horizon (not a distribution of fitted values).
observed_auc_summary <- tibble::tibble(
  statistic = factor(c("Minimum", "Median", "Maximum"),
                     levels = c("Minimum", "Median", "Maximum")),
  value = c(min(auc_data$auc_effect), stats::median(auc_data$auc_effect),
            max(auc_data$auc_effect))
)
observed_auc_plot <- ggplot2::ggplot(
  auc_data, ggplot2::aes(x = .data$auc_effect)
) +
  ggplot2::geom_histogram(bins = 80, fill = "#56B4E9", colour = "white") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.7) +
  ggplot2::geom_vline(
    data = observed_auc_summary,
    ggplot2::aes(xintercept = .data$value, colour = .data$statistic),
    linewidth = 0.75, linetype = 2
  ) +
  ggplot2::scale_colour_manual(
    values = c(Minimum = "#0072B2", Median = "#D55E00", Maximum = "#CC79A7")
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::labs(
    title = "Observed resident AUC-change distribution",
    subtitle = paste0(
      "Maximum-duration horizon; ", nrow(auc_data),
      " species–experiment effects. Min = ",
      sprintf("%.4f", min(auc_data$auc_effect)), "; median = ",
      sprintf("%.4f", stats::median(auc_data$auc_effect)), "; max = ",
      sprintf("%.4f", max(auc_data$auc_effect)), "."
    ),
    x = "AUC change (warmed minus origin)", y = "Number of effects",
    colour = "Summary"
  )

observed_auc_figure_stem <- file.path(
  "output", "frequentist_hierarchical_auc", "figures",
  "step_03_observed_auc_change_distribution"
)
ggplot2::ggsave(
  paste0(observed_auc_figure_stem, ".pdf"), observed_auc_plot,
  width = 9, height = 5.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(observed_auc_figure_stem, ".png"), observed_auc_plot,
  width = 9, height = 5.5, dpi = 180, bg = "white"
)

step3_summary <- tibble::tibble(
  horizon = horizon,
  converged_without_warning = step3_converged,
  warnings = paste(step3_attempt$warnings, collapse = " | "),
  error = step3_attempt$error,
  species_intercept_variance = step3_intercept_variance,
  species_warming_slope_variance_z = step3_slope_variance,
  species_warming_slope_variance_ci_low_z = step3_slope_ci[["ci.lb"]],
  species_warming_slope_variance_ci_high_z = step3_slope_ci[["ci.ub"]],
  species_warming_slope_sd_z = sqrt(step3_slope_variance),
  species_warming_slope_sd_per_C = sqrt(step3_slope_variance) / W_scale,
  species_warming_slope_identifiable = step3_slope_identifiable,
  step2_n_abs_z_above_3 = sum(abs(climate_standardized_residual) > 3),
  step3_n_abs_z_above_3 = sum(abs(step3_standardized_residual) > 3),
  step2_max_abs_z = max(abs(climate_standardized_residual)),
  step3_max_abs_z = max(abs(step3_standardized_residual)),
  previous_extremes_still_extreme = sum(
    abs(climate_standardized_residual) > 3 &
      abs(step3_standardized_residual) > 3
  ),
  new_step3_extremes = sum(
    abs(climate_standardized_residual) <= 3 &
      abs(step3_standardized_residual) > 3
  ),
  median_absolute_z_change_for_step2_extremes = stats::median(
    step3_diagnostics$absolute_z_change_from_step2[
      step3_diagnostics$step2_extreme_abs_z_above_3
    ]
  ),
  step2_max_cook_approx = max(climate_cooks_distance_approx),
  step3_max_cook_approx = max(step3_cooks_distance_approx),
  top_20_observation_overlap = length(intersect(
    step2_top_effects, step3_top_effects
  )),
  top_10_experiment_overlap = length(intersect(
    step2_top_experiments, step3_top_experiments
  )),
  surface_prediction_correlation = stats::cor(
    step3_surface_change$step2_estimate, step3_surface_change$step3_estimate
  ),
  surface_prediction_rmse = sqrt(mean(step3_surface_change$change^2)),
  surface_prediction_median_absolute_change = stats::median(
    abs(step3_surface_change$change)
  ),
  surface_prediction_max_absolute_change = max(abs(step3_surface_change$change)),
  observed_auc_change_min = min(auc_data$auc_effect),
  observed_auc_change_median = stats::median(auc_data$auc_effect),
  observed_auc_change_max = max(auc_data$auc_effect),
  surface_figure_pdf = paste0(step3_surface_stem, ".pdf"),
  slope_distribution_figure_pdf = paste0(step3_slope_figure_stem, ".pdf"),
  observed_distribution_figure_pdf = paste0(observed_auc_figure_stem, ".pdf")
)

invisible(write_csv_atomic(
  step3_summary, file.path(output_dir, "step3_model_summary.csv")
))
invisible(write_csv_atomic(
  step3_coefficients, file.path(output_dir, "step3_fixed_effects.csv")
))
invisible(write_csv_atomic(
  step3_diagnostics |>
    dplyr::arrange(dplyr::desc(.data$step3_cooks_distance_approx)),
  file.path(output_dir, "step3_observation_diagnostics.csv")
))
invisible(write_csv_atomic(
  step3_experiment_diagnostics,
  file.path(output_dir, "step3_experiment_diagnostics.csv")
))
invisible(write_csv_atomic(
  step3_predictions, file.path(output_dir, "step3_prediction_grid.csv")
))
invisible(write_csv_atomic(
  step3_surface_change,
  file.path(output_dir, "step3_surface_change_from_step2.csv")
))
invisible(write_csv_atomic(
  step3_species_slopes,
  file.path(output_dir, "step3_species_warming_slopes.csv")
))
invisible(write_csv_atomic(
  step3_slope_distribution_summary,
  file.path(output_dir, "step3_species_slope_distribution_summary.csv")
))
invisible(write_csv_atomic(
  observed_auc_summary,
  file.path(output_dir, "observed_auc_change_distribution_summary.csv")
))
invisible(write_rds_atomic(
  model_climate_species_slope,
  file.path(output_dir, "auc_change_climate_species_slope_model.rds")
))

cat("\nStep 3: Delta AUC ~ T_z * P_z * W_z with species intercept and W slope\n")
cat("Converged without warnings: ", step3_converged, "\n", sep = "")
cat(
  "Species warming-slope variance: ", sprintf("%.7f", step3_slope_variance),
  " (profile 95% CI ", sprintf("%.7f", step3_slope_ci[["ci.lb"]]),
  " to ", sprintf("%.7f", step3_slope_ci[["ci.ub"]]), "); identifiable: ",
  step3_slope_identifiable, "\n", sep = ""
)
cat(
  "Extreme |z| > 3: ", sum(abs(climate_standardized_residual) > 3),
  " in Step 2; ", sum(abs(step3_standardized_residual) > 3),
  " in Step 3.\n", sep = ""
)
cat(
  "Surface correlation with Step 2: ",
  sprintf("%.4f", step3_summary$surface_prediction_correlation),
  "; maximum absolute change = ",
  sprintf("%.5f", step3_summary$surface_prediction_max_absolute_change),
  ".\n", sep = ""
)
cat("Figures:\n", normalizePath(paste0(step3_surface_stem, ".pdf")), "\n",
    normalizePath(paste0(step3_slope_figure_stem, ".pdf")), "\n",
    normalizePath(paste0(observed_auc_figure_stem, ".pdf")), "\n", sep = "")
