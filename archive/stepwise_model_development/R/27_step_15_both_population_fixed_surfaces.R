# Step 15: population fixed-effect surfaces from the saved M_both -----------
#
# This script does not fit or alter a model. It loads the already fitted
# M_both and sets every species random effect to zero. Consequently, the
# predictions describe the hypothetical average species represented by the
# fixed effects.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "here", "metafor", "readr", "tibble",
  "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

model_dir <- file.path("output", "step_by_step_modeling", "species_both")
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
model_file <- file.path(model_dir, "auc_change_species_both_model.rds")
if (!file.exists(model_file)) {
  stop("The saved M_both model is required: ", model_file, call. = FALSE)
}
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load the fitted model and reconstruct only the analysis data/scaling. There
# is deliberately no rma.mv() call in this script.
M_both <- readRDS(model_file)
scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    TP_z = .data$T_z * .data$P_z
  )

expected_fixed_terms <- c(
  "intrcpt", "T_z", "P_z", "W_z", "T_z:P_z", "T_z:W_z", "P_z:W_z",
  "T_z:P_z:W_z"
)
if (!identical(colnames(M_both$X), expected_fixed_terms)) {
  stop("The saved M_both fixed-effect structure is not the expected T*P*W model.",
       call. = FALSE)
}

experiments <- auc_data |>
  dplyr::distinct(
    .data$experiment_id, .data$T_z, .data$P_z, .data$W_z,
    destination_T_JJA = .data$T,
    destination_P_JJA = .data$P_raw,
    warming_C = .data$W
  )
if (nrow(experiments) != 42L || dplyr::n_distinct(experiments$experiment_id) != 42L) {
  stop("Expected exactly 42 distinct experiment climates.", call. = FALSE)
}

scaling_value <- function(field, statistic) {
  value <- scaled_auc$scaling[[statistic]][
    scaled_auc$scaling$scaled_field == field
  ]
  if (length(value) != 1L || !is.finite(value)) {
    stop("Could not recover ", statistic, " for ", field, call. = FALSE)
  }
  value[[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
W_center <- scaling_value("W_z", "center")
W_scale <- scaling_value("W_z", "scale")

# A genuine three-dimensional convex-hull mask is evaluated separately at the
# experiment-weighted 25th, 50th and 75th percentiles of W. This is the same
# mask used by the resident hierarchical pipeline. Unsupported T-P-W cells
# remain in the mask CSV but never receive or display a prediction.
support_result <- build_auc_climate_grid(auc_data, points = 90L)
if (!all(c("supported", "support_method", "warming_level") %in%
         names(support_result$grid))) {
  stop("Three-dimensional support-mask construction failed.", call. = FALSE)
}
support_grid <- support_result$grid |>
  dplyr::mutate(
    destination_T_JJA = .data$T_z * .env$T_scale + .env$T_center,
    destination_P_JJA = exp(.data$P_z * .env$P_scale + .env$P_center),
    warming_C = .data$W_z * .env$W_scale + .env$W_center
  )
if (!identical(unique(support_grid$support_method),
               "three-dimensional convex hull")) {
  stop(
    "This figure requires the three-dimensional convex-hull support mask; ",
    "observed-point fallback was returned: ", support_result$diagnostic,
    call. = FALSE
  )
}

warming_lookup <- support_grid |>
  dplyr::distinct(.data$warming_level, .data$W_z, .data$warming_C) |>
  dplyr::arrange(.data$warming_level) |>
  dplyr::mutate(
    warming_panel = factor(
      paste0(
        dplyr::case_when(
          as.character(.data$warming_level) == "low (25th percentile)" ~ "Low warming",
          as.character(.data$warming_level) == "median" ~ "Median warming",
          as.character(.data$warming_level) == "high (75th percentile)" ~ "High warming",
          TRUE ~ as.character(.data$warming_level)
        ),
        "\nDelta T = ", sprintf("%.2f", .data$warming_C), " °C"
      ),
      levels = paste0(
        c("Low warming", "Median warming", "High warming"),
        "\nDelta T = ", sprintf("%.2f", .data$warming_C), " °C"
      )
    )
  )

supported_grid <- support_grid |>
  dplyr::filter(.data$supported) |>
  dplyr::left_join(warming_lookup, by = c("warming_level", "W_z", "warming_C"))

# Fixed-effect-only predictions: X beta. No Zb term is present, so all species
# random intercept, W, T, P and T*P deviations are exactly zero.
fixed_beta <- stats::coef(M_both)
fixed_vcov <- as.matrix(stats::vcov(M_both))
fixed_design <- auc_design_matrix(
  ~ T_z * P_z * W_z, supported_grid, names(fixed_beta)
)
fixed_estimate <- as.numeric(fixed_design %*% fixed_beta)
fixed_se <- sqrt(pmax(
  0, rowSums((fixed_design %*% fixed_vcov) * fixed_design)
))
fixed_predictions <- supported_grid |>
  dplyr::mutate(
    fixed_predicted_delta_auc = .env$fixed_estimate,
    model_se = .env$fixed_se,
    model_conf_low = .data$fixed_predicted_delta_auc -
      stats::qnorm(0.975) * .data$model_se,
    model_conf_high = .data$fixed_predicted_delta_auc +
      stats::qnorm(0.975) * .data$model_se,
    interpretation = paste(
      "Hypothetical average species; all species random effects set to zero."
    )
  )

# Save the complete mask separately and the actual supported prediction grid.
# This makes it explicit which candidate combinations were excluded.
invisible(write_csv_atomic(
  support_grid,
  file.path(model_dir, "population_fixed_auc_change_support_mask.csv")
))
invisible(write_csv_atomic(
  fixed_predictions,
  file.path(model_dir, "population_fixed_auc_change_prediction_grid.csv")
))

observed_experiments <- experiments |>
  dplyr::select(
    "experiment_id", "T_z", "P_z",
    "destination_T_JJA", "destination_P_JJA"
  )

auc_limit <- max(abs(fixed_predictions$fixed_predicted_delta_auc), na.rm = TRUE)
zero_contour_auc <- fixed_predictions |>
  dplyr::summarise(
    crosses_zero = min(.data$fixed_predicted_delta_auc) <= 0 &&
      max(.data$fixed_predicted_delta_auc) >= 0,
    .by = "warming_panel"
  ) |>
  dplyr::filter(.data$crosses_zero) |>
  dplyr::select("warming_panel") |>
  dplyr::inner_join(fixed_predictions, by = "warming_panel")

raw_T_labels <- function(x) sprintf("%.1f", x * T_scale + T_center)
raw_P_labels <- function(x) sprintf("%.0f", exp(x * P_scale + P_center))

fixed_surface_plot <- ggplot2::ggplot(
  fixed_predictions,
  ggplot2::aes(
    x = .data$T_z, y = .data$P_z,
    fill = .data$fixed_predicted_delta_auc
  )
) +
  ggplot2::geom_raster() +
  ggplot2::geom_contour(
    data = zero_contour_auc,
    ggplot2::aes(
      x = .data$T_z, y = .data$P_z,
      z = .data$fixed_predicted_delta_auc
    ),
    inherit.aes = FALSE, breaks = 0, colour = "grey15", linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data = observed_experiments,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey15",
    size = 1.7, stroke = 0.4
  ) +
  ggplot2::facet_wrap(~ warming_panel, nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-auc_limit, auc_limit)
  ) +
  ggplot2::scale_x_continuous(labels = raw_T_labels) +
  ggplot2::scale_y_continuous(labels = raw_P_labels) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Population fixed-effect Delta-AUC response",
    subtitle = paste(
      "Hypothetical average species; all species random effects are zero.",
      "White points are the 42 experiment climates; predictions use the observed 3D T-P-W support."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Predicted Delta AUC"
  )

surface_stem <- file.path(
  figure_dir, "step_15_both_population_fixed_auc_change_surface"
)
ggplot2::ggsave(
  paste0(surface_stem, ".pdf"), fixed_surface_plot,
  width = 11.5, height = 4.8, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(surface_stem, ".png"), fixed_surface_plot,
  width = 11.5, height = 4.8, dpi = 220, bg = "white"
)

# The fixed-effect derivative is evaluated on the median-W slice of the same
# 3D support. For this trilinear fixed formula it does not depend on W itself:
# d(Delta AUC)/dW_z = beta_W + beta_TW*T_z + beta_PW*P_z +
#                     beta_TPW*T_z*P_z.
# Dividing by the stored W_z scale converts it to Delta AUC per +1 degree C.
median_grid <- supported_grid |>
  dplyr::filter(as.character(.data$warming_level) == "median")
epsilon <- 1e-6
upper <- median_grid |>
  dplyr::mutate(W_z = .data$W_z + .env$epsilon)
lower <- median_grid |>
  dplyr::mutate(W_z = .data$W_z - .env$epsilon)
derivative_design_z <- (
  auc_design_matrix(~ T_z * P_z * W_z, upper, names(fixed_beta)) -
    auc_design_matrix(~ T_z * P_z * W_z, lower, names(fixed_beta))
) / (2 * epsilon)
derivative_design_C <- derivative_design_z / W_scale
warming_slope <- as.numeric(derivative_design_C %*% fixed_beta)
warming_slope_se <- sqrt(pmax(
  0,
  rowSums((derivative_design_C %*% fixed_vcov) * derivative_design_C)
))
warming_slope_grid <- median_grid |>
  dplyr::mutate(
    fixed_warming_slope_per_C = .env$warming_slope,
    model_se_per_C = .env$warming_slope_se,
    model_conf_low_per_C = .data$fixed_warming_slope_per_C -
      stats::qnorm(0.975) * .data$model_se_per_C,
    model_conf_high_per_C = .data$fixed_warming_slope_per_C +
      stats::qnorm(0.975) * .data$model_se_per_C,
    derivative_scale = paste(
      "Delta AUC per +1 degree C experimental warming;",
      "species random W slopes set to zero"
    ),
    support_slice = paste0(
      "Median warming (Delta T = ", sprintf("%.2f", .data$warming_C),
      " degree C) within the observed 3D T-P-W convex hull"
    )
  )
invisible(write_csv_atomic(
  warming_slope_grid,
  file.path(model_dir, "population_fixed_warming_slope_prediction_grid.csv")
))

slope_limit <- max(abs(warming_slope_grid$fixed_warming_slope_per_C), na.rm = TRUE)
zero_contour_slope <- if (
  min(warming_slope_grid$fixed_warming_slope_per_C) <= 0 &&
    max(warming_slope_grid$fixed_warming_slope_per_C) >= 0
) warming_slope_grid else warming_slope_grid[0, ]

warming_slope_plot <- ggplot2::ggplot(
  warming_slope_grid,
  ggplot2::aes(
    x = .data$T_z, y = .data$P_z,
    fill = .data$fixed_warming_slope_per_C
  )
) +
  ggplot2::geom_raster() +
  ggplot2::geom_contour(
    data = zero_contour_slope,
    ggplot2::aes(
      x = .data$T_z, y = .data$P_z,
      z = .data$fixed_warming_slope_per_C
    ),
    inherit.aes = FALSE, breaks = 0, colour = "grey15", linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data = observed_experiments,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey15",
    size = 1.7, stroke = 0.4
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-slope_limit, slope_limit)
  ) +
  ggplot2::scale_x_continuous(labels = raw_T_labels) +
  ggplot2::scale_y_continuous(labels = raw_P_labels) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::labs(
    title = "Population fixed-effect warming slope",
    subtitle = paste(
      "Delta AUC per +1°C additional experimental warming; species random W slopes are zero.\n",
      "Displayed within the median-W slice of the observed 3D climate support."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Delta AUC per +1°C"
  )

slope_stem <- file.path(
  figure_dir, "step_15_both_population_fixed_warming_slope_surface"
)
ggplot2::ggsave(
  paste0(slope_stem, ".pdf"), warming_slope_plot,
  width = 6.8, height = 5.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(slope_stem, ".png"), warming_slope_plot,
  width = 6.8, height = 5.5, dpi = 220, bg = "white"
)

summary_table <- warming_lookup |>
  dplyr::transmute(
    warming_level = as.character(.data$warming_level),
    warming_C = .data$warming_C,
    supported_grid_cells = vapply(
      .data$warming_level,
      function(level) sum(
        support_grid$supported & support_grid$warming_level == level
      ),
      integer(1)
    ),
    total_candidate_grid_cells = vapply(
      .data$warming_level,
      function(level) sum(support_grid$warming_level == level),
      integer(1)
    ),
    support_method = "three-dimensional convex hull"
  )
invisible(write_csv_atomic(
  summary_table,
  file.path(model_dir, "population_fixed_surface_support_summary.csv")
))

cat("Loaded M_both without refitting.\n")
cat("Fixed-effect predictions set every species random effect to zero.\n")
cat("Observed experiment climates overlaid: ", nrow(experiments), "\n", sep = "")
cat("3D support diagnostic: ", support_result$diagnostic, "\n", sep = "")
cat("W_z scaling SD used for per-degree derivative: ",
    sprintf("%.6f", W_scale), " °C\n", sep = "")
cat("AUC-change figure: ", normalizePath(paste0(surface_stem, ".pdf")), "\n",
    sep = "")
cat("Warming-slope figure: ", normalizePath(paste0(slope_stem, ".pdf")), "\n",
    sep = "")
