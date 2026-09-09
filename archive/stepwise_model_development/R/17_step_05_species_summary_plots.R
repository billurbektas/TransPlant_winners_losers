# Species summaries from M_warming and M_context -----------------------------
#
# This script loads the fitted Step-3 and Step-4 models and does not refit them.
# The response is Delta AUC = warmed AUC minus origin AUC.

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

model_dir <- file.path("output", "step_by_step_modeling")
output_dir <- file.path(model_dir, "species_summary_plots")
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

warming_model_file <- file.path(
  model_dir, "auc_intercept_species", "auc_change_climate_species_slope_model.rds"
)
context_model_file <- file.path(
  model_dir, "species_context", "auc_change_species_context_model.rds"
)
if (!file.exists(warming_model_file) || !file.exists(context_model_file)) {
  stop("Run Steps 3 and 4 before creating species summary plots.", call. = FALSE)
}
M_warming <- readRDS(warming_model_file)
M_context <- readRDS(context_model_file)

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    TP_z = .data$T_z * .data$P_z
  )

scaling_value <- function(field, statistic) {
  scaled_auc$scaling[[statistic]][
    scaled_auc$scaling$scaled_field == field
  ][[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
W_scale <- scaling_value("W_z", "scale")

replication_category_warming <- function(n) {
  factor(
    dplyr::case_when(n == 1L ~ "1 experiment",
                     n == 2L ~ "2 experiments",
                     TRUE ~ "≥3 experiments"),
    levels = c("1 experiment", "2 experiments", "≥3 experiments")
  )
}
replication_category_context <- function(n) {
  factor(
    dplyr::case_when(n == 1L ~ "1 experiment",
                     n == 2L ~ "2 experiments",
                     n == 3L ~ "3 experiments",
                     TRUE ~ "≥4 experiments"),
    levels = c("1 experiment", "2 experiments", "3 experiments",
               "≥4 experiments")
  )
}

# A. M_warming ---------------------------------------------------------------
# Conditional species warming deviations and their prediction-error variances.
X_warming <- as.matrix(M_warming$X)
M_warming_inverse_X <- Matrix::solve(M_warming$M, X_warming)
fixed_covariance_warming <- solve(crossprod(X_warming, M_warming_inverse_X))
warming_residual <- as.numeric(M_warming$yi - X_warming %*% M_warming$beta)
warming_weighted_residual <- Matrix::solve(M_warming$M, warming_residual)
warming_species_levels <- levels(M_warming$mf.g$outer)
Z_warming <- Matrix::Diagonal(
  x = as.numeric(M_warming$Z.G1[, 1])
) %*% M_warming$Z.G2
warming_variance <- as.numeric(M_warming$tau2[[1]])
random_warming_deviation <- as.numeric(
  warming_variance * Matrix::crossprod(Z_warming, warming_weighted_residual)
)

M_warming_inverse_Z <- Matrix::solve(M_warming$M, Z_warming)
XtMinvZ <- crossprod(X_warming, M_warming_inverse_Z)
diag_ZtMinvZ <- colSums(as.matrix(Z_warming) * as.matrix(M_warming_inverse_Z))
diag_fixed_adjustment <- colSums(
  as.matrix(XtMinvZ) *
    as.matrix(fixed_covariance_warming %*% as.matrix(XtMinvZ))
)
random_warming_prediction_variance <- pmax(
  0,
  warming_variance - warming_variance^2 *
    (diag_ZtMinvZ - diag_fixed_adjustment)
)
random_warming_prediction_se <- sqrt(random_warming_prediction_variance)
cross_fixed_random_error <-
  -fixed_covariance_warming %*% as.matrix(XtMinvZ) * warming_variance

species_climate_means <- auc_data |>
  dplyr::summarise(
    T_z_species_mean = mean(.data$T_z),
    P_z_species_mean = mean(.data$P_z),
    n_experiments = dplyr::n_distinct(.data$experiment_id),
    .by = "species_id"
  )

derivative_design <- function(T_z, P_z, epsilon = 1e-6) {
  reference <- tibble::tibble(T_z = T_z, P_z = P_z, W_z = 0)
  upper <- dplyr::mutate(reference, W_z = .data$W_z + epsilon)
  lower <- dplyr::mutate(reference, W_z = .data$W_z - epsilon)
  coefficient_names <- names(stats::coef(M_warming))
  (
    auc_design_matrix(~ T_z * P_z * W_z, upper, coefficient_names) -
      auc_design_matrix(~ T_z * P_z * W_z, lower, coefficient_names)
  ) / (2 * epsilon)
}

warming_species <- tibble::tibble(
  species_id = warming_species_levels,
  random_deviation_z = random_warming_deviation,
  random_deviation_se_z = random_warming_prediction_se
) |>
  dplyr::left_join(species_climate_means, by = "species_id") |>
  dplyr::mutate(replication = replication_category_warming(.data$n_experiments))

evaluate_total_warming_slope <- function(T_values, P_values, evaluation) {
  design <- derivative_design(T_values, P_values)
  fixed_slope <- as.numeric(design %*% stats::coef(M_warming))
  fixed_variance <- rowSums(
    (design %*% fixed_covariance_warming) * design
  )
  cross_term <- 2 * rowSums(
    design * t(cross_fixed_random_error)
  )
  total_prediction_variance <- pmax(
    0, fixed_variance + random_warming_prediction_variance + cross_term
  )
  total_slope_z <- fixed_slope + random_warming_deviation
  tibble::tibble(
    species_id = warming_species_levels,
    evaluation = evaluation,
    fixed_slope_z = fixed_slope,
    random_deviation_z = random_warming_deviation,
    random_deviation_se_z = random_warming_prediction_se,
    total_slope_z = total_slope_z,
    total_slope_se_z = sqrt(total_prediction_variance),
    fixed_slope_per_C = fixed_slope / W_scale,
    random_deviation_per_C = random_warming_deviation / W_scale,
    random_deviation_se_per_C = random_warming_prediction_se / W_scale,
    random_deviation_conf_low_per_C =
      (random_warming_deviation - stats::qnorm(0.975) *
         random_warming_prediction_se) / W_scale,
    random_deviation_conf_high_per_C =
      (random_warming_deviation + stats::qnorm(0.975) *
         random_warming_prediction_se) / W_scale,
    total_slope_per_C = total_slope_z / W_scale,
    total_slope_se_per_C = sqrt(total_prediction_variance) / W_scale,
    total_slope_conf_low_per_C =
      (total_slope_z - stats::qnorm(0.975) *
         sqrt(total_prediction_variance)) / W_scale,
    total_slope_conf_high_per_C =
      (total_slope_z + stats::qnorm(0.975) *
         sqrt(total_prediction_variance)) / W_scale
  )
}

warming_slopes <- dplyr::bind_rows(
  evaluate_total_warming_slope(
    rep(0, length(warming_species_levels)),
    rep(0, length(warming_species_levels)),
    "Mean standardized T and P"
  ),
  evaluate_total_warming_slope(
    warming_species$T_z_species_mean,
    warming_species$P_z_species_mean,
    "Each species' mean observed T and P"
  )
) |>
  dplyr::left_join(
    dplyr::select(warming_species, "species_id", "n_experiments", "replication",
                  "T_z_species_mean", "P_z_species_mean"),
    by = "species_id"
  ) |>
  dplyr::mutate(
    evaluation = factor(
      .data$evaluation,
      levels = c("Mean standardized T and P",
                 "Each species' mean observed T and P")
    )
  )

warming_order <- warming_slopes |>
  dplyr::filter(.data$evaluation == "Each species' mean observed T and P") |>
  dplyr::arrange(.data$total_slope_per_C) |>
  dplyr::pull(.data$species_id)
warming_slopes <- warming_slopes |>
  dplyr::mutate(species_ordered = factor(.data$species_id,
                                         levels = warming_order))

warming_colours <- c(
  "1 experiment" = "#999999", "2 experiments" = "#E69F00",
  "≥3 experiments" = "#0072B2"
)
warming_shapes <- c("1 experiment" = 1, "2 experiments" = 17,
                    "≥3 experiments" = 16)

warming_caterpillar <- ggplot2::ggplot(
  warming_slopes,
  ggplot2::aes(y = .data$species_ordered, colour = .data$replication)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.5) +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$total_slope_conf_low_per_C,
                 xend = .data$total_slope_conf_high_per_C,
                 yend = .data$species_ordered),
    linewidth = 0.35, alpha = 0.75
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = .data$total_slope_per_C, shape = .data$replication),
    size = 1.2
  ) +
  ggplot2::facet_wrap(~ evaluation, nrow = 1, scales = "free_x") +
  ggplot2::scale_colour_manual(values = warming_colours, drop = FALSE) +
  ggplot2::scale_shape_manual(values = warming_shapes, drop = FALSE) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 3.8),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  ) +
  ggplot2::labs(
    title = "Species-specific total warming slopes",
    subtitle = paste0(
      "Total slope = climate-dependent fixed warming slope + species BLUP; ",
      "intervals are model-based conditional 95% intervals."
    ),
    x = "Total warming slope (Delta AUC per °C)", y = "Species",
    colour = "Replication", shape = "Replication"
  )

warming_caterpillar_stem <- file.path(
  figure_dir, "step_05_warming_total_slope_caterpillar"
)
ggplot2::ggsave(
  paste0(warming_caterpillar_stem, ".pdf"), warming_caterpillar,
  width = 14, height = 40, device = grDevices::cairo_pdf, limitsize = FALSE
)
ggplot2::ggsave(
  paste0(warming_caterpillar_stem, ".png"), warming_caterpillar,
  width = 14, height = 40, dpi = 120, bg = "white", limitsize = FALSE
)

warming_distribution_data <- dplyr::bind_rows(
  warming_slopes |>
    dplyr::transmute(
      .data$species_id, panel = as.character(.data$evaluation),
      slope = .data$total_slope_per_C
    ),
  warming_species |>
    dplyr::transmute(
      .data$species_id, panel = "Random species deviation only",
      slope = .data$random_deviation_z / W_scale
    )
) |>
  dplyr::mutate(
    panel = factor(
      .data$panel,
      levels = c("Mean standardized T and P",
                 "Each species' mean observed T and P",
                 "Random species deviation only")
    )
  )
warming_distribution_summary <- warming_distribution_data |>
  dplyr::summarise(
    n_species = dplyr::n(), min = min(.data$slope),
    median = stats::median(.data$slope), mean = mean(.data$slope),
    max = max(.data$slope), .by = "panel"
  )

warming_distribution <- ggplot2::ggplot(
  warming_distribution_data, ggplot2::aes(x = .data$slope)
) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ggplot2::after_stat(density)),
    bins = 42, fill = "#56B4E9", colour = "white"
  ) +
  ggplot2::geom_density(colour = "#0072B2", linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.55) +
  ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::labs(
    title = "Distribution of total species warming slopes and BLUP deviations",
    x = "Slope (Delta AUC per °C)", y = "Density"
  )

warming_distribution_stem <- file.path(
  figure_dir, "step_05_warming_slope_distributions"
)
ggplot2::ggsave(
  paste0(warming_distribution_stem, ".pdf"), warming_distribution,
  width = 11, height = 4.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(warming_distribution_stem, ".png"), warming_distribution,
  width = 11, height = 4.5, dpi = 180, bg = "white"
)

# B. M_context ---------------------------------------------------------------
# Obtain conditional random coefficients directly from the fitted metafor
# matrices. The random intercept is independent of the 3x3 GEN context block.
context_residual <- as.numeric(
  M_context$yi - M_context$X %*% M_context$beta
)
context_weighted_residual <- Matrix::solve(M_context$M, context_residual)
context_species_levels <- levels(M_context$mf.g$outer)
n_context_species <- length(context_species_levels)
G_context <- as.matrix(M_context$G)
context_random_intercept <- as.numeric(
  M_context$sigma2[[1]] *
    Matrix::crossprod(M_context$Z.S[[1]], context_weighted_residual)
)
Z_context <- do.call(cbind, lapply(seq_len(ncol(G_context)), function(j) {
  Matrix::Diagonal(x = as.numeric(M_context$Z.G1[, j])) %*% M_context$Z.G2
}))
D_context <- Matrix::kronecker(
  G_context, Matrix::Diagonal(n_context_species)
)
context_random_vector <- as.numeric(
  D_context %*% Matrix::crossprod(Z_context, context_weighted_residual)
)
context_random_matrix <- matrix(
  context_random_vector, nrow = n_context_species, ncol = ncol(G_context),
  dimnames = list(context_species_levels, colnames(G_context))
)
context_random_coefficients <- tibble::tibble(
  species_id = context_species_levels,
  random_intercept = context_random_intercept,
  random_T = context_random_matrix[, "T_z"],
  random_P = context_random_matrix[, "P_z"],
  random_TP = context_random_matrix[, "TP_z"]
)

context_support <- dplyr::bind_rows(lapply(
  dplyr::group_split(auc_data, .data$species_id, .keep = TRUE),
  function(species) {
    points <- dplyr::distinct(species, .data$T_z, .data$P_z)
    design <- cbind(1, points$T_z, points$P_z, points$T_z * points$P_z)
    tibble::tibble(
      species_id = species$species_id[[1]],
      n_experiments = dplyr::n_distinct(species$experiment_id),
      n_distinct_TP = nrow(points),
      context_design_rank = qr(design)$rank,
      full_rank = qr(design)$rank == 4L
    )
  }
)) |>
  dplyr::mutate(replication = replication_category_context(.data$n_experiments))

# A convex hull can be an area, a line, or a point. Generate supported
# evaluation points for all three cases without extrapolating beyond the hull.
species_tp_support_grid <- function(species_name, area_points = 35L,
                                    line_points = 70L) {
  points <- auc_data |>
    dplyr::filter(.data$species_id == species_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  matrix_points <- as.matrix(points[, c("T_z", "P_z")])
  if (nrow(matrix_points) == 1L) {
    return(points)
  }
  centered <- scale(matrix_points, center = TRUE, scale = FALSE)
  dimension <- qr(centered)$rank
  if (dimension < 2L) {
    decomposition <- svd(centered)
    direction <- decomposition$v[, 1]
    scores <- as.numeric(centered %*% direction)
    grid_scores <- seq(min(scores), max(scores), length.out = line_points)
    candidate <- sweep(
      outer(grid_scores, direction), 2,
      attr(centered, "scaled:center"), "+"
    )
    return(tibble::as_tibble(candidate, .name_repair = ~ c("T_z", "P_z")))
  }
  candidate <- tidyr::expand_grid(
    T_z = seq(min(points$T_z), max(points$T_z), length.out = area_points),
    P_z = seq(min(points$P_z), max(points$P_z), length.out = area_points)
  )
  hull <- geometry::convhulln(matrix_points, options = "Qt")
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    ))
}

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
warming_quantiles <- stats::quantile(
  experiment_warming$W_z, c(0.25, 0.50, 0.75), names = FALSE
)
warming_raw_quantiles <- warming_quantiles * W_scale +
  scaling_value("W_z", "center")
warming_levels <- tibble::tibble(
  warming_level = factor(c("low", "median", "high"),
                         levels = c("low", "median", "high")),
  W_z = warming_quantiles,
  warming_C = warming_raw_quantiles
)
context_beta <- stats::coef(M_context)

context_prediction_one <- function(species_name) {
  support <- species_tp_support_grid(species_name)
  candidate <- tidyr::crossing(support, warming_levels) |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  design <- auc_design_matrix(
    ~ T_z * P_z * W_z, candidate, names(context_beta)
  )
  random <- context_random_coefficients |>
    dplyr::filter(.data$species_id == species_name)
  candidate |>
    dplyr::mutate(
      species_id = species_name,
      estimate = as.numeric(design %*% context_beta) +
        random$random_intercept[[1]] + random$random_T[[1]] * .data$T_z +
        random$random_P[[1]] * .data$P_z + random$random_TP[[1]] * .data$TP_z
    )
}

context_predictions <- dplyr::bind_rows(lapply(
  context_species_levels, context_prediction_one
))
context_interval_summary <- context_predictions |>
  dplyr::summarise(
    min_predicted_auc_change = min(.data$estimate),
    max_predicted_auc_change = max(.data$estimate),
    mean_predicted_auc_change = mean(.data$estimate),
    crosses_zero = .data$min_predicted_auc_change <= 0 &
      .data$max_predicted_auc_change >= 0,
    candidate_pattern = dplyr::case_when(
      .data$min_predicted_auc_change > 0 ~ "consistent winner candidate",
      .data$max_predicted_auc_change < 0 ~ "consistent loser candidate",
      TRUE ~ "context-dependent switcher candidate"
    ),
    n_grid_points = dplyr::n(),
    .by = c("species_id", "warming_level", "warming_C")
  ) |>
  dplyr::left_join(context_support, by = "species_id")

median_context <- context_interval_summary |>
  dplyr::filter(.data$warming_level == "median")
context_order <- median_context |>
  dplyr::arrange(.data$mean_predicted_auc_change) |>
  dplyr::pull(.data$species_id)
median_context <- median_context |>
  dplyr::mutate(
    species_ordered = factor(.data$species_id, levels = context_order),
    rank_support = factor(
      dplyr::if_else(.data$full_rank, "Full-rank T-P", "Partial-rank T-P"),
      levels = c("Partial-rank T-P", "Full-rank T-P")
    )
  )

context_colours <- c(
  "1 experiment" = "#999999", "2 experiments" = "#E69F00",
  "3 experiments" = "#56B4E9", "≥4 experiments" = "#0072B2"
)
rank_shapes <- c("Partial-rank T-P" = 1, "Full-rank T-P" = 16)

make_context_interval_plot <- function(data, title, text_size = 4) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(y = .data$species_ordered, colour = .data$replication)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.55) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$min_predicted_auc_change,
                   xend = .data$max_predicted_auc_change,
                   yend = .data$species_ordered),
      linewidth = 0.6, alpha = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$mean_predicted_auc_change,
                   shape = .data$rank_support),
      size = 1.5
    ) +
    ggplot2::scale_colour_manual(values = context_colours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = rank_shapes, drop = FALSE) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = text_size),
      legend.position = "top"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Line = min–max and point = mean predicted Delta AUC across each ",
        "species' observed T-P support at median warming (Delta T = ",
        sprintf("%.2f", warming_raw_quantiles[[2]]), " °C)."
      ),
      x = "Predicted Delta AUC", y = "Species",
      colour = "Replication", shape = "T-P support"
    )
}

context_all_plot <- make_context_interval_plot(
  median_context, "Species background-climate prediction ranges: all species",
  text_size = 3.8
)
context_all_stem <- file.path(
  figure_dir, "step_05_context_interval_all_species"
)
ggplot2::ggsave(
  paste0(context_all_stem, ".pdf"), context_all_plot,
  width = 11, height = 40, device = grDevices::cairo_pdf, limitsize = FALSE
)
ggplot2::ggsave(
  paste0(context_all_stem, ".png"), context_all_plot,
  width = 11, height = 40, dpi = 120, bg = "white", limitsize = FALSE
)

better_supported_context <- median_context |>
  dplyr::filter(.data$n_experiments >= 4L) |>
  dplyr::mutate(
    species_ordered = factor(.data$species_id,
                             levels = intersect(context_order, .data$species_id))
  )
context_better_plot <- make_context_interval_plot(
  better_supported_context,
  "Species background-climate prediction ranges: at least four experiments",
  text_size = 5.5
)
context_better_stem <- file.path(
  figure_dir, "step_05_context_interval_better_supported"
)
better_height <- max(12, 0.115 * nrow(better_supported_context))
ggplot2::ggsave(
  paste0(context_better_stem, ".pdf"), context_better_plot,
  width = 11, height = better_height, device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggplot2::ggsave(
  paste0(context_better_stem, ".png"), context_better_plot,
  width = 11, height = better_height, dpi = 150, bg = "white", limitsize = FALSE
)

# Low/median/high sensitivity for the better-supported species, preserving the
# median-warming order so changes across panels are directly comparable.
context_sensitivity <- context_interval_summary |>
  dplyr::filter(.data$n_experiments >= 4L) |>
  dplyr::mutate(
    species_ordered = factor(.data$species_id,
                             levels = intersect(context_order, .data$species_id)),
    rank_support = factor(
      dplyr::if_else(.data$full_rank, "Full-rank T-P", "Partial-rank T-P"),
      levels = c("Partial-rank T-P", "Full-rank T-P")
    ),
    warming_panel = factor(
      paste0(as.character(.data$warming_level), " (Delta T = ",
             sprintf("%.2f", .data$warming_C), " °C)"),
      levels = paste0(c("low", "median", "high"), " (Delta T = ",
                      sprintf("%.2f", warming_raw_quantiles), " °C)")
    )
  )
context_sensitivity_plot <- ggplot2::ggplot(
  context_sensitivity,
  ggplot2::aes(y = .data$species_ordered, colour = .data$replication)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.5) +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$min_predicted_auc_change,
                 xend = .data$max_predicted_auc_change,
                 yend = .data$species_ordered),
    linewidth = 0.5, alpha = 0.8
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = .data$mean_predicted_auc_change,
                 shape = .data$rank_support), size = 1.2
  ) +
  ggplot2::facet_wrap(~ warming_panel, nrow = 1, scales = "free_x") +
  ggplot2::scale_colour_manual(values = context_colours, drop = FALSE) +
  ggplot2::scale_shape_manual(values = rank_shapes, drop = FALSE) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 4.6),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  ) +
  ggplot2::labs(
    title = "Background-climate prediction ranges across warming levels",
    subtitle = "Species with at least four experiments; order fixed by median-warming mean.",
    x = "Predicted Delta AUC", y = "Species",
    colour = "Replication", shape = "T-P support"
  )
context_sensitivity_stem <- file.path(
  figure_dir, "step_05_context_interval_warming_sensitivity"
)
ggplot2::ggsave(
  paste0(context_sensitivity_stem, ".pdf"), context_sensitivity_plot,
  width = 15, height = better_height, device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggplot2::ggsave(
  paste0(context_sensitivity_stem, ".png"), context_sensitivity_plot,
  width = 15, height = better_height, dpi = 140, bg = "white", limitsize = FALSE
)

context_pattern_summary <- context_interval_summary |>
  dplyr::summarise(
    n_species = dplyr::n(),
    .by = c("warming_level", "warming_C", "candidate_pattern")
  ) |>
  dplyr::arrange(.data$warming_level, .data$candidate_pattern)

invisible(write_csv_atomic(
  warming_slopes, file.path(output_dir, "warming_species_total_slopes.csv")
))
invisible(write_csv_atomic(
  warming_distribution_summary,
  file.path(output_dir, "warming_slope_distribution_summary.csv")
))
invisible(write_csv_atomic(
  context_support, file.path(output_dir, "context_species_support.csv")
))
invisible(write_csv_atomic(
  context_interval_summary,
  file.path(output_dir, "context_species_interval_summaries.csv")
))
invisible(write_csv_atomic(
  context_pattern_summary,
  file.path(output_dir, "context_candidate_pattern_counts.csv")
))
invisible(write_csv_atomic(
  context_predictions,
  file.path(output_dir, "context_species_supported_predictions.csv")
))

cat("Species summary plots completed without refitting M_warming or M_context.\n")
cat("Warming caterpillar: ",
    normalizePath(paste0(warming_caterpillar_stem, ".pdf")), "\n", sep = "")
cat("Warming distributions: ",
    normalizePath(paste0(warming_distribution_stem, ".pdf")), "\n", sep = "")
cat("Context all species: ",
    normalizePath(paste0(context_all_stem, ".pdf")), "\n", sep = "")
cat("Context better supported: ",
    normalizePath(paste0(context_better_stem, ".pdf")), "\n", sep = "")
