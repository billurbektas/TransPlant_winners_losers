# Step 10: reproduce M_context classification figures using M_both -----------
#
# Loads the already fitted M_both. No model is fitted in this script.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "ggrepel", "Matrix", "metafor",
  "patchwork", "readr", "scales", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

both_dir <- file.path("output", "step_by_step_modeling", "species_both")
context_dir <- file.path("output", "step_by_step_modeling", "species_context")
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
both_model_file <- file.path(both_dir, "auc_change_species_both_model.rds")
context_model_file <- file.path(context_dir, "auc_change_species_context_model.rds")
context_classification_file <- file.path(
  context_dir, "context_full_rank_minmax_classification.csv"
)
context_random_file <- file.path(
  context_dir, "species_random_context_coefficients.csv"
)
selection_file <- file.path(context_dir, "context_step08_selected_species.csv")
required_files <- c(
  both_model_file, context_model_file, context_classification_file,
  context_random_file, selection_file
)
if (any(!file.exists(required_files))) {
  stop("Missing saved prerequisite: ",
       paste(required_files[!file.exists(required_files)], collapse = ", "),
       call. = FALSE)
}
dir.create(both_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

M_both <- readRDS(both_model_file)
M_context <- readRDS(context_model_file)
context_classification <- readr::read_csv(
  context_classification_file, show_col_types = FALSE
) |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    fitted_class = factor(
      .data$fitted_class,
      levels = c("Consistent loser", "Context-dependent", "Consistent winner")
    )
  )
if (nrow(context_classification) != 118L ||
    any(context_classification$n_experiments < 4L) ||
    any(context_classification$context_design_rank != 4L)) {
  stop("The saved M_context eligible set is not the required 118 species.",
       call. = FALSE)
}

surface_selection <- readr::read_csv(selection_file, show_col_types = FALSE) |>
  dplyr::arrange(.data$display_order)
selected_species <- surface_selection$species_id
selected_panel_levels <- paste0(
  surface_selection$species_id, "\n(", surface_selection$n_experiments,
  " experiments)"
)

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    TP_z = .data$T_z * .data$P_z
  )

# Conditional species coefficients for the five independent/correlated random
# components in M_both: intercept, W, and the T/P/T*P block.
X_model <- as.matrix(M_both$X)
fixed_beta <- stats::coef(M_both)
model_residual <- as.numeric(M_both$yi - X_model %*% M_both$beta)
weighted_residual <- Matrix::solve(M_both$M, model_residual)
species_levels <- levels(M_both$mf.g$outer)
context_species_levels <- levels(M_both$mf.h$outer)
if (!identical(species_levels, context_species_levels) ||
    !identical(colnames(M_both$G), "W_z") ||
    !identical(colnames(M_both$H), c("T_z", "P_z", "TP_z"))) {
  stop("M_both random-effect blocks or species ordering are unexpected.",
       call. = FALSE)
}
n_species <- length(species_levels)
Z_intercept <- M_both$Z.S[[1]]
Z_warming <- Matrix::Diagonal(x = as.numeric(M_both$Z.G1[, 1])) %*%
  M_both$Z.G2
Z_context <- do.call(cbind, lapply(seq_len(ncol(M_both$H)), function(j) {
  Matrix::Diagonal(x = as.numeric(M_both$Z.H1[, j])) %*% M_both$Z.H2
}))

random_intercept <- as.numeric(
  M_both$sigma2[[1]] * Matrix::crossprod(Z_intercept, weighted_residual)
)
random_W <- as.numeric(
  M_both$tau2[[1]] * Matrix::crossprod(Z_warming, weighted_residual)
)
D_context <- Matrix::kronecker(
  as.matrix(M_both$H), Matrix::Diagonal(n_species)
)
random_context_vector <- as.numeric(
  D_context %*% Matrix::crossprod(Z_context, weighted_residual)
)
random_context_matrix <- matrix(
  random_context_vector, nrow = n_species, ncol = 3L,
  dimnames = list(species_levels, c("T_z", "P_z", "TP_z"))
)
both_random_coefficients <- tibble::tibble(
  species_id = species_levels,
  random_intercept = random_intercept,
  random_W = random_W,
  random_T = random_context_matrix[, "T_z"],
  random_P = random_context_matrix[, "P_z"],
  random_TP = random_context_matrix[, "TP_z"]
)
invisible(write_csv_atomic(
  both_random_coefficients,
  file.path(both_dir, "species_random_both_coefficients.csv")
))

# Conditional prediction-error covariance for total species fitted means.
M_inverse_X <- Matrix::solve(M_both$M, X_model)
fixed_covariance <- solve(crossprod(X_model, M_inverse_X))
D_species <- matrix(
  0, nrow = 5L, ncol = 5L,
  dimnames = list(
    c("intercept", "W_z", "T_z", "P_z", "TP_z"),
    c("intercept", "W_z", "T_z", "P_z", "TP_z")
  )
)
D_species["intercept", "intercept"] <- M_both$sigma2[[1]]
D_species["W_z", "W_z"] <- M_both$tau2[[1]]
D_species[c("T_z", "P_z", "TP_z"),
          c("T_z", "P_z", "TP_z")] <- as.matrix(M_both$H)

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
W_z_median <- stats::median(experiment_warming$W_z)
W_median <- stats::median(experiment_warming$W)

species_hull_grid <- function(species_name, points = 81L) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == species_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  )
  hull <- geometry::convhulln(
    as.matrix(observed[, c("T_z", "P_z")]), options = "Qt"
  )
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::mutate(support_point_type = "regular hull grid") |>
    dplyr::bind_rows(
      dplyr::mutate(observed, support_point_type = "observed hull vertex")
    ) |>
    dplyr::distinct(.data$T_z, .data$P_z, .keep_all = TRUE) |>
    dplyr::mutate(W_z = W_z_median, TP_z = .data$T_z * .data$P_z)
}

prediction_one_species <- function(species_name) {
  species_index <- match(species_name, species_levels)
  grid <- species_hull_grid(species_name)
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, grid, names(fixed_beta)
  )
  random_design <- cbind(
    intercept = 1, W_z = grid$W_z, T_z = grid$T_z,
    P_z = grid$P_z, TP_z = grid$TP_z
  )
  random <- both_random_coefficients |>
    dplyr::filter(.data$species_id == species_name)
  random_vector <- unname(as.numeric(random[
    1, c("random_intercept", "random_W", "random_T", "random_P", "random_TP")
  ]))
  Z_species <- cbind(
    Z_intercept[, species_index, drop = FALSE],
    Z_warming[, species_index, drop = FALSE],
    Z_context[, species_index, drop = FALSE],
    Z_context[, n_species + species_index, drop = FALSE],
    Z_context[, 2L * n_species + species_index, drop = FALSE]
  )
  M_inverse_Z <- Matrix::solve(M_both$M, Z_species)
  XtM_inverse_Z <- as.matrix(crossprod(X_model, M_inverse_Z))
  ZtM_inverse_Z <- as.matrix(crossprod(Z_species, M_inverse_Z))
  random_prediction_error_covariance <-
    D_species - D_species %*% ZtM_inverse_Z %*% D_species +
    D_species %*% t(XtM_inverse_Z) %*% fixed_covariance %*%
      XtM_inverse_Z %*% D_species
  fixed_random_error_covariance <-
    -fixed_covariance %*% XtM_inverse_Z %*% D_species
  fixed_variance <- rowSums(
    (fixed_design %*% fixed_covariance) * fixed_design
  )
  random_variance <- rowSums(
    (random_design %*% random_prediction_error_covariance) * random_design
  )
  cross_covariance <- 2 * rowSums(
    (fixed_design %*% fixed_random_error_covariance) * random_design
  )
  prediction_variance <- pmax(
    0, fixed_variance + random_variance + cross_covariance
  )
  estimate <- as.numeric(fixed_design %*% fixed_beta) +
    as.numeric(random_design %*% random_vector)
  grid |>
    dplyr::mutate(
      species_id = species_name,
      estimate = estimate,
      prediction_se = sqrt(prediction_variance),
      conf_low = .data$estimate - stats::qnorm(0.975) * .data$prediction_se,
      conf_high = .data$estimate + stats::qnorm(0.975) * .data$prediction_se
    )
}

both_predictions <- dplyr::bind_rows(lapply(
  context_classification$species_id, prediction_one_species
))
both_classification <- both_predictions |>
  dplyr::group_split(.data$species_id, .keep = TRUE) |>
  lapply(function(species) {
    minimum <- species[which.min(species$estimate), , drop = FALSE]
    maximum <- species[which.max(species$estimate), , drop = FALSE]
    tibble::tibble(
      species_id = species$species_id[[1]],
      min_fitted_delta_auc = minimum$estimate[[1]],
      min_prediction_se = minimum$prediction_se[[1]],
      min_conf_low = minimum$conf_low[[1]],
      min_conf_high = minimum$conf_high[[1]],
      min_T_z = minimum$T_z[[1]], min_P_z = minimum$P_z[[1]],
      max_fitted_delta_auc = maximum$estimate[[1]],
      max_prediction_se = maximum$prediction_se[[1]],
      max_conf_low = maximum$conf_low[[1]],
      max_conf_high = maximum$conf_high[[1]],
      max_T_z = maximum$T_z[[1]], max_P_z = maximum$P_z[[1]],
      n_support_grid_points = nrow(species)
    )
  }) |>
  dplyr::bind_rows() |>
  dplyr::left_join(
    dplyr::select(
      context_classification, "species_id", "n_experiments",
      "n_distinct_TP", "context_design_rank"
    ), by = "species_id"
  ) |>
  dplyr::mutate(
    fitted_class = dplyr::case_when(
      .data$min_fitted_delta_auc > 0 ~ "Consistent winner",
      .data$max_fitted_delta_auc < 0 ~ "Consistent loser",
      TRUE ~ "Context-dependent"
    ),
    fitted_class = factor(
      .data$fitted_class,
      levels = c("Consistent loser", "Context-dependent", "Consistent winner")
    ),
    warming_C = W_median,
    warming_z = W_z_median
  )
classification_valid <- with(
  both_classification,
  (fitted_class == "Consistent winner" & min_fitted_delta_auc > 0) |
    (fitted_class == "Consistent loser" & max_fitted_delta_auc < 0) |
    (fitted_class == "Context-dependent" &
       min_fitted_delta_auc <= 0 & max_fitted_delta_auc >= 0)
)
if (!all(classification_valid) || nrow(both_classification) != 118L) {
  stop("M_both classification failed its rule/set validation.", call. = FALSE)
}
invisible(write_csv_atomic(
  both_classification,
  file.path(both_dir, "both_full_rank_minmax_classification.csv")
))

classification_comparison <- context_classification |>
  dplyr::transmute(
    species = .data$species_id,
    `M_context classification` = as.character(.data$fitted_class)
  ) |>
  dplyr::left_join(
    both_classification |>
      dplyr::transmute(
        species = .data$species_id,
        `M_both classification` = as.character(.data$fitted_class)
      ),
    by = "species"
  ) |>
  dplyr::mutate(
    `changed? yes/no` = dplyr::if_else(
      .data$`M_context classification` == .data$`M_both classification`,
      "no", "yes"
    )
  )
invisible(write_csv_atomic(
  classification_comparison,
  file.path(both_dir, "context_vs_both_species_classification.csv")
))

class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#5E3C99",
  "Consistent winner" = "#1B9E77"
)
context_label_ids <- dplyr::bind_rows(
  context_classification |>
    dplyr::filter(.data$fitted_class == "Consistent winner") |>
    dplyr::slice_max(.data$min_fitted_delta_auc, n = 3L, with_ties = FALSE),
  context_classification |>
    dplyr::filter(.data$fitted_class == "Consistent loser") |>
    dplyr::slice_min(.data$max_fitted_delta_auc, n = 3L, with_ties = FALSE),
  context_classification |>
    dplyr::filter(.data$fitted_class == "Context-dependent") |>
    dplyr::mutate(response_span = .data$max_fitted_delta_auc -
                    .data$min_fitted_delta_auc) |>
    dplyr::slice_max(.data$response_span, n = 4L, with_ties = FALSE)
) |>
  dplyr::distinct(.data$species_id) |>
  dplyr::pull(.data$species_id)

make_classification_plot <- function(label_ids, subtitle) {
  labels <- both_classification |>
    dplyr::filter(.data$species_id %in% label_ids)
  x_range <- range(c(
    both_classification$min_conf_low, both_classification$min_conf_high,
    both_classification$min_fitted_delta_auc, 0
  ), na.rm = TRUE)
  y_range <- range(c(
    both_classification$max_conf_low, both_classification$max_conf_high,
    both_classification$max_fitted_delta_auc, 0
  ), na.rm = TRUE)
  x_limits <- x_range + c(-0.08, 0.08) * diff(x_range)
  y_limits <- y_range + c(-0.08, 0.08) * diff(y_range)
  ggplot2::ggplot(both_classification) +
    ggplot2::annotate(
      "rect", xmin = x_limits[[1]], xmax = 0,
      ymin = y_limits[[1]], ymax = 0,
      fill = class_colours[["Consistent loser"]], alpha = 0.07
    ) +
    ggplot2::annotate(
      "rect", xmin = x_limits[[1]], xmax = 0,
      ymin = 0, ymax = y_limits[[2]],
      fill = class_colours[["Context-dependent"]], alpha = 0.07
    ) +
    ggplot2::annotate(
      "rect", xmin = 0, xmax = x_limits[[2]],
      ymin = 0, ymax = y_limits[[2]],
      fill = class_colours[["Consistent winner"]], alpha = 0.07
    ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      colour = "grey55", linewidth = 0.4
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$min_conf_low, xend = .data$min_conf_high,
        y = .data$max_fitted_delta_auc, yend = .data$max_fitted_delta_auc,
        colour = .data$fitted_class
      ), linewidth = 0.4, alpha = 0.22
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$min_fitted_delta_auc, xend = .data$min_fitted_delta_auc,
        y = .data$max_conf_low, yend = .data$max_conf_high,
        colour = .data$fitted_class
      ), linewidth = 0.4, alpha = 0.22
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$min_fitted_delta_auc,
        y = .data$max_fitted_delta_auc,
        colour = .data$fitted_class,
        fill = .data$n_experiments
      ), shape = 21, size = 3, stroke = 0.85, alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(
        x = .data$min_fitted_delta_auc,
        y = .data$max_fitted_delta_auc,
        label = .data$species_id,
        colour = .data$fitted_class
      ), size = 3, show.legend = FALSE, min.segment.length = 0,
      segment.colour = "grey60", segment.size = 0.25,
      box.padding = 0.35, point.padding = 0.25, seed = 20260907
    ) +
    ggplot2::scale_colour_manual(values = class_colours, drop = FALSE) +
    ggplot2::scale_fill_viridis_c(
      option = "C", begin = 0.12, end = 0.88,
      name = "Number of\nexperiments"
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(), legend.position = "top",
      plot.margin = ggplot2::margin(6, 15, 6, 6)
    ) +
    ggplot2::labs(
      title = "M_both: minimum versus maximum fitted Delta AUC",
      subtitle = subtitle,
      x = "Minimum fitted Delta AUC across observed T-P support",
      y = "Maximum fitted Delta AUC across observed T-P support",
      colour = "Fitted classification",
      caption = paste0(
        "Pointwise model-based 95% intervals at the fitted extrema; not ",
        "simultaneous intervals for the whole surface."
      )
    )
}

classification_plot <- make_classification_plot(
  context_label_ids,
  paste0(
    "Same 118 species and M_context label set; median warming (Delta T = ",
    sprintf("%.2f", W_median), " °C)."
  )
)
classification_selected_plot <- make_classification_plot(
  selected_species,
  "Only species displayed in the matched surface grid are labelled."
)
classification_stem <- file.path(
  figure_dir, "step_09_both_context_minmax_classification"
)
ggplot2::ggsave(
  paste0(classification_stem, ".pdf"), classification_plot,
  width = 11, height = 9, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(classification_stem, ".png"), classification_plot,
  width = 11, height = 9, dpi = 220, bg = "white"
)

# Matched median/high selected-species surfaces. Recompute the saved M_context
# surface range so M_both uses exactly the same color limits.
warming_z <- stats::quantile(
  experiment_warming$W_z, c(0.50, 0.75), names = FALSE
)
warming_C <- stats::quantile(
  experiment_warming$W, c(0.50, 0.75), names = FALSE
)
warming_levels <- tibble::tibble(
  warming_panel = factor(
    paste0(c("Median warming", "High warming"), "\nDelta T = ",
           sprintf("%.2f", warming_C), " °C"),
    levels = paste0(c("Median warming", "High warming"), "\nDelta T = ",
                    sprintf("%.2f", warming_C), " °C")
  ),
  W_z = warming_z
)
context_random <- readr::read_csv(
  context_random_file, show_col_types = FALSE
)

surface_grid_one <- function(species_name, grid_points = 71L) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == species_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z),
              length.out = grid_points),
    P_z = seq(min(observed$P_z), max(observed$P_z),
              length.out = grid_points)
  )
  hull <- geometry::convhulln(
    as.matrix(observed[, c("T_z", "P_z")]), options = "Qt"
  )
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    tidyr::crossing(warming_levels) |>
    dplyr::mutate(
      species_id = species_name,
      TP_z = .data$T_z * .data$P_z,
      tile_width = diff(range(observed$T_z)) / (grid_points - 1L) * 1.02,
      tile_height = diff(range(observed$P_z)) / (grid_points - 1L) * 1.02
    )
}
surface_grid <- dplyr::bind_rows(lapply(selected_species, surface_grid_one)) |>
  dplyr::left_join(
    dplyr::transmute(
      surface_selection, species_id = .data$species_id,
      species_panel = paste0(.data$species_id, "\n(", .data$n_experiments,
                             " experiments)")
    ), by = "species_id"
  ) |>
  dplyr::mutate(
    species_panel = factor(.data$species_panel,
                           levels = selected_panel_levels)
  )

predict_surfaces <- function(model, random_data, include_W) {
  design <- auc_design_matrix(
    ~ T_z * P_z * W_z, surface_grid, names(stats::coef(model))
  )
  joined <- surface_grid |>
    dplyr::left_join(random_data, by = "species_id")
  random_part <- joined$random_intercept + joined$random_T * joined$T_z +
    joined$random_P * joined$P_z + joined$random_TP * joined$TP_z
  if (include_W) {
    random_part <- random_part + joined$random_W * joined$W_z
  }
  joined |>
    dplyr::mutate(
      estimate = as.numeric(design %*% stats::coef(model)) + random_part
    )
}
context_surfaces <- predict_surfaces(M_context, context_random, FALSE)
both_surfaces <- predict_surfaces(M_both, both_random_coefficients, TRUE)
context_surface_limit <- max(abs(context_surfaces$estimate), na.rm = TRUE)

surface_observations <- auc_data |>
  dplyr::filter(.data$species_id %in% selected_species) |>
  dplyr::distinct(.data$species_id, .data$experiment_id, .data$T_z, .data$P_z) |>
  dplyr::left_join(
    dplyr::transmute(
      surface_selection, species_id = .data$species_id,
      species_panel = paste0(.data$species_id, "\n(", .data$n_experiments,
                             " experiments)")
    ), by = "species_id"
  ) |>
  tidyr::crossing(dplyr::select(warming_levels, "warming_panel")) |>
  dplyr::mutate(
    species_panel = factor(.data$species_panel,
                           levels = selected_panel_levels)
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

surface_plot <- ggplot2::ggplot(
  both_surfaces,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_tile(
    ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
  ) +
  ggplot2::geom_point(
    data = surface_observations,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey20",
    size = 0.9, stroke = 0.3
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$warming_panel),
    cols = ggplot2::vars(.data$species_panel),
    scales = "fixed", space = "fixed"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-context_surface_limit, context_surface_limit),
    oob = scales::squish
  ) +
  ggplot2::scale_x_continuous(
    labels = function(x) sprintf("%.1f", x * T_scale + T_center)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(face = "bold", size = 6.8),
    strip.text.y = ggplot2::element_text(face = "bold", size = 7),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "M_both: selected full-rank species surfaces",
    subtitle = paste0(
      "Same species, support masks, fixed T-P axes, layout, and M_context ",
      "color limits."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Fitted Delta AUC"
  )

surface_stem <- file.path(
  figure_dir, "step_09_both_selected_species_tp_surfaces_median_high"
)
combined_stem <- file.path(
  figure_dir, "step_09_both_classification_selected_surfaces_median_high"
)
combined_plot <- patchwork::wrap_plots(
  classification_selected_plot, surface_plot,
  ncol = 2, widths = c(0.92, 2.08)
) +
  patchwork::plot_annotation(
    title = "M_both species classification and T-P response surfaces"
  )
ggplot2::ggsave(
  paste0(surface_stem, ".pdf"), surface_plot,
  width = 23, height = 8.5, device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggplot2::ggsave(
  paste0(surface_stem, ".png"), surface_plot,
  width = 23, height = 8.5, dpi = 170, bg = "white", limitsize = FALSE
)
ggplot2::ggsave(
  paste0(combined_stem, ".pdf"), combined_plot,
  width = 31, height = 10.5, device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggplot2::ggsave(
  paste0(combined_stem, ".png"), combined_plot,
  width = 31, height = 10.5, dpi = 160, bg = "white", limitsize = FALSE
)

class_counts <- both_classification |>
  dplyr::count(.data$fitted_class, name = "n_species", .drop = FALSE) |>
  dplyr::mutate(percent = 100 * .data$n_species / sum(.data$n_species))
changed_count <- sum(classification_comparison$`changed? yes/no` == "yes")
cat("M_both figures completed without refitting.\n")
cat("Classification counts:\n")
print(class_counts)
cat("Changed versus M_context: ", changed_count, "/118 species\n", sep = "")
cat("Classification figure: ",
    normalizePath(paste0(classification_stem, ".pdf")), "\n", sep = "")
cat("Surface figure: ", normalizePath(paste0(surface_stem, ".pdf")),
    "\n", sep = "")
cat("Combined figure: ", normalizePath(paste0(combined_stem, ".pdf")),
    "\n", sep = "")
