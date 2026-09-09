# Step 7: classify full-rank species by fitted min/max Delta AUC --------------
#
# Loads the already fitted M_context model; no model is fitted here.
# For species with >=4 experiments and rank 4 for [1, T, P, T*P], predictions
# are evaluated at median warming only within each species' observed T-P hull.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "ggrepel", "Matrix", "metafor",
  "patchwork", "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

model_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "auc_change_species_context_model.rds"
)
output_dir <- file.path(
  "output", "step_by_step_modeling", "species_context"
)
figure_dir <- file.path(
  "output", "frequentist_hierarchical_auc", "figures"
)
figure_stem <- file.path(
  figure_dir, "step_07_context_minmax_classification"
)
summary_file <- file.path(
  output_dir, "context_full_rank_minmax_classification.csv"
)
count_file <- file.path(
  output_dir, "context_full_rank_minmax_classification_counts.csv"
)

if (!file.exists(model_file)) {
  stop("M_context was not found. Run R/16_step_04_species_context_model.R first.",
       call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

M_context <- readRDS(model_file)
fixed_beta <- stats::coef(M_context)
X_model <- as.matrix(M_context$X)

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(species_id = as.character(.data$species_id))

# Conditional random coefficients from D Z' M^-1 (y - X beta).
model_residual <- as.numeric(M_context$yi - X_model %*% M_context$beta)
weighted_residual <- Matrix::solve(M_context$M, model_residual)
species_levels <- levels(M_context$mf.g$outer)
n_species <- length(species_levels)
G_context <- as.matrix(M_context$G)
if (!identical(colnames(G_context), c("T_z", "P_z", "TP_z"))) {
  stop("Unexpected M_context random-slope columns.", call. = FALSE)
}

Z_intercept <- M_context$Z.S[[1]]
Z_context <- do.call(cbind, lapply(seq_len(ncol(G_context)), function(j) {
  Matrix::Diagonal(x = as.numeric(M_context$Z.G1[, j])) %*%
    M_context$Z.G2
}))
D_context <- Matrix::kronecker(G_context, Matrix::Diagonal(n_species))
random_intercept <- as.numeric(
  M_context$sigma2[[1]] *
    Matrix::crossprod(Z_intercept, weighted_residual)
)
random_context_vector <- as.numeric(
  D_context %*% Matrix::crossprod(Z_context, weighted_residual)
)
random_context_matrix <- matrix(
  random_context_vector,
  nrow = n_species,
  ncol = ncol(G_context),
  dimnames = list(species_levels, colnames(G_context))
)
random_coefficients <- tibble::tibble(
  species_id = species_levels,
  random_intercept = random_intercept,
  random_T = random_context_matrix[, "T_z"],
  random_P = random_context_matrix[, "P_z"],
  random_TP = random_context_matrix[, "TP_z"]
)

# Fixed-effect covariance and conditional prediction-error covariance for each
# species' four random coefficients. These intervals include fixed-effect and
# BLUP prediction uncertainty, conditional on the fitted variance parameters.
M_inverse_X <- Matrix::solve(M_context$M, X_model)
fixed_covariance <- solve(crossprod(X_model, M_inverse_X))
D_species <- matrix(
  0, nrow = 4L, ncol = 4L,
  dimnames = list(
    c("intercept", "T_z", "P_z", "TP_z"),
    c("intercept", "T_z", "P_z", "TP_z")
  )
)
D_species[1, 1] <- as.numeric(M_context$sigma2[[1]])
D_species[2:4, 2:4] <- G_context

species_support <- auc_data |>
  dplyr::group_split(.data$species_id, .keep = TRUE) |>
  lapply(function(species) {
    contexts <- dplyr::distinct(species, .data$T_z, .data$P_z)
    design <- cbind(
      intercept = 1,
      T_z = contexts$T_z,
      P_z = contexts$P_z,
      TP_z = contexts$T_z * contexts$P_z
    )
    tibble::tibble(
      species_id = species$species_id[[1]],
      n_experiments = dplyr::n_distinct(species$experiment_id),
      n_distinct_TP = nrow(contexts),
      context_design_rank = qr(design)$rank
    )
  }) |>
  dplyr::bind_rows() |>
  dplyr::filter(
    .data$n_experiments >= 4L,
    .data$context_design_rank == 4L
  ) |>
  dplyr::left_join(random_coefficients, by = "species_id")

if (!nrow(species_support)) {
  stop("No species have both >=4 experiments and full-rank T-P support.",
       call. = FALSE)
}

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
W_z_display <- stats::median(experiment_warming$W_z)
W_display <- stats::median(experiment_warming$W)

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
    dplyr::mutate(W_z = W_z_display, TP_z = .data$T_z * .data$P_z)
}

prediction_one_species <- function(species_name) {
  species_index <- match(species_name, species_levels)
  if (is.na(species_index)) {
    stop("Species is absent from the M_context random-effect levels: ",
         species_name, call. = FALSE)
  }
  grid <- species_hull_grid(species_name)
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, grid, names(fixed_beta)
  )
  random_design <- cbind(
    intercept = 1,
    T_z = grid$T_z,
    P_z = grid$P_z,
    TP_z = grid$TP_z
  )
  random_row <- random_coefficients |>
    dplyr::filter(.data$species_id == species_name)
  random_vector <- unname(as.numeric(random_row[
    1, c("random_intercept", "random_T", "random_P", "random_TP")
  ]))

  # Columns for this species in the intercept and coefficient-major slope
  # design matrices.
  Z_species <- cbind(
    Z_intercept[, species_index, drop = FALSE],
    Z_context[, species_index, drop = FALSE],
    Z_context[, n_species + species_index, drop = FALSE],
    Z_context[, 2L * n_species + species_index, drop = FALSE]
  )
  M_inverse_Z <- Matrix::solve(M_context$M, Z_species)
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

species_predictions <- dplyr::bind_rows(lapply(
  species_support$species_id, prediction_one_species
))

classification <- species_predictions |>
  dplyr::group_split(.data$species_id, .keep = TRUE) |>
  lapply(function(species) {
    minimum_row <- species[which.min(species$estimate), , drop = FALSE]
    maximum_row <- species[which.max(species$estimate), , drop = FALSE]
    tibble::tibble(
      species_id = species$species_id[[1]],
      min_fitted_delta_auc = minimum_row$estimate[[1]],
      min_prediction_se = minimum_row$prediction_se[[1]],
      min_conf_low = minimum_row$conf_low[[1]],
      min_conf_high = minimum_row$conf_high[[1]],
      min_T_z = minimum_row$T_z[[1]],
      min_P_z = minimum_row$P_z[[1]],
      max_fitted_delta_auc = maximum_row$estimate[[1]],
      max_prediction_se = maximum_row$prediction_se[[1]],
      max_conf_low = maximum_row$conf_low[[1]],
      max_conf_high = maximum_row$conf_high[[1]],
      max_T_z = maximum_row$T_z[[1]],
      max_P_z = maximum_row$P_z[[1]],
      n_support_grid_points = nrow(species)
    )
  }) |>
  dplyr::bind_rows() |>
  dplyr::left_join(
    dplyr::select(
      species_support, "species_id", "n_experiments", "n_distinct_TP",
      "context_design_rank"
    ),
    by = "species_id"
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
    warming_C = W_display,
    warming_z = W_z_display
  )

class_counts <- classification |>
  dplyr::count(.data$fitted_class, name = "n_species", .drop = FALSE) |>
  dplyr::mutate(percent = 100 * .data$n_species / sum(.data$n_species))

# Enforce the biological classification rules before plotting.
classification_valid <- with(
  classification,
  (fitted_class == "Consistent winner" & min_fitted_delta_auc > 0) |
    (fitted_class == "Consistent loser" & max_fitted_delta_auc < 0) |
    (fitted_class == "Context-dependent" &
       min_fitted_delta_auc <= 0 & max_fitted_delta_auc >= 0)
)
if (!all(classification_valid) ||
    any(classification$min_fitted_delta_auc >
          classification$max_fitted_delta_auc)) {
  stop("The fitted min/max classification rules failed validation.",
       call. = FALSE)
}

# Label only biologically prominent cases: three strongest winners, three
# strongest losers, and four switchers with the widest fitted response range.
label_species <- dplyr::bind_rows(
  classification |>
    dplyr::filter(.data$fitted_class == "Consistent winner") |>
    dplyr::slice_max(.data$min_fitted_delta_auc, n = 3L, with_ties = FALSE),
  classification |>
    dplyr::filter(.data$fitted_class == "Consistent loser") |>
    dplyr::slice_min(.data$max_fitted_delta_auc, n = 3L, with_ties = FALSE),
  classification |>
    dplyr::filter(.data$fitted_class == "Context-dependent") |>
    dplyr::mutate(response_span = .data$max_fitted_delta_auc -
                    .data$min_fitted_delta_auc) |>
    dplyr::slice_max(.data$response_span, n = 4L, with_ties = FALSE)
) |>
  dplyr::distinct(.data$species_id, .keep_all = TRUE)

# Set independent axes from the complete relevant CI ranges: horizontal
# minimum-response intervals determine x and vertical maximum-response
# intervals determine y. Thus no uncertainty whisker is clipped.
x_data_range <- range(c(
  classification$min_conf_low, classification$min_conf_high,
  classification$min_fitted_delta_auc, 0
), na.rm = TRUE)
y_data_range <- range(c(
  classification$max_conf_low, classification$max_conf_high,
  classification$max_fitted_delta_auc, 0
), na.rm = TRUE)
x_padding <- 0.08 * diff(x_data_range)
y_padding <- 0.08 * diff(y_data_range)
x_limits <- x_data_range + c(-x_padding, x_padding)
y_limits <- y_data_range + c(-y_padding, y_padding)

# Exact category colors used by the upstream figure pipeline in
# R/pipeline/figures.R. Purple is reused for the new context-dependent class.
class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#5E3C99",
  "Consistent winner" = "#1B9E77"
)

classification_plot <- ggplot2::ggplot(classification) +
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
    colour = "grey55", linewidth = 0.45
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.55) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.55) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_conf_low, xend = .data$min_conf_high,
      y = .data$max_fitted_delta_auc, yend = .data$max_fitted_delta_auc,
      colour = .data$fitted_class
    ),
    linewidth = 0.45, alpha = 0.25
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_fitted_delta_auc, xend = .data$min_fitted_delta_auc,
      y = .data$max_conf_low, yend = .data$max_conf_high,
      colour = .data$fitted_class
    ),
    linewidth = 0.45, alpha = 0.25
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      colour = .data$fitted_class,
      fill = .data$n_experiments
    ),
    shape = 21, size = 3.1, stroke = 0.9, alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = label_species,
    ggplot2::aes(
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      label = .data$species_id,
      colour = .data$fitted_class
    ),
    size = 3.1, show.legend = FALSE, min.segment.length = 0,
    segment.colour = "grey60", segment.size = 0.3,
    box.padding = 0.35, point.padding = 0.25,
    seed = 20260907
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[1]] + 0.17 * diff(x_limits),
    y = y_limits[[1]] + 0.13 * diff(y_limits),
    label = "CONSISTENT\nLOSERS", colour = class_colours[["Consistent loser"]],
    fontface = "bold", size = 3.6, alpha = 0.75
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[1]] + 0.17 * diff(x_limits),
    y = y_limits[[2]] - 0.10 * diff(y_limits),
    label = "CONTEXT-\nDEPENDENT",
    colour = class_colours[["Context-dependent"]],
    fontface = "bold", size = 3.6, alpha = 0.75
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[2]] - 0.17 * diff(x_limits),
    y = y_limits[[2]] - 0.10 * diff(y_limits),
    label = "CONSISTENT\nWINNERS",
    colour = class_colours[["Consistent winner"]],
    fontface = "bold", size = 3.6, alpha = 0.75
  ) +
  ggplot2::scale_colour_manual(values = class_colours, drop = FALSE) +
  ggplot2::scale_fill_viridis_c(
    option = "C", begin = 0.12, end = 0.88,
    name = "Number of\nexperiments"
  ) +
  ggplot2::coord_cartesian(
    xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    plot.margin = ggplot2::margin(8, 18, 8, 8)
  ) +
  ggplot2::labs(
    title = "Species classification by fitted minimum and maximum Delta AUC",
    subtitle = paste0(
      nrow(classification), " species with >=4 experiments and full-rank T-P ",
      "support; median warming (Delta T = ", sprintf("%.2f", W_display),
      " °C)."
    ),
    x = "Minimum fitted Delta AUC across observed T-P support",
    y = "Maximum fitted Delta AUC across observed T-P support",
    colour = "Fitted classification",
    caption = paste0(
      "Faint horizontal and vertical whiskers are pointwise model-based 95% ",
      "intervals at the fitted minimum and maximum locations.\nThey are not ",
      "simultaneous confidence intervals for the whole surface."
    )
  )

invisible(write_csv_atomic(classification, summary_file))
invisible(write_csv_atomic(class_counts, count_file))
ggplot2::ggsave(
  paste0(figure_stem, ".pdf"), classification_plot,
  width = 11, height = 9, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(figure_stem, ".png"), classification_plot,
  width = 11, height = 9, dpi = 220, bg = "white"
)

# Companion surfaces requested from the labelled upper-left cluster. The three
# otherwise-unlabelled screenshot points are Poa alpina, Agrostis capillaris,
# and Achillea millefolium. Potentilla "erracta" is normalized to the accepted
# name Potentilla erecta used by the analysis data.
surface_species <- c(
  "Silene acaulis", "Nardus stricta", "Festuca rubra", "Carex capillaris",
  "Potentilla erecta", "Poa alpina", "Agrostis capillaris",
  "Achillea millefolium"
)
missing_surface_species <- setdiff(surface_species, classification$species_id)
if (length(missing_surface_species)) {
  stop(
    "Requested surface species are not eligible: ",
    paste(missing_surface_species, collapse = ", "), call. = FALSE
  )
}

surface_metadata <- classification |>
  dplyr::filter(.data$species_id %in% surface_species) |>
  dplyr::mutate(
    species_order = match(.data$species_id, surface_species),
    panel = paste0(.data$species_id, " (", .data$n_experiments,
                   " experiments)")
  ) |>
  dplyr::arrange(.data$species_order)
panel_levels <- surface_metadata$panel

selected_surfaces <- species_predictions |>
  dplyr::filter(
    .data$species_id %in% surface_species,
    .data$support_point_type == "regular hull grid"
  ) |>
  dplyr::left_join(
    dplyr::select(surface_metadata, "species_id", "panel", "n_experiments"),
    by = "species_id"
  ) |>
  dplyr::mutate(panel = factor(.data$panel, levels = panel_levels)) |>
  dplyr::mutate(
    tile_width = diff(range(.data$T_z)) / 80 * 1.02,
    tile_height = diff(range(.data$P_z)) / 80 * 1.02,
    .by = "species_id"
  )

surface_observations <- auc_data |>
  dplyr::filter(.data$species_id %in% surface_species) |>
  dplyr::distinct(.data$species_id, .data$experiment_id, .data$T_z, .data$P_z) |>
  dplyr::left_join(
    dplyr::select(surface_metadata, "species_id", "panel"),
    by = "species_id"
  ) |>
  dplyr::mutate(panel = factor(.data$panel, levels = panel_levels))

scaling_value <- function(field, statistic) {
  scaled_auc$scaling[[statistic]][
    scaled_auc$scaling$scaled_field == field
  ][[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
surface_limit <- max(abs(selected_surfaces$estimate), na.rm = TRUE)

selected_surface_plot <- ggplot2::ggplot(
  selected_surfaces,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_tile(
    ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
  ) +
  ggplot2::geom_point(
    data = surface_observations,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey20",
    size = 1.25, stroke = 0.35
  ) +
  ggplot2::facet_wrap(~ panel, ncol = 3, scales = "free") +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-surface_limit, surface_limit)
  ) +
  ggplot2::scale_x_continuous(
    labels = function(x) sprintf("%.1f", x * T_scale + T_center)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
  ) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7.5),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Selected species T-P response surfaces",
    subtitle = paste0(
      "Median warming (Delta T = ", sprintf("%.2f", W_display),
      " °C); restricted to each species' observed T-P convex hull."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Fitted Delta AUC"
  )

surface_figure_stem <- file.path(
  figure_dir, "step_07_context_selected_species_tp_surfaces"
)
ggplot2::ggsave(
  paste0(surface_figure_stem, ".pdf"), selected_surface_plot,
  width = 13, height = 10, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(surface_figure_stem, ".png"), selected_surface_plot,
  width = 13, height = 10, dpi = 200, bg = "white"
)

classification_with_surfaces <- patchwork::wrap_plots(
  classification_plot, selected_surface_plot,
  ncol = 2, widths = c(1.05, 1.25)
) +
  patchwork::plot_annotation(
    title = "Species min-max classification and selected climate surfaces"
  )
combined_figure_stem <- file.path(
  figure_dir, "step_07_context_classification_with_species_surfaces"
)
ggplot2::ggsave(
  paste0(combined_figure_stem, ".pdf"), classification_with_surfaces,
  width = 25, height = 11, device = grDevices::cairo_pdf, limitsize = FALSE
)
ggplot2::ggsave(
  paste0(combined_figure_stem, ".png"), classification_with_surfaces,
  width = 25, height = 11, dpi = 170, bg = "white", limitsize = FALSE
)

cat("M_context was loaded without refitting.\n")
cat("Median warming: ", sprintf("%.2f", W_display), " °C\n", sep = "")
cat("Eligible species: ", nrow(classification), "\n", sep = "")
print(class_counts)
cat("Figure: ", normalizePath(paste0(figure_stem, ".pdf")), "\n", sep = "")
cat("Selected surfaces: ",
    normalizePath(paste0(surface_figure_stem, ".pdf")), "\n", sep = "")
cat("Combined figure: ",
    normalizePath(paste0(combined_figure_stem, ".pdf")), "\n", sep = "")
cat("Table: ", normalizePath(summary_file), "\n", sep = "")
