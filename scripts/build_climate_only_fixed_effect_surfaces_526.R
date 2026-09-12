#!/usr/bin/env Rscript

# Figure 5 population-average prediction surfaces from the canonical
# climate-only M_both_526 model. This script only reads existing exported
# results and does not fit models or invoke the targets pipeline.

required_packages <- c(
  "cowplot", "dplyr", "geometry", "ggplot2", "ggrepel", "metafor",
  "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
if (!file.exists("_targets.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}

suppressPackageStartupMessages(library(metafor))
source("R/functions/auc_model_helpers.R")

arguments <- commandArgs(trailingOnly = TRUE)
output_stem <- if (length(arguments) >= 1L) {
  arguments[[1]]
} else {
  file.path(
    "figures", "candidates",
    "climate_only_fixed_effect_prediction_surfaces_526"
  )
}
output_root <- if (length(arguments) >= 2L) arguments[[2]] else "output"

pdf_file <- paste0(output_stem, ".pdf")
png_file <- paste0(output_stem, ".png")
output_name <- basename(output_stem)
data_file <- file.path(output_root, "data", paste0(output_name, "_data.csv"))
quantile_file <- file.path(
  output_root, "tables", paste0(output_name, "_quantiles.csv")
)
validation_file <- file.path(
  output_root, "audits", paste0(output_name, "_validation.csv")
)

model_file <- "output/models/M_both_526.rds"
coefficient_file <- "output/tables/M_both_526_fixed_coefficients.csv"
scaling_file <- "output/data/downstream_climate_scaling_526.rds"
experiment_file <- "output/tables/experiments.csv"
site_file <- "output/tables/sites.csv"
input_files <- c(
  model_file, coefficient_file, scaling_file, experiment_file, site_file
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    "Required canonical input(s) missing: ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

model <- readRDS(model_file)
coefficient_export <- readr::read_csv(
  coefficient_file, show_col_types = FALSE
)
climate_scaling <- readRDS(scaling_file)
experiments <- readr::read_csv(experiment_file, show_col_types = FALSE)
sites <- readr::read_csv(site_file, show_col_types = FALSE)
if (!inherits(model, "rma.mv") || !is.data.frame(model$data)) {
  stop("M_both_526.rds is not the expected fitted rma.mv export.", call. = FALSE)
}

model_formula <- stats::as.formula("~ T_z * P_z * W_z")
expected_terms <- c(
  "intrcpt", "T_z", "P_z", "W_z", "T_z:P_z", "T_z:W_z",
  "P_z:W_z", "T_z:P_z:W_z"
)
if (!identical(names(stats::coef(model)), expected_terms)) {
  stop("The canonical climate-only coefficient terms changed.", call. = FALSE)
}
coefficient_export <- coefficient_export |>
  dplyr::arrange(match(.data$term, .env$expected_terms))
coefficient_error <- max(abs(
  as.numeric(stats::coef(model)) - coefficient_export$estimate
))
if (!is.finite(coefficient_error) || coefficient_error > 1e-12) {
  stop("The model and exported fixed-coefficient table disagree.", call. = FALSE)
}

single_scaling_value <- function(field, statistic) {
  value <- climate_scaling[[statistic]][
    climate_scaling$scaled_field == field
  ]
  if (length(value) != 1L || !is.finite(value)) {
    stop("Missing unique scaling value for ", field, ".", call. = FALSE)
  }
  as.numeric(value[[1]])
}

T_center <- single_scaling_value("T_z", "center")
T_scale <- single_scaling_value("T_z", "scale")
P_center <- single_scaling_value("P_z", "center")
P_scale <- single_scaling_value("P_z", "scale")
W_center <- single_scaling_value("W_z", "center")
W_scale <- single_scaling_value("W_z", "scale")

reference_experiments <- model$data |>
  dplyr::distinct(
    .data$experiment_id, .data$T, .data$P, .data$P_raw, .data$W,
    .data$T_z, .data$P_z, .data$W_z
  ) |>
  dplyr::mutate(experiment_id = as.character(.data$experiment_id)) |>
  dplyr::arrange(.data$experiment_id)
if (nrow(reference_experiments) != 42L) {
  stop("The canonical climate model no longer contains 42 experiments.",
       call. = FALSE)
}

climate_scaling_error <- max(c(
  abs(reference_experiments$T_z -
        (reference_experiments$T - T_center) / T_scale),
  abs(reference_experiments$P_z -
        (reference_experiments$P - P_center) / P_scale),
  abs(reference_experiments$W_z -
        (reference_experiments$W - W_center) / W_scale)
))
if (!is.finite(climate_scaling_error) || climate_scaling_error > 1e-10) {
  stop("Stored climate columns do not match the canonical scaling export.",
       call. = FALSE)
}

origin_sites <- sites |>
  dplyr::select("region", "site_id", "elevation") |>
  dplyr::rename(
    origin_site = "site_id",
    origin_elevation = "elevation"
  )
destination_sites <- sites |>
  dplyr::select("region", "site_id", "elevation") |>
  dplyr::rename(
    destination_site = "site_id",
    destination_elevation = "elevation"
  )
experiment_elevations <- experiments |>
  dplyr::mutate(experiment_id = as.character(.data$experiment_id)) |>
  dplyr::semi_join(
    reference_experiments |> dplyr::select("experiment_id"),
    by = "experiment_id"
  ) |>
  dplyr::select(
    "experiment_id", "region", "origin_site", "destination_site"
  ) |>
  dplyr::left_join(
    origin_sites, by = c("region", "origin_site")
  ) |>
  dplyr::left_join(
    destination_sites, by = c("region", "destination_site")
  ) |>
  dplyr::mutate(
    experiment_label = paste0(
      .data$region, "\n",
      sprintf("%.0f", .data$origin_elevation), " → ",
      sprintf("%.0f", .data$destination_elevation)
    )
  )
if (nrow(experiment_elevations) != 42L ||
    anyDuplicated(experiment_elevations$experiment_id) ||
    any(!is.finite(experiment_elevations$origin_elevation)) ||
    any(!is.finite(experiment_elevations$destination_elevation))) {
  stop(
    "Every canonical experiment must resolve to one finite elevation pair.",
    call. = FALSE
  )
}

warming_raw <- stats::quantile(
  reference_experiments$W,
  probs = c(0.25, 0.75),
  names = FALSE,
  type = 7
)
warming_table <- tibble::tibble(
  warming_level = factor(
    c("Low warming", "High warming"),
    levels = c("Low warming", "High warming")
  ),
  quantile = c("Q25", "Q75"),
  W_raw = warming_raw,
  W_z = (warming_raw - W_center) / W_scale,
  raw_unit = "degrees C",
  n_experiments = nrow(reference_experiments),
  weighting = "42 distinct experiments; each experiment contributes once"
)

T_grid <- seq(
  min(reference_experiments$T_z), max(reference_experiments$T_z),
  length.out = 81L
)
P_grid <- seq(
  min(reference_experiments$P_z), max(reference_experiments$P_z),
  length.out = 81L
)
support_experiments <- reference_experiments |>
  dplyr::select("experiment_id", "T_z", "P_z", "W_z")

make_surface <- function(warming_index) {
  warming_level <- warming_table$warming_level[[warming_index]]
  W_raw <- warming_table$W_raw[[warming_index]]
  W_z <- warming_table$W_z[[warming_index]]
  candidate_grid <- tidyr::expand_grid(
    T_z = T_grid,
    P_z = P_grid,
    W_z = W_z
  )
  support_result <- auc_support_mask(support_experiments, candidate_grid)
  if (nrow(support_result$grid) != nrow(candidate_grid) ||
      !identical(
        support_result$grid$support_method,
        rep("three-dimensional convex hull", nrow(candidate_grid))
      )) {
    stop("The climate support mask could not be constructed.", call. = FALSE)
  }

  grid <- support_result$grid |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  design <- auc_design_matrix(
    model_formula, grid, names(stats::coef(model))
  )
  fixed_prediction <- as.numeric(design %*% stats::coef(model))
  grid |>
    dplyr::mutate(
      row_type = "fixed_effect_prediction_grid",
      model = "M_both_526",
      model_display = factor(
        "Climate-only\nmodel", levels = "Climate-only\nmodel"
      ),
      warming_level = factor(
        warming_level, levels = levels(warming_table$warming_level)
      ),
      W_raw = W_raw,
      destination_temperature = .data$T_z * T_scale + T_center,
      destination_precipitation = exp(.data$P_z * P_scale + P_center),
      predicted_winning_losing_strength = dplyr::if_else(
        .data$supported, .env$fixed_prediction, NA_real_
      ),
      prediction_component = "fixed effects only",
      support_scope = "observed three-dimensional climate-warming hull"
    )
}

predictions <- dplyr::bind_rows(
  lapply(seq_len(nrow(warming_table)), make_surface)
)
make_observed_points <- function(warming_index) {
  warming_level <- warming_table$warming_level[[warming_index]]
  W_raw <- warming_table$W_raw[[warming_index]]
  W_z <- warming_table$W_z[[warming_index]]
  panel_points <- reference_experiments |>
    dplyr::left_join(experiment_elevations, by = "experiment_id") |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      experiment_label = .data$experiment_label,
      origin_elevation = .data$origin_elevation,
      destination_elevation = .data$destination_elevation,
      destination_temperature = .data$T,
      destination_precipitation = .data$P_raw,
      T_z = .data$T_z,
      P_z = .data$P_z,
      W_z = .env$W_z,
      observed_experiment_W_raw = .data$W
    )
  point_support <- auc_support_mask(
    support_experiments,
    panel_points |> dplyr::select("T_z", "P_z", "W_z")
  )$grid$supported
  panel_points |>
    dplyr::mutate(
      row_type = "observed_experiment",
      model = "M_both_526",
      model_display = factor(
        "Climate-only\nmodel", levels = "Climate-only\nmodel"
      ),
      warming_level = factor(
        .env$warming_level, levels = levels(warming_table$warming_level)
      ),
      W_raw = .env$W_raw,
      TP_z = .data$T_z * .data$P_z,
      predicted_winning_losing_strength = NA_real_,
      supported = TRUE,
      supported_at_panel_warming = point_support,
      support_method = "panel-specific three-dimensional convex-hull test",
      prediction_component = "fixed effects only",
      support_scope = "observed three-dimensional climate-warming hull"
    )
}

observed <- dplyr::bind_rows(
  lapply(seq_len(nrow(warming_table)), make_observed_points)
)

surface_plot_data <- predictions |>
  dplyr::filter(
    .data$supported,
    is.finite(.data$predicted_winning_losing_strength)
  )
global_limit <- max(abs(surface_plot_data$predicted_winning_losing_strength))
if (!is.finite(global_limit) || global_limit <= 0) {
  stop("The supported fixed-effect prediction range is invalid.", call. = FALSE)
}

T_breaks <- pretty(range(T_grid), n = 3)
T_breaks <- T_breaks[T_breaks >= min(T_grid) & T_breaks <= max(T_grid)]
P_breaks <- pretty(range(P_grid), n = 3)
P_breaks <- P_breaks[P_breaks >= min(P_grid) & P_breaks <= max(P_grid)]

surface_core <- ggplot2::ggplot(
  predictions,
  ggplot2::aes(x = .data$T_z, y = .data$P_z)
) +
  ggplot2::geom_raster(
    data = surface_plot_data,
    ggplot2::aes(fill = .data$predicted_winning_losing_strength)
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(observed, !.data$supported_at_panel_warming),
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "grey78",
    alpha = 0.55,
    size = 0.95,
    stroke = 0.30
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(observed, .data$supported_at_panel_warming),
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "grey20",
    alpha = 0.85,
    size = 1.20,
    stroke = 0.38
  ) +
  ggrepel::geom_text_repel(
    data = observed,
    ggplot2::aes(
      x = .data$T_z,
      y = .data$P_z,
      label = .data$experiment_label,
      colour = .data$supported_at_panel_warming
    ),
    inherit.aes = FALSE,
    size = 4.7,
    lineheight = 0.88,
    seed = 20260824,
    box.padding = 0.35,
    point.padding = 0.18,
    min.segment.length = 0,
    segment.colour = "grey72",
    segment.size = 0.28,
    force = 5,
    force_pull = 0.12,
    max.time = 12,
    max.iter = 200000,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::facet_grid(cols = ggplot2::vars(.data$warming_level)) +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "grey57", `TRUE` = "grey25"),
    guide = "none"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#E66101",
    mid = "white",
    high = "#1B9E77",
    midpoint = 0,
    limits = c(-global_limit, global_limit),
    na.value = "white",
    name = "Predicted winning-losing strength"
  ) +
  ggplot2::scale_x_continuous(
    breaks = T_breaks,
    labels = function(x) sprintf("%.1f", x * T_scale + T_center),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = P_breaks,
    labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center)),
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::labs(
    x = "Destination temperature (°C)",
    y = "Destination precipitation (mm)"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(14.5, "cm"),
      barheight = grid::unit(0.70, "cm")
    )
  ) +
  ggplot2::theme_bw(base_size = 23) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.spacing = grid::unit(1.6, "mm"),
    panel.border = ggplot2::element_rect(colour = "grey35", linewidth = 0.65),
    strip.placement = "outside",
    strip.text.x = ggplot2::element_blank(),
    strip.background.x = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 27),
    axis.text = ggplot2::element_text(size = 19, colour = "black"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 25, face = "bold"),
    legend.text = ggplot2::element_text(size = 21),
    plot.margin = ggplot2::margin(t = 27, r = 8, b = 8, l = 8)
  )

surface_figure <- cowplot::ggdraw() +
  cowplot::draw_plot(
    surface_core, x = 0, y = 0, width = 1, height = 0.982
  ) +
  cowplot::draw_label(
    "Low warming", x = 0.295, y = 0.983,
    fontface = "bold", size = 27
  ) +
  cowplot::draw_label(
    "High warming", x = 0.725, y = 0.983,
    fontface = "bold", size = 27
  )

validation <- tibble::tribble(
  ~check, ~pass, ~detail,
  "No model refit", TRUE,
  "Read output/models/M_both_526.rds; the script contains no fitting call and does not invoke targets",
  "Canonical coefficient identity", coefficient_error <= 1e-12,
  paste("Maximum difference from exported fixed coefficients:",
        format(coefficient_error, scientific = TRUE)),
  "Canonical analysis scope",
  model$k == 1604L &&
    dplyr::n_distinct(model$data$species_id) == 526L &&
    nrow(reference_experiments) == 42L,
  paste(model$k, "effects;", dplyr::n_distinct(model$data$species_id),
        "species;", nrow(reference_experiments), "experiments"),
  "Fixed-effect predictions", TRUE,
  "Every surface value is model.matrix(~ T_z * P_z * W_z) multiplied by the eight stored fixed coefficients; no random effect is added",
  "Climate scaling identity", climate_scaling_error <= 1e-10,
  paste("Maximum standardized-climate reconstruction error:",
        format(climate_scaling_error, scientific = TRUE)),
  "Canonical elevation labels",
  nrow(experiment_elevations) == 42L &&
    !anyDuplicated(experiment_elevations$experiment_id) &&
    all(is.finite(experiment_elevations$origin_elevation)) &&
    all(is.finite(experiment_elevations$destination_elevation)),
  "All 42 experiment labels use unique origin and destination elevations from output/tables/sites.csv; elevation units are omitted from the labels",
  "Support-mask enforcement",
  all(is.finite(
    predictions$predicted_winning_losing_strength[predictions$supported]
  )) && all(is.na(
    predictions$predicted_winning_losing_strength[!predictions$supported]
  )),
  "Only cells inside the observed three-dimensional temperature-precipitation-warming hull receive predictions",
  "Panel-specific experiment styling",
  all(c(TRUE, FALSE) %in% observed$supported_at_panel_warming),
  paste(
    sum(observed$supported_at_panel_warming),
    "projected experiment points are supported and",
    sum(!observed$supported_at_panel_warming),
    "are outside their displayed fixed-warming slice"
  ),
  "All experiments annotated",
  all(
    dplyr::count(observed, .data$warming_level)$n == 42L
  ) && all(!is.na(observed$experiment_label)) &&
    all(nzchar(observed$experiment_label)),
  "Each warming panel contains a non-empty two-line site and elevation-change label for all 42 canonical experiments"
)
if (!all(validation$pass)) {
  stop("One or more candidate-figure validation checks failed.", call. = FALSE)
}

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(data_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(quantile_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)

ggplot2::ggsave(
  pdf_file, surface_figure,
  width = 24, height = 13, units = "in", device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  png_file, surface_figure,
  width = 24, height = 13, units = "in", dpi = 300, bg = "white"
)
readr::write_csv(
  dplyr::bind_rows(predictions, observed), data_file, na = ""
)
readr::write_csv(warming_table, quantile_file, na = "")
readr::write_csv(validation, validation_file, na = "")
invisible(lapply(
  c(pdf_file, png_file),
  function(path) {
    writeLines(
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      paste0(path, ".generated")
    )
  }
))

cat("Created:", pdf_file, "\n")
cat("Created:", png_file, "\n")
cat("Created:", data_file, "\n")
cat("Created:", quantile_file, "\n")
cat("Created:", validation_file, "\n")
cat(
  "Supported prediction range:",
  sprintf(
    "%.4f to %.4f",
    min(surface_plot_data$predicted_winning_losing_strength),
    max(surface_plot_data$predicted_winning_losing_strength)
  ),
  "\n"
)
