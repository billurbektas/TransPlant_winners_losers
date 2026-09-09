#!/usr/bin/env Rscript

# Step 23: updated 526-species classification with linked response surfaces
#
# Recreates the established neutral Step 12 linked-callout design from the
# already fitted M_both_526. No model is fitted or altered here. The four
# previously selected biological examples are retained and every surface is
# restricted to that species' own observed destination T-P convex hull.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "cowplot", "dplyr", "geometry", "ggplot2", "here", "metafor", "readr",
  "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

root <- here::here()
output_dir <- file.path(root, "output", "step_by_step_modeling", "downstream_526")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "models")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  model = file.path(model_dir, "M_both_526.rds"),
  data = file.path(table_dir, "auc_data_526.csv"),
  classification = file.path(table_dir, "species_minmax_classification_526.csv"),
  random = file.path(table_dir, "M_both_526_species_random_coefficients.csv"),
  scaling = file.path(table_dir, "climate_scaling_526.csv")
)
if (any(!file.exists(paths))) {
  stop("Run Step 22 first; missing: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
}

M_both_526 <- readRDS(paths[["model"]])
auc_data <- readr::read_csv(paths[["data"]], show_col_types = FALSE) |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    experiment_id = as.character(.data$experiment_id)
  )
classification <- readr::read_csv(
  paths[["classification"]], show_col_types = FALSE
) |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    fitted_class = factor(
      .data$fitted_class,
      levels = c("Consistent loser", "Context-dependent", "Consistent winner")
    ),
    response_span = .data$max_fitted_delta_auc - .data$min_fitted_delta_auc
  )
random_coefficients <- readr::read_csv(paths[["random"]], show_col_types = FALSE) |>
  dplyr::mutate(species_id = as.character(.data$species_id))
scaling <- readr::read_csv(paths[["scaling"]], show_col_types = FALSE)

if (nrow(classification) != 116L ||
    any(classification$n_experiments < 4L) ||
    any(classification$context_design_rank != 4L) ||
    anyNA(classification$random_W)) {
  stop("The updated eligible-species classification is incomplete.", call. = FALSE)
}

# These are the final examples explicitly retained for the Step 12 design.
# Under M_both_526, Potentilla is the highest-minimum winner, Plantago is one
# of the two lowest-maximum losers, and Elymus/Silene have the two widest
# context-dependent fitted spans.
requested_species <- c(
  "Potentilla erecta", "Plantago lanceolata", "Elymus repens",
  "Silene acaulis"
)
selection <- classification |>
  dplyr::filter(.data$species_id %in% .env$requested_species) |>
  dplyr::mutate(display_order = match(.data$species_id, .env$requested_species)) |>
  dplyr::arrange(.data$display_order)
if (!identical(selection$species_id, requested_species)) {
  stop("One or more established callout species are absent from the updated eligible set.",
       call. = FALSE)
}
expected_classes <- c(
  "Consistent winner", "Consistent loser", "Context-dependent",
  "Context-dependent"
)
if (!identical(as.character(selection$fitted_class), expected_classes)) {
  stop("An established callout species changed fitted class; do not silently relabel it.",
       call. = FALSE)
}

class_rank <- classification |>
  dplyr::group_by(.data$fitted_class) |>
  dplyr::arrange(
    dplyr::case_when(
      .data$fitted_class == "Consistent winner" ~ -.data$min_fitted_delta_auc,
      .data$fitted_class == "Consistent loser" ~ .data$max_fitted_delta_auc,
      TRUE ~ -.data$response_span
    ),
    .data$species_id,
    .by_group = TRUE
  ) |>
  dplyr::mutate(selection_rank_within_class = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select("species_id", "selection_rank_within_class")
selection <- selection |>
  dplyr::left_join(class_rank, by = "species_id") |>
  dplyr::mutate(
    selection_rule = dplyr::case_when(
      .data$fitted_class == "Consistent winner" ~ "highest fitted minimum",
      .data$fitted_class == "Consistent loser" ~ "lowest fitted maximum",
      TRUE ~ "widest fitted min-max span"
    )
  )
write_csv_atomic(
  selection,
  file.path(table_dir, "species_minmax_callout_examples_526.csv")
)

fixed_beta <- stats::coef(M_both_526)
warming_experiments <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
warming_lookup <- tibble::tibble(
  warming_name = c("Low warming", "High warming"),
  probability = c(0.25, 0.75),
  W_z = as.numeric(stats::quantile(warming_experiments$W_z, c(0.25, 0.75))),
  W_C = as.numeric(stats::quantile(warming_experiments$W, c(0.25, 0.75)))
)
warming_panel_levels <- paste0(
  warming_lookup$warming_name, "\nDelta T = ",
  sprintf("%.2f", warming_lookup$W_C), " °C"
)

surface_one <- function(species_name, W_z_value, grid_points = 81L) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == .env$species_name) |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
  if (nrow(observed) < 4L ||
      qr(cbind(1, observed$T_z, observed$P_z,
               observed$T_z * observed$P_z))$rank != 4L) {
    stop(species_name, " no longer has full-rank T-P support.", call. = FALSE)
  }
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = grid_points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = grid_points)
  )
  hull <- geometry::convhulln(as.matrix(observed[, c("T_z", "P_z")]), options = "Qt")
  candidate <- candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::mutate(W_z = W_z_value, TP_z = .data$T_z * .data$P_z)
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, candidate, names(fixed_beta)
  )
  random <- random_coefficients |>
    dplyr::filter(.data$species_id == .env$species_name)
  if (nrow(random) != 1L) stop("Random coefficients unavailable for ", species_name)
  candidate |>
    dplyr::mutate(
      estimate = as.numeric(fixed_design %*% fixed_beta) +
        random$random_intercept[[1]] + random$random_W[[1]] * .data$W_z +
        random$random_T[[1]] * .data$T_z + random$random_P[[1]] * .data$P_z +
        random$random_TP[[1]] * .data$TP_z,
      tile_width = diff(range(observed$T_z)) / (grid_points - 1L) * 1.025,
      tile_height = diff(range(observed$P_z)) / (grid_points - 1L) * 1.025
    ) |>
    list(grid = _, observed = observed)
}

surface_sets <- stats::setNames(lapply(requested_species, function(species_name) {
  dplyr::bind_rows(lapply(seq_len(nrow(warming_lookup)), function(index) {
    surface_one(species_name, warming_lookup$W_z[[index]])$grid |>
      dplyr::mutate(
        warming_panel = factor(
          warming_panel_levels[[index]], levels = warming_panel_levels
        )
      )
  }))
}), requested_species)
observation_sets <- stats::setNames(lapply(requested_species, function(species_name) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == .env$species_name) |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
  tidyr::crossing(
    observed,
    warming_panel = factor(warming_panel_levels, levels = warming_panel_levels)
  )
}), requested_species)

surface_limit <- max(abs(unlist(lapply(surface_sets, `[[`, "estimate"))), na.rm = TRUE)
warming_limit <- max(abs(classification$random_W), na.rm = TRUE)
class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#0072B2",
  "Consistent winner" = "#1B9E77"
)
neutral_colour <- "#3F3F3F"

x_range <- range(c(
  classification$min_CI_low, classification$min_CI_high,
  classification$min_fitted_delta_auc, 0
), na.rm = TRUE)
y_range <- range(c(
  classification$max_CI_low, classification$max_CI_high,
  classification$max_fitted_delta_auc, 0
), na.rm = TRUE)
x_limits <- x_range + c(-0.08, 0.08) * diff(x_range)
y_limits <- y_range + c(-0.08, 0.08) * diff(y_range)

class_counts <- classification |>
  dplyr::count(.data$fitted_class, name = "n_species")
class_n <- function(class_name) {
  class_counts$n_species[as.character(class_counts$fitted_class) == class_name]
}

# Central panel geometry matches the final neutral Step 12 version. Connector
# segments originate at the exact data coordinates of highlighted point
# centers, reach the panel border in data coordinates, then continue to the
# callout frames in the outer cowplot canvas.
class_box <- c(x = 0.265, y = 0.155, width = 0.47, height = 0.69)
inset_boxes <- list(
  "Elymus repens" = c(x = 0.012, y = 0.575, width = 0.225, height = 0.39),
  "Silene acaulis" = c(x = 0.012, y = 0.075, width = 0.225, height = 0.39),
  "Potentilla erecta" = c(x = 0.763, y = 0.535, width = 0.225, height = 0.39),
  "Plantago lanceolata" = c(x = 0.763, y = 0.075, width = 0.225, height = 0.39)
)
panel_box <- c(left = 0.294, right = 0.734, bottom = 0.231, top = 0.843)
callout_centres_y <- vapply(
  inset_boxes, function(box) box[["y"]] + box[["height"]] / 2, numeric(1)
)
callout_centres_y <- pmin(
  pmax(callout_centres_y, panel_box[["bottom"]]), panel_box[["top"]]
)
canvas_y_to_data <- function(canvas_y) {
  y_limits[[1]] +
    (canvas_y - panel_box[["bottom"]]) /
    (panel_box[["top"]] - panel_box[["bottom"]]) * diff(y_limits)
}
connector_data <- selection |>
  dplyr::transmute(
    species_id = .data$species_id,
    x = .data$min_fitted_delta_auc,
    y = .data$max_fitted_delta_auc,
    xend = ifelse(
      .data$species_id %in% c("Elymus repens", "Silene acaulis"),
      x_limits[[1]], x_limits[[2]]
    ),
    yend = vapply(
      .data$species_id,
      function(species_name) canvas_y_to_data(callout_centres_y[[species_name]]),
      numeric(1)
    )
  )

classification_plot <- ggplot2::ggplot(classification) +
  ggplot2::annotate(
    "rect", xmin = x_limits[[1]], xmax = 0, ymin = y_limits[[1]], ymax = 0,
    fill = class_colours[["Consistent loser"]], alpha = 0.07
  ) +
  ggplot2::annotate(
    "rect", xmin = x_limits[[1]], xmax = 0, ymin = 0, ymax = y_limits[[2]],
    fill = class_colours[["Context-dependent"]], alpha = 0.07
  ) +
  ggplot2::annotate(
    "rect", xmin = 0, xmax = x_limits[[2]], ymin = 0, ymax = y_limits[[2]],
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
      x = .data$min_CI_low, xend = .data$min_CI_high,
      y = .data$max_fitted_delta_auc, yend = .data$max_fitted_delta_auc,
      colour = .data$fitted_class
    ), linewidth = 0.35, alpha = 0.17
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_fitted_delta_auc, xend = .data$min_fitted_delta_auc,
      y = .data$max_CI_low, yend = .data$max_CI_high,
      colour = .data$fitted_class
    ), linewidth = 0.35, alpha = 0.17
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$min_fitted_delta_auc, y = .data$max_fitted_delta_auc,
      colour = .data$fitted_class, fill = .data$random_W
    ), shape = 21, size = 2.5, stroke = 0.75, alpha = 0.94
  ) +
  ggplot2::geom_segment(
    data = connector_data,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    inherit.aes = FALSE, colour = neutral_colour,
    linewidth = 0.85, lineend = "round"
  ) +
  ggplot2::geom_point(
    data = selection,
    ggplot2::aes(
      x = .data$min_fitted_delta_auc, y = .data$max_fitted_delta_auc,
      fill = .data$random_W
    ), inherit.aes = FALSE, shape = 21, size = 4.2, stroke = 1.25,
    colour = neutral_colour, show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[1]] + 0.16 * diff(x_limits),
    y = y_limits[[1]] + 0.09 * diff(y_limits),
    label = paste0("CONSISTENT LOSER\nn = ", class_n("Consistent loser")),
    colour = class_colours[["Consistent loser"]], fontface = "bold",
    size = 3.1, hjust = 0
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[1]] + 0.05 * diff(x_limits),
    y = y_limits[[2]] - 0.07 * diff(y_limits),
    label = paste0("CONTEXT-DEPENDENT\nn = ", class_n("Context-dependent")),
    colour = class_colours[["Context-dependent"]], fontface = "bold",
    size = 3.1, hjust = 0
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[2]] - 0.04 * diff(x_limits),
    y = y_limits[[2]] - 0.07 * diff(y_limits),
    label = paste0("CONSISTENT WINNER\nn = ", class_n("Consistent winner")),
    colour = class_colours[["Consistent winner"]], fontface = "bold",
    size = 3.1, hjust = 1
  ) +
  ggplot2::scale_colour_manual(values = class_colours, drop = FALSE, guide = "none") +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-warming_limit, warming_limit),
    name = "Species random W\nslope deviation"
  ) +
  ggplot2::coord_cartesian(
    xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom", legend.box = "vertical",
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  ) +
  ggplot2::labs(
    x = "Minimum fitted Delta AUC",
    y = "Maximum fitted Delta AUC",
    colour = NULL
  )

scaling_value <- function(field, statistic) {
  scaling[[statistic]][scaling$scaled_field == field][[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
raw_T_labels <- function(x) sprintf("%.1f", x * T_scale + T_center)
raw_P_labels <- function(x) sprintf("%.0f", exp(x * P_scale + P_center))

surface_plot <- function(species_name) {
  ggplot2::ggplot(
    surface_sets[[species_name]],
    ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
    ) +
    ggplot2::geom_point(
      data = observation_sets[[species_name]],
      ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
      shape = 21, fill = "white", colour = "grey15", size = 1.05, stroke = 0.3
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$warming_panel)) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
      limits = c(-surface_limit, surface_limit),
      guide = ggplot2::guide_colourbar(
        title.position = "top", direction = "horizontal",
        barwidth = grid::unit(35, "mm"), barheight = grid::unit(3, "mm")
      )
    ) +
    ggplot2::scale_x_continuous(labels = raw_T_labels) +
    ggplot2::scale_y_continuous(labels = raw_P_labels) +
    ggplot2::theme_bw(base_size = 7.1) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(face = "bold", size = 6.5),
      plot.title = ggplot2::element_text(
        face = "bold.italic", size = 9, colour = neutral_colour, hjust = 0.5
      ),
      plot.background = ggplot2::element_rect(
        colour = neutral_colour, linewidth = 1, fill = "white"
      ),
      plot.margin = ggplot2::margin(4, 4, 3, 4),
      axis.title = ggplot2::element_text(size = 6.7),
      axis.text = ggplot2::element_text(size = 6.1),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.7),
      legend.text = ggplot2::element_text(size = 6.1),
      legend.margin = ggplot2::margin(-2, 0, -2, 0)
    ) +
    ggplot2::labs(
      title = paste0(
        species_name, " (",
        dplyr::filter(selection, .data$species_id == .env$species_name)$n_experiments,
        " experiments)"
      ),
      x = "Destination temperature (°C)",
      y = "Destination precipitation",
      fill = "Fitted AUC change"
    )
}

mini_plots <- stats::setNames(
  lapply(requested_species, surface_plot), requested_species
)

connector_target <- function(species_name) {
  box <- inset_boxes[[species_name]]
  if (box[["x"]] < class_box[["x"]]) {
    c(x = box[["x"]] + box[["width"]], y = box[["y"]] + box[["height"]] / 2)
  } else {
    c(x = box[["x"]], y = box[["y"]] + box[["height"]] / 2)
  }
}

canvas <- cowplot::ggdraw() +
  cowplot::draw_plot(
    classification_plot,
    x = class_box[["x"]], y = class_box[["y"]],
    width = class_box[["width"]], height = class_box[["height"]]
  )
for (species_name in names(inset_boxes)) {
  is_left <- species_name %in% c("Elymus repens", "Silene acaulis")
  start <- c(
    x = if (is_left) panel_box[["left"]] else panel_box[["right"]],
    y = callout_centres_y[[species_name]]
  )
  finish <- connector_target(species_name)
  canvas <- canvas + cowplot::draw_line(
    x = c(start[["x"]], finish[["x"]]),
    y = c(start[["y"]], finish[["y"]]),
    colour = neutral_colour, linewidth = 0.85, alpha = 0.9
  )
}
for (species_name in names(inset_boxes)) {
  box <- inset_boxes[[species_name]]
  canvas <- canvas + cowplot::draw_plot(
    mini_plots[[species_name]],
    x = box[["x"]], y = box[["y"]],
    width = box[["width"]], height = box[["height"]]
  )
}

output_stem <- file.path(figure_dir, "species_minmax_classification_526")
ggplot2::ggsave(
  paste0(output_stem, ".pdf"), canvas,
  width = 16, height = 12.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(output_stem, ".png"), canvas,
  width = 16, height = 12.5, dpi = 190, bg = "white"
)

message("Updated linked-callout classification figure written without refitting M_both_526.")
