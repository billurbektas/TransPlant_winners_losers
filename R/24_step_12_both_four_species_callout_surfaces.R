# Step 12: M_both classification with four linked surface callouts ----------
#
# Loads the fitted M_both object. No model is fitted here. The classification
# plot is central and each requested species is linked to a species-specific
# fitted T-P surface restricted to its own observed convex hull.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "cowplot", "dplyr", "geometry", "ggplot2", "metafor", "readr",
  "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

both_dir <- file.path("output", "step_by_step_modeling", "species_both")
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
model_file <- file.path(both_dir, "auc_change_species_both_model.rds")
classification_file <- file.path(
  both_dir, "both_full_rank_minmax_classification.csv"
)
random_file <- file.path(both_dir, "species_random_both_coefficients.csv")
required_files <- c(model_file, classification_file, random_file)
if (any(!file.exists(required_files))) {
  stop("Run Steps 9 and 10 before this plotting script.", call. = FALSE)
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

M_both <- readRDS(model_file)
classification <- readr::read_csv(
  classification_file, show_col_types = FALSE
) |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    fitted_class = factor(
      .data$fitted_class,
      levels = c("Consistent loser", "Context-dependent", "Consistent winner")
    )
  )
random_coefficients <- readr::read_csv(
  random_file, show_col_types = FALSE
) |>
  dplyr::mutate(species_id = as.character(.data$species_id))
classification <- classification |>
  dplyr::left_join(
    dplyr::select(random_coefficients, "species_id", "random_W"),
    by = "species_id"
  )

requested_species <- c(
  "Potentilla erecta", "Plantago lanceolata", "Elymus repens",
  "Silene acaulis"
)
selection <- classification |>
  dplyr::filter(.data$species_id %in% requested_species) |>
  dplyr::mutate(display_order = match(.data$species_id, requested_species)) |>
  dplyr::arrange(.data$display_order)
if (!identical(selection$species_id, requested_species)) {
  stop("One or more requested species are absent from the 118-species set.",
       call. = FALSE)
}

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    TP_z = .data$T_z * .data$P_z
  )
fixed_beta <- stats::coef(M_both)

scaling_value <- function(field, statistic) {
  scaled_auc$scaling[[statistic]][
    scaled_auc$scaling$scaled_field == field
  ][[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
raw_T_labels <- function(x) sprintf("%.1f", x * T_scale + T_center)
raw_P_labels <- function(x) sprintf("%.0f", exp(x * P_scale + P_center))

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
warming_lookup <- tibble::tibble(
  warming_name = c("median", "high"),
  probability = c(0.50, 0.75),
  W_z = as.numeric(stats::quantile(
    experiment_warming$W_z, c(0.50, 0.75)
  )),
  W_C = as.numeric(stats::quantile(
    experiment_warming$W, c(0.50, 0.75)
  ))
)

surface_one <- function(species_name, W_z_value, grid_points = 81L) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == species_name) |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = grid_points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = grid_points)
  )
  hull <- geometry::convhulln(
    as.matrix(observed[, c("T_z", "P_z")]), options = "Qt"
  )
  candidate <- candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::mutate(
      W_z = W_z_value,
      TP_z = .data$T_z * .data$P_z
    )
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, candidate, names(fixed_beta)
  )
  random <- random_coefficients |>
    dplyr::filter(.data$species_id == species_name)
  candidate |>
    dplyr::mutate(
      estimate = as.numeric(fixed_design %*% fixed_beta) +
        random$random_intercept[[1]] + random$random_W[[1]] * .data$W_z +
        random$random_T[[1]] * .data$T_z +
        random$random_P[[1]] * .data$P_z +
        random$random_TP[[1]] * .data$TP_z,
      species_id = species_name,
      tile_width = diff(range(observed$T_z)) / (grid_points - 1L) * 1.025,
      tile_height = diff(range(observed$P_z)) / (grid_points - 1L) * 1.025
    ) |>
    list(grid = _, observed = observed)
}

surface_sets <- lapply(seq_len(nrow(warming_lookup)), function(index) {
  stats::setNames(
    lapply(requested_species, surface_one,
           W_z_value = warming_lookup$W_z[[index]]),
    requested_species
  )
})
names(surface_sets) <- warming_lookup$warming_name
all_surface_estimates <- unlist(lapply(
  surface_sets,
  function(level) unlist(lapply(level, function(item) item$grid$estimate))
))
surface_limit <- max(abs(all_surface_estimates), na.rm = TRUE)

class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#5E3C99",
  "Consistent winner" = "#1B9E77"
)
species_colours <- c(
  "Potentilla erecta" = "#1B9E77",
  "Plantago lanceolata" = "#D95F02",
  "Elymus repens" = "#7570B3",
  "Silene acaulis" = "#E7298A"
)

x_range <- range(c(
  classification$min_conf_low, classification$min_conf_high,
  classification$min_fitted_delta_auc, 0
), na.rm = TRUE)
y_range <- range(c(
  classification$max_conf_low, classification$max_conf_high,
  classification$max_fitted_delta_auc, 0
), na.rm = TRUE)
x_limits <- x_range + c(-0.08, 0.08) * diff(x_range)
y_limits <- y_range + c(-0.08, 0.08) * diff(y_range)
warming_limit <- max(abs(classification$random_W), na.rm = TRUE)

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
    colour = "grey55", linewidth = 0.4
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.5) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_conf_low, xend = .data$min_conf_high,
      y = .data$max_fitted_delta_auc, yend = .data$max_fitted_delta_auc,
      colour = .data$fitted_class
    ), linewidth = 0.35, alpha = 0.17
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_fitted_delta_auc, xend = .data$min_fitted_delta_auc,
      y = .data$max_conf_low, yend = .data$max_conf_high,
      colour = .data$fitted_class
    ), linewidth = 0.35, alpha = 0.17
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      colour = .data$fitted_class,
      fill = .data$random_W
    ), shape = 21, size = 2.5, stroke = 0.75, alpha = 0.94
  ) +
  ggplot2::geom_point(
    data = selection,
    ggplot2::aes(
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      fill = .data$random_W
    ), shape = 21, size = 4.2, stroke = 1.25,
    colour = unname(species_colours[selection$species_id]),
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = class_colours, drop = FALSE) +
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
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  ) +
  ggplot2::labs(
    x = "Minimum fitted Delta AUC",
    y = "Maximum fitted Delta AUC",
    colour = "Fitted classification"
  )

surface_plot_one <- function(surface_item, species_name, W_C_value) {
  ggplot2::ggplot(
    surface_item$grid,
    ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
    ) +
    ggplot2::geom_point(
      data = surface_item$observed,
      ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
      shape = 21, fill = "white", colour = "grey15", size = 1.4,
      stroke = 0.35
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
      limits = c(-surface_limit, surface_limit), guide = "none"
    ) +
    ggplot2::scale_x_continuous(labels = raw_T_labels) +
    ggplot2::scale_y_continuous(labels = raw_P_labels) +
    ggplot2::theme_bw(base_size = 7.4) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold.italic", size = 9,
        colour = species_colours[[species_name]], hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0.5),
      plot.background = ggplot2::element_rect(
        colour = species_colours[[species_name]], linewidth = 1.1,
        fill = "white"
      ),
      plot.margin = ggplot2::margin(5, 5, 4, 5),
      axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 6.5)
    ) +
    ggplot2::labs(
      title = species_name,
      subtitle = paste0(
        dplyr::filter(selection, .data$species_id == species_name)$n_experiments,
        " experiments; Delta T = ", sprintf("%.2f", W_C_value), " °C"
      ),
      x = "Destination temperature (°C)",
      y = "Destination precipitation"
    )
}

# Canvas positions: two callouts on each side of the central classification.
class_box <- c(x = 0.265, y = 0.155, width = 0.47, height = 0.69)
inset_boxes <- list(
  "Elymus repens" = c(x = 0.012, y = 0.535, width = 0.225, height = 0.39),
  "Silene acaulis" = c(x = 0.012, y = 0.075, width = 0.225, height = 0.39),
  "Potentilla erecta" = c(x = 0.763, y = 0.535, width = 0.225, height = 0.39),
  "Plantago lanceolata" = c(x = 0.763, y = 0.075, width = 0.225, height = 0.39)
)

# Approximate conversion from classification data coordinates to page
# coordinates. The fractions account for axes and the bottom legends within
# the central ggplot box; the subsequent image check guards the layout.
point_on_canvas <- function(species_name) {
  row <- dplyr::filter(selection, .data$species_id == species_name)
  x_fraction <- (row$min_fitted_delta_auc - x_limits[[1]]) /
    diff(x_limits)
  y_fraction <- (row$max_fitted_delta_auc - y_limits[[1]]) /
    diff(y_limits)
  point <- c(
    x = class_box[["x"]] + class_box[["width"]] *
      (0.115 + 0.855 * x_fraction),
    y = class_box[["y"]] + class_box[["height"]] *
      (0.205 + 0.745 * y_fraction)
  )
  # Small species-specific registration offsets align connector endpoints with
  # the highlighted circle centers after ggplot's axes/legend are laid out.
  registration <- list(
    "Elymus repens" = c(x = -0.010, y = -0.024),
    "Silene acaulis" = c(x = -0.018, y = -0.007),
    "Potentilla erecta" = c(x = 0, y = 0),
    "Plantago lanceolata" = c(x = -0.009, y = 0.006)
  )
  point + registration[[species_name]]
}

connector_target <- function(species_name) {
  box <- inset_boxes[[species_name]]
  if (box[["x"]] < class_box[["x"]]) {
    c(x = box[["x"]] + box[["width"]],
      y = box[["y"]] + box[["height"]] / 2)
  } else {
    c(x = box[["x"]], y = box[["y"]] + box[["height"]] / 2)
  }
}

build_callout <- function(warming_name) {
  warming_row <- dplyr::filter(
    warming_lookup, .data$warming_name == .env$warming_name
  )
  mini_plots <- stats::setNames(lapply(requested_species, function(species_name) {
    surface_plot_one(
      surface_sets[[warming_name]][[species_name]], species_name,
      warming_row$W_C[[1]]
    )
  }), requested_species)

  canvas <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "M_both species classification with linked T-P response surfaces — ",
        tools::toTitleCase(warming_name), " warming"
      ),
      x = 0.5, y = 0.987, hjust = 0.5, vjust = 1,
      fontface = "bold", size = 15
    ) +
    cowplot::draw_label(
      if (warming_name == "median") {
        paste0(
          "Central classification and surface fill use Delta T = ",
          sprintf("%.2f", warming_row$W_C[[1]]),
          " °C; every surface stays inside that species' observed T-P hull."
        )
      } else {
        paste0(
          "Central positions retain the median-warming classification; surface ",
          "fill uses Delta T = ", sprintf("%.2f", warming_row$W_C[[1]]),
          " °C within each species' observed T-P hull."
        )
      },
      x = 0.5, y = 0.958, hjust = 0.5, vjust = 1, size = 9
    ) +
    cowplot::draw_plot(
      classification_plot,
      x = class_box[["x"]], y = class_box[["y"]],
      width = class_box[["width"]], height = class_box[["height"]]
    )

  # Draw the connectors above the central plot and below the callout panels.
  for (species_name in names(inset_boxes)) {
    start <- point_on_canvas(species_name)
    finish <- connector_target(species_name)
    canvas <- canvas + cowplot::draw_line(
      x = c(start[["x"]], finish[["x"]]),
      y = c(start[["y"]], finish[["y"]]),
      colour = species_colours[[species_name]], linewidth = 0.8,
      alpha = 0.9
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
  canvas + cowplot::draw_label(
    paste0(
      "Circle fill in the central panel is the species random warming-slope ",
      "deviation; surface fill is fitted Delta AUC."
    ),
    x = 0.5, y = 0.018, hjust = 0.5, vjust = 0, size = 8
  )
}

median_callout <- build_callout("median")
high_callout <- build_callout("high")
median_stem <- file.path(
  figure_dir, "step_11_both_four_species_linked_callouts_median"
)
high_stem <- file.path(
  figure_dir, "step_11_both_four_species_linked_callouts_high"
)
ggplot2::ggsave(
  paste0(median_stem, ".pdf"), median_callout,
  width = 16, height = 12, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(median_stem, ".png"), median_callout,
  width = 16, height = 12, dpi = 190, bg = "white"
)
ggplot2::ggsave(
  paste0(high_stem, ".pdf"), high_callout,
  width = 16, height = 12, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(high_stem, ".png"), high_callout,
  width = 16, height = 12, dpi = 190, bg = "white"
)

cat("M_both was loaded without refitting.\n")
cat("Linked callout species:\n")
print(dplyr::select(
  selection, "species_id", "fitted_class", "n_experiments", "random_W"
))
cat("Median-warming callout: ", normalizePath(paste0(median_stem, ".pdf")),
    "\n", sep = "")
cat("High-warming callout: ", normalizePath(paste0(high_stem, ".pdf")),
    "\n", sep = "")
