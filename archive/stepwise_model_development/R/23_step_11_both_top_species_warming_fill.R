# Step 11: M_both top-two species per fitted class ---------------------------
#
# Loads the fitted M_both and saved classification. No model is fitted here.
# Circle fill is the conditional species random warming-slope deviation.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "ggrepel", "metafor", "patchwork",
  "readr", "tibble", "tidyr"
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
  stop("Run R/21_step_09_both_model.R and R/22_step_10_both_classification_figures.R first.",
       call. = FALSE)
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
    ),
    response_span = .data$max_fitted_delta_auc - .data$min_fitted_delta_auc
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
if (nrow(classification) != 118L || anyNA(classification$random_W)) {
  stop("The M_both classification/random-slope join is incomplete.",
       call. = FALSE)
}

# "Best" is defined on the fitted biological response scale, without using
# uncertainty or p-values: highest minimum for winners, lowest maximum for
# losers, and widest supported fitted range for context-dependent species.
top_winners <- classification |>
  dplyr::filter(.data$fitted_class == "Consistent winner") |>
  dplyr::arrange(
    dplyr::desc(.data$min_fitted_delta_auc),
    dplyr::desc(.data$max_fitted_delta_auc), .data$species_id
  ) |>
  dplyr::slice_head(n = 2L) |>
  dplyr::mutate(selection_group = "Top consistent winner")
top_losers <- classification |>
  dplyr::filter(.data$fitted_class == "Consistent loser") |>
  dplyr::arrange(
    .data$max_fitted_delta_auc, .data$min_fitted_delta_auc, .data$species_id
  ) |>
  dplyr::slice_head(n = 2L) |>
  dplyr::mutate(selection_group = "Top consistent loser")
top_context <- classification |>
  dplyr::filter(.data$fitted_class == "Context-dependent") |>
  dplyr::arrange(dplyr::desc(.data$response_span), .data$species_id) |>
  dplyr::slice_head(n = 2L) |>
  dplyr::mutate(selection_group = "Top context-dependent")
selection <- dplyr::bind_rows(top_winners, top_losers, top_context) |>
  dplyr::mutate(
    display_order = dplyr::row_number(),
    species_panel = paste0(.data$species_id, "\n(", .data$n_experiments,
                           " experiments)")
  )
selected_species <- selection$species_id
panel_levels <- selection$species_panel
invisible(write_csv_atomic(
  selection,
  file.path(both_dir, "both_top_two_species_per_class.csv")
))

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    TP_z = .data$T_z * .data$P_z
  )
experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
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
fixed_beta <- stats::coef(M_both)

surface_one <- function(species_name, grid_points = 71L) {
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
  candidate <- candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    tidyr::crossing(warming_levels) |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, candidate, names(fixed_beta)
  )
  random <- random_coefficients |>
    dplyr::filter(.data$species_id == species_name)
  metadata <- selection |>
    dplyr::filter(.data$species_id == species_name)
  candidate |>
    dplyr::mutate(
      estimate = as.numeric(fixed_design %*% fixed_beta) +
        random$random_intercept[[1]] + random$random_W[[1]] * .data$W_z +
        random$random_T[[1]] * .data$T_z +
        random$random_P[[1]] * .data$P_z +
        random$random_TP[[1]] * .data$TP_z,
      species_id = species_name,
      species_panel = metadata$species_panel[[1]],
      tile_width = diff(range(observed$T_z)) / (grid_points - 1L) * 1.02,
      tile_height = diff(range(observed$P_z)) / (grid_points - 1L) * 1.02
    )
}
surfaces <- dplyr::bind_rows(lapply(selected_species, surface_one)) |>
  dplyr::mutate(
    species_panel = factor(.data$species_panel, levels = panel_levels),
    warming_panel = factor(
      .data$warming_panel, levels = levels(warming_levels$warming_panel)
    )
  )
observations <- auc_data |>
  dplyr::filter(.data$species_id %in% selected_species) |>
  dplyr::distinct(.data$species_id, .data$experiment_id, .data$T_z, .data$P_z) |>
  dplyr::left_join(
    dplyr::select(selection, "species_id", "species_panel"),
    by = "species_id"
  ) |>
  tidyr::crossing(dplyr::select(warming_levels, "warming_panel")) |>
  dplyr::mutate(
    species_panel = factor(.data$species_panel, levels = panel_levels),
    warming_panel = factor(
      .data$warming_panel, levels = levels(warming_levels$warming_panel)
    )
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
surface_limit <- max(abs(surfaces$estimate), na.rm = TRUE)

surface_plot <- ggplot2::ggplot(
  surfaces,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_tile(
    ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
  ) +
  ggplot2::geom_point(
    data = observations,
    ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
    shape = 21, fill = "white", colour = "grey20", size = 0.9, stroke = 0.3
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$warming_panel),
    cols = ggplot2::vars(.data$species_panel),
    scales = "fixed", space = "fixed"
  ) +
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
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(face = "bold", size = 7),
    strip.text.y = ggplot2::element_text(face = "bold", size = 7),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "M_both: two strongest species per fitted class",
    subtitle = paste0(
      "Fixed T-P axes and shared fitted-Delta-AUC scale; predictions stay ",
      "inside species-specific observed T-P hulls."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Fitted Delta AUC"
  )

# Classification panel: class remains the outline; random W BLUP is the fill.
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
class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#5E3C99",
  "Consistent winner" = "#1B9E77"
)
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
      fill = .data$random_W
    ), shape = 21, size = 3, stroke = 0.9, alpha = 0.96
  ) +
  ggrepel::geom_text_repel(
    data = selection,
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
    panel.grid.minor = ggplot2::element_blank(), legend.position = "top",
    plot.margin = ggplot2::margin(6, 15, 6, 6)
  ) +
  ggplot2::labs(
    title = "M_both classification with random warming-slope fill",
    subtitle = "Only the six surface species are labelled.",
    x = "Minimum fitted Delta AUC",
    y = "Maximum fitted Delta AUC",
    colour = "Fitted classification",
    caption = paste0(
      "Pointwise 95% intervals at fitted extrema. Circle fill is the ",
      "conditional random W-slope deviation, not the total warming slope."
    )
  )

classification_stem <- file.path(
  figure_dir, "step_10_both_classification_top2_random_warming_fill"
)
surface_stem <- file.path(
  figure_dir, "step_10_both_top2_species_surfaces_median_high"
)
combined_stem <- file.path(
  figure_dir,
  "step_10_both_classification_random_warming_fill_top2_surfaces_median_high"
)
combined_plot <- patchwork::wrap_plots(
  classification_plot, surface_plot,
  ncol = 2, widths = c(0.95, 2.05)
) +
  patchwork::plot_annotation(
    title = "M_both classification and top-two species per fitted class"
  )
ggplot2::ggsave(
  paste0(classification_stem, ".pdf"), classification_plot,
  width = 11, height = 9, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(classification_stem, ".png"), classification_plot,
  width = 11, height = 9, dpi = 220, bg = "white"
)
ggplot2::ggsave(
  paste0(surface_stem, ".pdf"), surface_plot,
  width = 20, height = 8.5, device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggplot2::ggsave(
  paste0(surface_stem, ".png"), surface_plot,
  width = 20, height = 8.5, dpi = 170, bg = "white", limitsize = FALSE
)
ggplot2::ggsave(
  paste0(combined_stem, ".pdf"), combined_plot,
  width = 29, height = 10.5, device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggplot2::ggsave(
  paste0(combined_stem, ".png"), combined_plot,
  width = 29, height = 10.5, dpi = 160, bg = "white", limitsize = FALSE
)

cat("M_both was loaded without refitting.\n")
cat("Selected species:\n")
print(dplyr::select(
  selection, "selection_group", "species_id", "n_experiments",
  "min_fitted_delta_auc", "max_fitted_delta_auc", "random_W"
))
cat("Combined figure: ", normalizePath(paste0(combined_stem, ".pdf")),
    "\n", sep = "")
