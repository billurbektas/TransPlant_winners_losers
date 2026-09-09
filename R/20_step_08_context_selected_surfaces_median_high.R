# Step 8: selected full-rank species at median and high warming ---------------
#
# This script reads the already fitted M_context model and the Step-7 min/max
# table. It does not fit a new model. All T-P surface facets use identical x,
# y, and fitted-Delta-AUC scales.

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

model_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "auc_change_species_context_model.rds"
)
random_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "species_random_context_coefficients.csv"
)
classification_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "context_full_rank_minmax_classification.csv"
)
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
surface_stem <- file.path(
  figure_dir, "step_08_context_selected_species_tp_surfaces_median_high"
)
combined_stem <- file.path(
  figure_dir, "step_08_context_classification_selected_surfaces_median_high"
)
selection_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "context_step08_selected_species.csv"
)

required_files <- c(model_file, random_file, classification_file)
if (any(!file.exists(required_files))) {
  stop(
    "Required saved Step-4/Step-7 outputs are missing: ",
    paste(required_files[!file.exists(required_files)], collapse = ", "),
    call. = FALSE
  )
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

M_context <- readRDS(model_file)
fixed_beta <- stats::coef(M_context)
random_coefficients <- readr::read_csv(random_file, show_col_types = FALSE) |>
  dplyr::mutate(species_id = as.character(.data$species_id))
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

# Representative winner: strongest positive lower bound, then closest to 1:1.
# Representative loser: strongest negative upper bound, then closest to 1:1.
representative_winner <- classification |>
  dplyr::filter(.data$fitted_class == "Consistent winner") |>
  dplyr::arrange(
    dplyr::desc(.data$min_fitted_delta_auc), .data$response_span,
    .data$species_id
  ) |>
  dplyr::slice(1L) |>
  dplyr::pull(.data$species_id)
representative_loser <- classification |>
  dplyr::filter(.data$fitted_class == "Consistent loser") |>
  dplyr::arrange(
    .data$max_fitted_delta_auc, .data$response_span, .data$species_id
  ) |>
  dplyr::slice(1L) |>
  dplyr::pull(.data$species_id)

selected_species <- unique(c(
  representative_winner,
  representative_loser,
  "Nardus stricta",
  "Silene acaulis",
  "Poa alpina",
  "Achillea millefolium",
  "Agrostis capillaris"
))
missing_species <- setdiff(selected_species, classification$species_id)
if (length(missing_species)) {
  stop("Selected species are absent from the eligible full-rank set: ",
       paste(missing_species, collapse = ", "), call. = FALSE)
}

selection <- classification |>
  dplyr::filter(.data$species_id %in% selected_species) |>
  dplyr::mutate(
    display_order = match(.data$species_id, selected_species),
    selection_reason = dplyr::case_when(
      .data$species_id == representative_winner ~
        "representative consistent winner: highest fitted minimum",
      .data$species_id == representative_loser ~
        "representative consistent loser: lowest fitted maximum",
      TRUE ~ "requested context-dependent species"
    ),
    species_panel = paste0(
      .data$species_id, "\n(", .data$n_experiments, " experiments)"
    )
  ) |>
  dplyr::arrange(.data$display_order)
panel_levels <- selection$species_panel
invisible(write_csv_atomic(selection, selection_file))

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(species_id = as.character(.data$species_id))

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
warming_z <- stats::quantile(
  experiment_warming$W_z, c(0.50, 0.75), names = FALSE
)
warming_C <- stats::quantile(
  experiment_warming$W, c(0.50, 0.75), names = FALSE
)
warming_levels <- tibble::tibble(
  warming = factor(
    c("Median warming", "High warming"),
    levels = c("Median warming", "High warming")
  ),
  W_z = warming_z,
  warming_C = warming_C,
  warming_panel = factor(
    paste0(c("Median warming", "High warming"), "\nDelta T = ",
           sprintf("%.2f", warming_C), " °C"),
    levels = paste0(c("Median warming", "High warming"), "\nDelta T = ",
                    sprintf("%.2f", warming_C), " °C")
  )
)

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
        random$random_intercept[[1]] + random$random_T[[1]] * .data$T_z +
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
  tidyr::crossing(
    dplyr::select(warming_levels, "warming_panel")
  ) |>
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
    strip.text.x = ggplot2::element_text(face = "bold", size = 6.8),
    strip.text.y = ggplot2::element_text(face = "bold", size = 7),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Selected full-rank species across median and high warming",
    subtitle = paste0(
      "Every facet has identical T, P, and color scales; surfaces are ",
      "restricted to species-specific observed T-P convex hulls."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Fitted Delta AUC"
  )

# Recreate the Step-7 classification panel, labelling only displayed species.
classification_valid <- with(
  classification,
  (fitted_class == "Consistent winner" & min_fitted_delta_auc > 0) |
    (fitted_class == "Consistent loser" & max_fitted_delta_auc < 0) |
    (fitted_class == "Context-dependent" &
       min_fitted_delta_auc <= 0 & max_fitted_delta_auc >= 0)
)
if (!all(classification_valid)) {
  stop("The saved fitted classifications do not obey the min/max rules.",
       call. = FALSE)
}
x_data_range <- range(c(
  classification$min_conf_low, classification$min_conf_high,
  classification$min_fitted_delta_auc, 0
), na.rm = TRUE)
y_data_range <- range(c(
  classification$max_conf_low, classification$max_conf_high,
  classification$max_fitted_delta_auc, 0
), na.rm = TRUE)
x_limits <- x_data_range + c(-0.08, 0.08) * diff(x_data_range)
y_limits <- y_data_range + c(-0.08, 0.08) * diff(y_data_range)
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
    ),
    shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = selection,
    ggplot2::aes(
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      label = .data$species_id,
      colour = .data$fitted_class
    ),
    size = 2.8, show.legend = FALSE, min.segment.length = 0,
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
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    plot.margin = ggplot2::margin(6, 15, 6, 6)
  ) +
  ggplot2::labs(
    title = "Minimum versus maximum fitted Delta AUC",
    subtitle = "Only the seven species displayed at right are labelled.",
    x = "Minimum fitted Delta AUC",
    y = "Maximum fitted Delta AUC",
    colour = "Fitted classification"
  )

combined_plot <- patchwork::wrap_plots(
  classification_plot, surface_plot,
  ncol = 2, widths = c(0.92, 2.08)
) +
  patchwork::plot_annotation(
    title = "Selected species classification and T-P response surfaces"
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

cat("M_context was loaded without refitting.\n")
cat("Representative winner: ", representative_winner, "\n", sep = "")
cat("Representative loser: ", representative_loser, "\n", sep = "")
cat("Median warming: ", sprintf("%.2f", warming_C[[1]]), " °C\n", sep = "")
cat("High warming: ", sprintf("%.2f", warming_C[[2]]), " °C\n", sep = "")
cat("All surface panels use fixed T-P scales and one shared color scale.\n")
cat("Surfaces: ", normalizePath(paste0(surface_stem, ".pdf")), "\n", sep = "")
cat("Combined: ", normalizePath(paste0(combined_stem, ".pdf")), "\n", sep = "")
