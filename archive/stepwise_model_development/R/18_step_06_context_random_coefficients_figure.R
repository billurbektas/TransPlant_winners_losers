# Step 6: visualize the species random coefficients from M_context -----------
#
# This script reads the already fitted model and does not fit a new model.
# Response: Delta AUC = warmed AUC minus origin AUC, maximum-duration horizon.

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
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

model_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "auc_change_species_context_model.rds"
)
coefficient_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "context_random_coefficient_scatter_data.csv"
)
figure_stem <- file.path(
  "output", "frequentist_hierarchical_auc", "figures",
  "step_06_context_random_coefficient_map"
)

if (!file.exists(model_file)) {
  stop("M_context was not found. Run R/16_step_04_species_context_model.R first.",
       call. = FALSE)
}
dir.create(dirname(coefficient_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_stem), recursive = TRUE, showWarnings = FALSE)

M_context <- readRDS(model_file)

# Reconstruct the conditional species random coefficients (BLUPs) directly
# from the fitted metafor object: D Z' M^-1 (y - X beta).
model_residual <- as.numeric(
  M_context$yi - M_context$X %*% M_context$beta
)
weighted_residual <- Matrix::solve(M_context$M, model_residual)
species_levels <- levels(M_context$mf.g$outer)
n_species <- length(species_levels)

random_intercept <- as.numeric(
  M_context$sigma2[[1]] *
    Matrix::crossprod(M_context$Z.S[[1]], weighted_residual)
)

G_context <- as.matrix(M_context$G)
if (!identical(colnames(G_context), c("T_z", "P_z", "TP_z"))) {
  stop(
    "Unexpected M_context random-slope columns: ",
    paste(colnames(G_context), collapse = ", "),
    call. = FALSE
  )
}
Z_context <- do.call(cbind, lapply(seq_len(ncol(G_context)), function(j) {
  Matrix::Diagonal(x = as.numeric(M_context$Z.G1[, j])) %*%
    M_context$Z.G2
}))
D_context <- Matrix::kronecker(
  G_context, Matrix::Diagonal(n_species)
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

# Rank is calculated from each species' observed [1, T, P, T*P] design.
scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(species_id = as.character(.data$species_id))

support_rank <- auc_data |>
  dplyr::group_split(.data$species_id, .keep = TRUE) |>
  lapply(function(species) {
    contexts <- dplyr::distinct(species, .data$T_z, .data$P_z)
    rank <- qr(cbind(
      intercept = 1,
      T_z = contexts$T_z,
      P_z = contexts$P_z,
      TP_z = contexts$T_z * contexts$P_z
    ))$rank
    tibble::tibble(
      species_id = species$species_id[[1]],
      n_experiments = dplyr::n_distinct(species$experiment_id),
      context_design_rank = rank,
      rank_support = if (rank == 4L) "Full rank (4)" else "Partial rank (1-3)"
    )
  }) |>
  dplyr::bind_rows()

plot_data <- tibble::tibble(
  species_id = species_levels,
  random_intercept = random_intercept,
  random_T = random_context_matrix[, "T_z"],
  random_P = random_context_matrix[, "P_z"],
  random_TP = random_context_matrix[, "TP_z"]
) |>
  dplyr::left_join(support_rank, by = "species_id") |>
  dplyr::mutate(
    rank_support = factor(
      .data$rank_support,
      levels = c("Partial rank (1-3)", "Full rank (4)")
    )
  )

if (anyNA(plot_data$context_design_rank)) {
  stop("Some fitted species could not be matched to their T-P support rank.",
       call. = FALSE)
}

coefficient_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$random_intercept,
    y = .data$random_TP,
    fill = .data$random_T,
    colour = .data$random_P,
    shape = .data$rank_support
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.35) +
  ggplot2::geom_point(size = 3.2, stroke = 1.05, alpha = 0.9) +
  ggplot2::scale_shape_manual(
    values = c("Partial rank (1-3)" = 21, "Full rank (4)" = 24),
    drop = FALSE
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
    midpoint = 0
  ) +
  ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "#4D4D4D", high = "#B2182B",
    midpoint = 0
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Species random coefficients from M_context",
    subtitle = paste0(
      "Each point is one species. Coefficients are conditional random ",
      "deviations (BLUPs), not total fixed + random effects."
    ),
    x = "Species random intercept deviation (Delta AUC)",
    y = "Species random T x P coefficient",
    fill = "Random T\ncoefficient",
    colour = "Random P\ncoefficient",
    shape = "Species T-P\ndesign support",
    caption = paste0(
      "Full rank means rank 4 for [1, T, P, T x P] within that species; ",
      "partial rank combines ranks 1-3."
    )
  )

invisible(write_csv_atomic(plot_data, coefficient_file))
ggplot2::ggsave(
  paste0(figure_stem, ".pdf"), coefficient_plot,
  width = 11, height = 7.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(figure_stem, ".png"), coefficient_plot,
  width = 11, height = 7.5, dpi = 220, bg = "white"
)

cat("M_context was loaded without refitting.\n")
cat("Species plotted: ", nrow(plot_data), "\n", sep = "")
cat("Full rank: ", sum(plot_data$context_design_rank == 4L), "\n", sep = "")
cat("Partial rank: ", sum(plot_data$context_design_rank < 4L), "\n", sep = "")
cat("Figure: ", normalizePath(paste0(figure_stem, ".pdf")), "\n", sep = "")
cat("Data: ", normalizePath(coefficient_file), "\n", sep = "")

# Full-rank version with every species labelled ------------------------------
full_rank_data <- plot_data |>
  dplyr::filter(.data$context_design_rank == 4L)

full_rank_plot <- ggplot2::ggplot(
  full_rank_data,
  ggplot2::aes(
    x = .data$random_intercept,
    y = .data$random_TP,
    fill = .data$random_T,
    colour = .data$random_P
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.35) +
  ggplot2::geom_point(shape = 21, size = 3.2, stroke = 1.05, alpha = 0.9) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = .data$species_id),
    size = 2.45, colour = "grey15", min.segment.length = 0,
    segment.colour = "grey68", segment.size = 0.25,
    box.padding = 0.28, point.padding = 0.18, max.overlaps = Inf,
    force = 4, force_pull = 0.01, max.time = 12, max.iter = 1000000,
    seed = 20260907
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
    midpoint = 0
  ) +
  ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "#4D4D4D", high = "#B2182B",
    midpoint = 0
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.20)) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 35, 10, 35),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Full-rank species random coefficients from M_context",
    subtitle = paste0(
      "Only species with rank 4 for [1, T, P, T x P]; every species is ",
      "labelled. Values are conditional random deviations (BLUPs)."
    ),
    x = "Species random intercept deviation (Delta AUC)",
    y = "Species random T x P coefficient",
    fill = "Random T\ncoefficient",
    colour = "Random P\ncoefficient"
  )

full_rank_stem <- file.path(
  dirname(figure_stem), "step_06_context_random_coefficient_map_full_rank"
)
ggplot2::ggsave(
  paste0(full_rank_stem, ".pdf"), full_rank_plot,
  width = 20, height = 15, device = grDevices::cairo_pdf, limitsize = FALSE
)
ggplot2::ggsave(
  paste0(full_rank_stem, ".png"), full_rank_plot,
  width = 20, height = 15, dpi = 170, bg = "white", limitsize = FALSE
)

# Select the extrema of the two plotted axes. A species can occupy more than
# one extreme, in which case all selection reasons are retained in its label.
extreme_rows <- dplyr::bind_rows(
  full_rank_data |>
    dplyr::arrange(.data$random_intercept, .data$species_id) |>
    dplyr::slice(1L) |>
    dplyr::mutate(extreme = "minimum random intercept"),
  full_rank_data |>
    dplyr::arrange(dplyr::desc(.data$random_intercept), .data$species_id) |>
    dplyr::slice(1L) |>
    dplyr::mutate(extreme = "maximum random intercept"),
  full_rank_data |>
    dplyr::arrange(.data$random_TP, .data$species_id) |>
    dplyr::slice(1L) |>
    dplyr::mutate(extreme = "minimum random T x P"),
  full_rank_data |>
    dplyr::arrange(dplyr::desc(.data$random_TP), .data$species_id) |>
    dplyr::slice(1L) |>
    dplyr::mutate(extreme = "maximum random T x P")
)
extreme_species <- extreme_rows |>
  dplyr::summarise(
    extreme = paste(.data$extreme, collapse = "; "),
    dplyr::across(
      c("random_intercept", "random_T", "random_P", "random_TP",
        "n_experiments", "context_design_rank"),
      dplyr::first
    ),
    .by = "species_id"
  )

experiment_warming <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
W_z_display <- stats::median(experiment_warming$W_z)
W_display <- stats::median(experiment_warming$W)
fixed_beta <- stats::coef(M_context)

scaling_value <- function(field, statistic) {
  scaled_auc$scaling[[statistic]][
    scaled_auc$scaling$scaled_field == field
  ][[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")

extreme_surface_one <- function(species_name, points = 65L) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == species_name) |>
    dplyr::distinct(
      .data$experiment_id, .data$T_z, .data$P_z, .data$T, .data$P_raw
    )
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  )
  hull <- geometry::convhulln(
    as.matrix(observed[, c("T_z", "P_z")]), options = "Qt"
  )
  candidate <- candidate |>
    dplyr::mutate(
      tile_width = diff(range(observed$T_z)) / (points - 1L) * 1.02,
      tile_height = diff(range(observed$P_z)) / (points - 1L) * 1.02,
      supported = geometry::inhulln(
        hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
      )
    ) |>
    dplyr::filter(.data$supported) |>
    dplyr::mutate(W_z = W_z_display, TP_z = .data$T_z * .data$P_z)

  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, candidate, names(fixed_beta)
  )
  species_row <- extreme_species |>
    dplyr::filter(.data$species_id == species_name)
  candidate |>
    dplyr::mutate(
      estimate = as.numeric(fixed_design %*% fixed_beta) +
        species_row$random_intercept[[1]] +
        species_row$random_T[[1]] * .data$T_z +
        species_row$random_P[[1]] * .data$P_z +
        species_row$random_TP[[1]] * .data$TP_z,
      species_id = species_name,
      n_experiments = species_row$n_experiments[[1]],
      extreme = species_row$extreme[[1]],
      panel = paste0(
        species_name, " (", .data$n_experiments, " experiments)\n",
        .data$extreme
      )
    )
}

extreme_surfaces <- dplyr::bind_rows(lapply(
  extreme_species$species_id, extreme_surface_one
))
panel_levels <- vapply(extreme_species$species_id, function(species_name) {
  row <- extreme_species[extreme_species$species_id == species_name, ]
  paste0(species_name, " (", row$n_experiments, " experiments)\n", row$extreme)
}, character(1))
extreme_surfaces <- extreme_surfaces |>
  dplyr::mutate(panel = factor(.data$panel, levels = panel_levels))

extreme_observations <- auc_data |>
  dplyr::filter(.data$species_id %in% extreme_species$species_id) |>
  dplyr::mutate(
    panel = factor(
      vapply(.data$species_id, function(species_name) {
        row <- extreme_species[extreme_species$species_id == species_name, ]
        paste0(
          species_name, " (", row$n_experiments, " experiments)\n", row$extreme
        )
      }, character(1)),
      levels = panel_levels
    )
  )

surface_limit <- max(abs(extreme_surfaces$estimate), na.rm = TRUE)
extreme_surface_plot <- ggplot2::ggplot(
  extreme_surfaces,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_tile(
    ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
  ) +
  ggplot2::geom_point(
    data = extreme_observations,
    ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
    shape = 21, fill = "white", colour = "grey20", size = 1.6, stroke = 0.4
  ) +
  ggplot2::facet_wrap(~ panel, ncol = 2, scales = "free") +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
    limits = c(-surface_limit, surface_limit)
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
    strip.text = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::labs(
    title = "T-P surfaces for full-rank coefficient extremes",
    subtitle = paste0(
      "Predicted Delta AUC at median warming (Delta T = ",
      sprintf("%.2f", W_display), " °C); predictions are restricted to each ",
      "species' observed T-P convex hull."
    ),
    x = "Destination JJA temperature (°C)",
    y = "Destination JJA precipitation",
    fill = "Predicted Delta AUC"
  )

extreme_surface_stem <- file.path(
  dirname(figure_stem), "step_06_context_extreme_species_tp_surfaces"
)
ggplot2::ggsave(
  paste0(extreme_surface_stem, ".pdf"), extreme_surface_plot,
  width = 12, height = 9, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(extreme_surface_stem, ".png"), extreme_surface_plot,
  width = 12, height = 9, dpi = 180, bg = "white"
)

combined_plot <- patchwork::wrap_plots(
  full_rank_plot, extreme_surface_plot, ncol = 1, heights = c(1.25, 1)
) +
  patchwork::plot_annotation(
    title = "Full-rank M_context coefficients and their extreme T-P surfaces"
  )
combined_stem <- file.path(
  dirname(figure_stem), "step_06_context_full_rank_extremes_with_surfaces"
)
ggplot2::ggsave(
  paste0(combined_stem, ".pdf"), combined_plot,
  width = 17, height = 22, device = grDevices::cairo_pdf, limitsize = FALSE
)
ggplot2::ggsave(
  paste0(combined_stem, ".png"), combined_plot,
  width = 17, height = 22, dpi = 150, bg = "white", limitsize = FALSE
)

extreme_species_file <- file.path(
  dirname(coefficient_file), "context_full_rank_extreme_species.csv"
)
extreme_surface_file <- file.path(
  dirname(coefficient_file), "context_full_rank_extreme_species_surfaces.csv"
)
invisible(write_csv_atomic(extreme_species, extreme_species_file))
invisible(write_csv_atomic(extreme_surfaces, extreme_surface_file))

cat("Full-rank labelled figure: ",
    normalizePath(paste0(full_rank_stem, ".pdf")), "\n", sep = "")
cat("Extreme-species surfaces: ",
    normalizePath(paste0(extreme_surface_stem, ".pdf")), "\n", sep = "")
cat("Combined figure: ",
    normalizePath(paste0(combined_stem, ".pdf")), "\n", sep = "")
cat("Selected extreme species:\n")
print(dplyr::select(
  extreme_species, "species_id", "extreme", "n_experiments",
  "random_intercept", "random_TP"
))
